/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_interpolant.h"
#include "impl.h"

/*
    Divide and conquer on D, the number of columns of V, which are sorted by
    point (see impl.h). A node where every point has the same number n of
    columns goes to pmintbasis, on the m x n blocks of V (windows, no copy);
    a node with few columns goes to the direct Krylov function, with n the
    average number of columns at one point; otherwise the first half is the
    first c = floor(D/2) columns (the columns of a point may fall in both
    halves) and the second half is the residual of the other columns,

        V2[:, b] = P1(pts[b]) V[:, b],

    computed from the coefficients of P1 (_nmod_poly_mat_evaluate_columns)
    when there are few columns per point, and with one evaluation of P1 for
    each point otherwise. The basis is P2 P1, as in pmintbasis.
*/

void
_nmod_poly_mat_interpolant_basis(nmod_poly_mat_t P, slong * shift,
                                 const nmod_mat_t V, nn_srcptr pts,
                                 const _nmod_poly_mat_basis_dac_param_t prm)
{
    const slong m = V->r, D = V->c;
    slong ndist = 0, nmin = WORD_MAX, nmax = 0;
    slong * start;

    if (D == 0)
    {
        nmod_poly_mat_one(P);
        return;
    }

    /* groups of columns with the same point: start[k] .. start[k+1] - 1 */
    start = flint_malloc((D + 1) * sizeof(slong));
    for (slong b = 0; b < D; b++)
        if (b == 0 || pts[b] != pts[b - 1])
            start[ndist++] = b;
    start[ndist] = D;
    for (slong k = 0; k < ndist; k++)
    {
        nmin = FLINT_MIN(nmin, start[k + 1] - start[k]);
        nmax = FLINT_MAX(nmax, start[k + 1] - start[k]);
    }

    /* few columns: direct Krylov call, with ceil(D / ndist) columns per
     * point on average (only if the counts differ: one count goes to
     * pmintbasis and its Krylov leaf) */
    if (nmin < nmax
        && _nmod_poly_mat_krylov_wanted(m, (D + ndist - 1) / ndist, D, prm->direct_cutoff))
    {
        _nmod_poly_mat_interpolant_basis_krylov_leaf(P, shift, V, pts);
        flint_free(start);
        return;
    }

    /* as many columns at every point: pmintbasis, with its Krylov leaf, on
     * windows of V */
    if (nmin == nmax)
    {
        nmod_mat_struct * E = flint_malloc(ndist * sizeof(nmod_mat_struct));
        ulong * dpts = flint_malloc(ndist * sizeof(ulong));

        for (slong k = 0; k < ndist; k++)
        {
            nmod_mat_window_init(E + k, V, 0, start[k], m, start[k + 1]);
            dpts[k] = pts[start[k]];
        }

        _nmod_poly_mat_pmintbasis(P, shift, dpts, E, ndist, prm->leaf_cutoff);

        for (slong k = 0; k < ndist; k++)
            nmod_mat_window_clear(E + k);
        flint_free(dpts);
        flint_free(E);
        flint_free(start);
        return;
    }

    /* split after c = floor(D/2) columns */
    {
        const slong c = D / 2;
        nmod_mat_t V1, V2, W, W2, E;
        nmod_poly_mat_t P1, P2;

        nmod_poly_mat_init(P1, m, m, V->mod.n);
        nmod_mat_window_init(V1, V, 0, 0, m, c);
        _nmod_poly_mat_interpolant_basis(P1, shift, V1, pts, prm);
        nmod_mat_window_clear(V1);

        /* residual: from the coefficients of P1 when the second half has few
         * columns per point (n^2 <= m and 8 n <= m for its average n, as in
         * pmintbasis), otherwise one evaluation of P1 for each of its
         * points */
        nmod_mat_init(V2, m, D - c, V->mod.n);
        {
            slong np2 = 0, n2;
            for (slong k = 0; k < ndist; k++)
                if (start[k] >= c)
                    np2++;
            if (start[ndist - np2] > c)  /* c inside a group */
                np2++;
            n2 = (D - c + np2 - 1) / np2;

            if (n2 * n2 <= m && 8 * n2 <= m)
            {
                nmod_mat_window_init(W, V, 0, c, m, D);
                _nmod_poly_mat_evaluate_columns(V2, P1, W, pts + c);
                nmod_mat_window_clear(W);
            }
            else
            {
                nmod_mat_init(E, m, m, V->mod.n);
                for (slong b = c; b < D; )
                {
                    slong e = b + 1;
                    while (e < D && pts[e] == pts[b])
                        e++;
                    nmod_poly_mat_evaluate_nmod(E, P1, pts[b]);
                    nmod_mat_window_init(W, V, 0, b, m, e);
                    nmod_mat_window_init(W2, V2, 0, b - c, m, e - c);
                    nmod_mat_mul(W2, E, W);
                    nmod_mat_window_clear(W2);
                    nmod_mat_window_clear(W);
                    b = e;
                }
                nmod_mat_clear(E);
            }
        }

        nmod_poly_mat_init(P2, m, m, V->mod.n);
        _nmod_poly_mat_interpolant_basis(P2, shift, V2, pts + c, prm);
        nmod_poly_mat_multiply(P, P2, P1);

        nmod_poly_mat_clear(P2);
        nmod_mat_clear(V2);
        nmod_poly_mat_clear(P1);
    }

    flint_free(start);
}

typedef struct
{
    const nmod_mat_struct * V;
    nn_srcptr pts;
    const _nmod_poly_mat_basis_dac_param_struct * prm;
}
_interpolant_basis_data_struct;

static void
_interpolant_basis_callback(nmod_poly_mat_t P, slong * shift, const void * data)
{
    const _interpolant_basis_data_struct * dt = data;
    _nmod_poly_mat_interpolant_basis(P, shift, dt->V, dt->pts, dt->prm);
}

typedef struct
{
    ulong pt;
    slong col;
}
_point_column_struct;

static int
_point_column_cmp(const void * a, const void * b)
{
    const _point_column_struct * x = a, * y = b;
    if (x->pt != y->pt)
        return (x->pt > y->pt) - (x->pt < y->pt);
    return (x->col > y->col) - (x->col < y->col);
}

void
nmod_poly_mat_interpolant_basis(nmod_poly_mat_t P, slong * rdeg,
                                const nmod_mat_t V, nn_srcptr pts,
                                const slong * shift, poly_mat_form_t form)
{
    const slong m = V->r, D = V->c;
    slong * s;
    int sorted = 1;
    nmod_mat_t Vs;
    nn_ptr spts;
    _nmod_poly_mat_basis_dac_param_t prm;

    if (P->r != m || P->c != m)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_interpolant_basis: P must be V->r x V->r\n");
    if (form > POPOV)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_interpolant_basis: form must be at most POPOV\n");

    if (m == 0)
        return;

    /* columns sorted by point (the module does not depend on their order) */
    for (slong b = 1; b < D && sorted; b++)
        if (pts[b] < pts[b - 1])
            sorted = 0;
    spts = flint_malloc((D + 1) * sizeof(ulong));
    if (sorted)
    {
        nmod_mat_window_init(Vs, V, 0, 0, m, D);
        for (slong b = 0; b < D; b++)
            spts[b] = pts[b];
    }
    else
    {
        _point_column_struct * pc = flint_malloc(D * sizeof(_point_column_struct));
        for (slong b = 0; b < D; b++)
        {
            pc[b].pt = pts[b];
            pc[b].col = b;
        }
        qsort(pc, D, sizeof(_point_column_struct), _point_column_cmp);
        nmod_mat_init(Vs, m, D, V->mod.n);
        for (slong b = 0; b < D; b++)
        {
            spts[b] = pc[b].pt;
            for (slong i = 0; i < m; i++)
                nmod_mat_entry(Vs, i, b) = nmod_mat_entry(V, i, pc[b].col);
        }
        flint_free(pc);
    }

    /* s: shift on input, shifted row degree on output; shift and rdeg are
     * read and written only here, so they may alias */
    s = flint_malloc(m * sizeof(slong));
    for (slong i = 0; i < m; i++)
        s[i] = shift ? shift[i] : 0;

    _nmod_poly_mat_basis_dac_param_init(prm);
    _nmod_poly_mat_interpolant_basis(P, s, Vs, spts, prm);

    if (form == POPOV)
    {
        slong * s0 = flint_malloc(m * sizeof(slong));
        for (slong i = 0; i < m; i++)
            s0[i] = shift ? shift[i] : 0;

        /* TODO this test is useful to avoid useless computations when the base case
         * (Krylov leaves) already provides a Popov form. But a common situation is
         * also when the obtained basis P is already Popov "up to a constant transformation", 
         * would it make sense to add such a test (based on the -delta-leading matrix being
         * invertible), to avoid re-doing a whole basis computation in these cases?
         * (also: check the similar call for approximant bases) */
        if (!nmod_poly_mat_is_popov(P, s0, ROW_LOWER))
        {
            _interpolant_basis_data_struct dt = {Vs, spts, prm};
            _nmod_poly_mat_popov_from_second_basis(P, s0, s,
                                                   _interpolant_basis_callback, &dt);
        }
        flint_free(s0);
    }

    if (rdeg != NULL)
        for (slong i = 0; i < m; i++)
            rdeg[i] = s[i];

    if (sorted)
        nmod_mat_window_clear(Vs);
    else
        nmod_mat_clear(Vs);
    flint_free(spts);
    flint_free(s);
}
