/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_approximant.h"
#include "impl.h"

/*
    Divide and conquer on D, the sum of the orders (see impl.h). At a node
    whose active columns (positive order) do not all have the same order and
    which is too large for a direct Krylov call, the first half is the
    problem at orders t_b <= d_b and the second half the residual problem at
    orders d_b - t_b, for the columns with d_b > t_b:

        G[:,b] = (x^(-t_b) P1 F[:,b]) mod x^(d_b - t_b),

    computed by middle products, one for each value of t_b and each group of
    windows d_b - t_b within a factor 2 (for the quasi-uniform cut
    t_b = min(d_b, tau), all t_b are equal to tau). The basis is P2 P1, as in
    pmbasis: it is shift-ordered weak Popov if P1 is, and P2 is for the
    shifted row degree of P1.
*/

/* tau = smallest value with sum_b min(d_b, tau) >= D/2 */
static slong
_quasi_uniform_cut(const slong * orders, slong n, slong D, slong dmax)
{
    slong lo = 1, hi = dmax;
    while (lo < hi)
    {
        slong mid = lo + (hi - lo) / 2, s = 0;
        for (slong b = 0; b < n; b++)
            s += FLINT_MIN(orders[b], mid);
        if (2 * s >= D)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

void
_nmod_poly_mat_approximant_basis(nmod_poly_mat_t P, slong * shift,
                                 const nmod_poly_mat_t F, const slong * orders,
                                 const _nmod_poly_mat_basis_dac_param_t prm)
{
    const slong m = F->r, n = F->c;
    slong nact = 0, D = 0, dmin = WORD_MAX, dmax = 0;
    slong * act;

    act = flint_malloc((n + 1) * sizeof(slong));
    for (slong b = 0; b < n; b++)
    {
        if (orders[b] > 0)
        {
            act[nact++] = b;
            D += orders[b];
            dmin = FLINT_MIN(dmin, orders[b]);
            dmax = FLINT_MAX(dmax, orders[b]);
        }
    }

    if (nact == 0)
    {
        nmod_poly_mat_one(P);
        flint_free(act);
        return;
    }

    /* one order: pmbasis, with its own Krylov leaf */
    if (dmin == dmax)
    {
        if (nact == n)
            _nmod_poly_mat_pmbasis(P, shift, F, dmax, prm->leaf_cutoff);
        else
        {
            nmod_poly_mat_t view;
            _nmod_poly_mat_view_init_submatrix(view, F, NULL, m, act, nact);
            _nmod_poly_mat_pmbasis(P, shift, view, dmax, prm->leaf_cutoff);
            _nmod_poly_mat_view_clear(view);
        }
        flint_free(act);
        return;
    }

    /* few conditions: direct Krylov call (columns of order 0 are skipped),
     * with ceil(D / dmax) conditions per level on average */
    if (_nmod_poly_mat_krylov_wanted(m, (D + dmax - 1) / dmax, D, prm->direct_cutoff))
    {
        _nmod_poly_mat_approximant_basis_krylov_leaf(P, shift, F, orders);
        flint_free(act);
        return;
    }

    /* close orders: pmbasis on x^(dmax - d_b) F[:,b], at order dmax */
    if (prm->raise_slack >= 0.0 && (double) nact * dmax <= (1.0 + prm->raise_slack) * D)
    {
        nmod_poly_mat_t G;
        nmod_poly_mat_init(G, m, nact, F->modulus);
        for (slong k = 0; k < nact; k++)
        {
            const slong b = act[k];
            for (slong i = 0; i < m; i++)
            {
                nmod_poly_struct * g = nmod_poly_mat_entry(G, i, k);
                nmod_poly_set_trunc(g, nmod_poly_mat_entry(F, i, b), orders[b]);
                nmod_poly_shift_left(g, g, dmax - orders[b]);
            }
        }
        _nmod_poly_mat_pmbasis(P, shift, G, dmax, prm->leaf_cutoff);
        nmod_poly_mat_clear(G);
        flint_free(act);
        return;
    }

    /* split: first half at orders t_b, second half on the residual */
    {
        slong * t = flint_malloc(n * sizeof(slong));
        slong * rem = flint_malloc(nact * sizeof(slong));
        slong * r = flint_malloc(nact * sizeof(slong));
        slong nrem = 0;
        nmod_poly_mat_t P1, P2, G;

        if (prm->split == 0)
        {
            const slong tau = _quasi_uniform_cut(orders, n, D, dmax);
            for (slong b = 0; b < n; b++)
                t[b] = FLINT_MIN(FLINT_MAX(orders[b], 0), tau);
        }
        else
        {
            for (slong b = 0; b < n; b++)
                t[b] = FLINT_MAX(orders[b], 0) / 2;
        }

        nmod_poly_mat_init(P1, m, m, F->modulus);
        _nmod_poly_mat_approximant_basis(P1, shift, F, t, prm);

        /* remaining columns, by increasing t_b and then decreasing d_b; one
         * middle product for each group of contiguous columns with the same
         * t_b whose windows d_b - t_b are at least half the first one (the
         * windows are computed at the length of the first one) */
        for (slong k = 0; k < nact; k++)
            if (orders[act[k]] > t[act[k]])
                rem[nrem++] = act[k];
        for (slong k = 1; k < nrem; k++)
        {
            const slong b = rem[k];
            slong l = k;
            while (l > 0 && (t[rem[l - 1]] > t[b]
                             || (t[rem[l - 1]] == t[b] && orders[rem[l - 1]] < orders[b])))
            {
                rem[l] = rem[l - 1];
                l--;
            }
            rem[l] = b;
        }

        nmod_poly_mat_init(G, m, nrem, F->modulus);
        for (slong k = 0; k < nrem; )
        {
            const slong tk = t[rem[k]], nhi = orders[rem[k]];
            slong l = k + 1;
            nmod_poly_mat_t view, R;

            while (l < nrem && t[rem[l]] == tk
                   && 2 * (orders[rem[l]] - tk) >= nhi - tk)
                l++;

            _nmod_poly_mat_view_init_submatrix(view, F, NULL, m, rem + k, l - k);
            nmod_poly_mat_init(R, m, l - k, F->modulus);
            nmod_poly_mat_mulmid(R, P1, view, tk, nhi);
            for (slong j = k; j < l; j++)
            {
                r[j] = orders[rem[j]] - tk;
                for (slong i = 0; i < m; i++)
                    nmod_poly_swap(nmod_poly_mat_entry(G, i, j),
                                   nmod_poly_mat_entry(R, i, j - k));
            }
            nmod_poly_mat_clear(R);
            _nmod_poly_mat_view_clear(view);
            k = l;
        }

        nmod_poly_mat_init(P2, m, m, F->modulus);
        _nmod_poly_mat_approximant_basis(P2, shift, G, r, prm);
        nmod_poly_mat_multiply(P, P2, P1);

        nmod_poly_mat_clear(G);
        nmod_poly_mat_clear(P2);
        nmod_poly_mat_clear(P1);
        flint_free(r);
        flint_free(rem);
        flint_free(t);
    }

    flint_free(act);
}

/* Popov form, from a second basis computed with the shift -delta */
void
_nmod_poly_mat_popov_from_second_basis(nmod_poly_mat_t P,
        const slong * shift, const slong * rdeg,
        void (* basis)(nmod_poly_mat_t, slong *, const void *), const void * data)
{
    const slong m = P->r;
    slong * mdelta = flint_malloc((m + 1) * sizeof(slong));
    slong * tmp = flint_malloc((m + 1) * sizeof(slong));
    nmod_poly_mat_t R;
    nmod_mat_t L, C;
    nmod_mat_poly_t Rp;

    /* P is shift-ordered weak Popov: the pivot of row i is P[i,i], of degree
     * delta_i = rdeg_i - shift_i */
    for (slong i = 0; i < m; i++)
    {
        mdelta[i] = (shift ? shift[i] : 0) - rdeg[i];
        tmp[i] = mdelta[i];
    }

    nmod_poly_mat_init(R, m, m, P->modulus);
    basis(R, tmp, data);

    nmod_mat_init(L, m, m, P->modulus);
    nmod_mat_init(C, m, m, P->modulus);
    nmod_poly_mat_leading_matrix(C, R, mdelta, ROW_LOWER);
    if (!nmod_mat_inv(L, C))
        flint_throw(FLINT_ERROR, "_nmod_poly_mat_popov_from_second_basis: singular leading matrix\n");

    /* P = L^(-1) R, coefficient by coefficient */
    nmod_mat_poly_init(Rp, m, m, P->modulus);
    nmod_mat_poly_set_from_poly_mat(Rp, R);
    for (slong k = 0; k < Rp->length; k++)
    {
        nmod_mat_t Rk;  /* view, writes go to Rp */
        nmod_mat_poly_coeff_attach(Rk, Rp, k);
        nmod_mat_mul(C, L, Rk);
        nmod_mat_set(Rk, C);
    }
    nmod_poly_mat_set_from_mat_poly(P, Rp);

    nmod_mat_clear(C);
    nmod_mat_poly_clear(Rp);
    nmod_mat_clear(L);
    nmod_poly_mat_clear(R);
    flint_free(tmp);
    flint_free(mdelta);
}

typedef struct
{
    const nmod_poly_mat_struct * F;
    const slong * orders;
    const _nmod_poly_mat_basis_dac_param_struct * prm;
}
_approximant_basis_data_struct;

static void
_approximant_basis_callback(nmod_poly_mat_t P, slong * shift, const void * data)
{
    const _approximant_basis_data_struct * dt = data;
    _nmod_poly_mat_approximant_basis(P, shift, dt->F, dt->orders, dt->prm);
}

void
nmod_poly_mat_approximant_basis(nmod_poly_mat_t P, slong * rdeg,
                                const nmod_poly_mat_t F, const slong * orders,
                                const slong * shift, poly_mat_form_t form)
{
    const slong m = F->r;
    slong * s;
    _nmod_poly_mat_basis_dac_param_t prm;

    if (P->r != m || P->c != m)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_approximant_basis: P must be F->r x F->r\n");
    if (form > POPOV)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_approximant_basis: form must be at most POPOV\n");
    for (slong b = 0; b < F->c; b++)
        if (orders[b] < 0)
            flint_throw(FLINT_ERROR, "nmod_poly_mat_approximant_basis: negative order\n");

    if (m == 0)
        return;

    /* s: shift on input, shifted row degree on output; shift and rdeg are
     * read and written only here, so they may alias */
    s = flint_malloc(m * sizeof(slong));
    for (slong i = 0; i < m; i++)
        s[i] = shift ? shift[i] : 0;

    _nmod_poly_mat_basis_dac_param_init(prm);
    _nmod_poly_mat_approximant_basis(P, s, F, orders, prm);

    if (form == POPOV)
    {
        slong * s0 = flint_malloc(m * sizeof(slong));
        for (slong i = 0; i < m; i++)
            s0[i] = shift ? shift[i] : 0;

        if (!nmod_poly_mat_is_popov(P, s0, ROW_LOWER))
        {
            _approximant_basis_data_struct dt = {F, orders, prm};
            _nmod_poly_mat_popov_from_second_basis(P, s0, s,
                                                   _approximant_basis_callback, &dt);
        }
        flint_free(s0);
    }

    if (rdeg != NULL)
        for (slong i = 0; i < m; i++)
            rdeg[i] = s[i];

    flint_free(s);
}
