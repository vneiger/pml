/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_interpolant.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_extra/impl.h"

/* P(pts[b]) V[:,b] = 0 for every column b */
static int
_ib_is_interpolant(const nmod_poly_mat_t P, const nmod_mat_t V, nn_srcptr pts)
{
    const slong m = V->r;
    nmod_mat_t E, v, w;
    int ok = 1;

    nmod_mat_init(E, m, m, V->mod.n);
    nmod_mat_init(w, m, 1, V->mod.n);
    for (slong b = 0; b < V->c && ok; b++)
    {
        nmod_poly_mat_evaluate_nmod(E, P, pts[b]);
        nmod_mat_window_init(v, V, 0, b, m, b + 1);
        nmod_mat_mul(w, E, v);
        nmod_mat_window_clear(v);
        if (!nmod_mat_is_zero(w))
            ok = 0;
    }
    nmod_mat_clear(w);
    nmod_mat_clear(E);
    return ok;
}

/* P is an interpolant basis, shift-ordered weak Popov, with shifted row
 * degree rdeg equal to rdeg_ref (see t-approximant_basis.c) */
static int
_ib_check(const nmod_poly_mat_t P, const slong * rdeg, const nmod_mat_t V,
          nn_srcptr pts, const slong * shift, const slong * rdeg_ref)
{
    const slong m = V->r;
    slong * rd = flint_malloc((m + 1) * sizeof(slong));
    int ok = 1;

    if (!nmod_poly_mat_is_ordered_weak_popov(P, shift, ROW_LOWER))
        ok = 0;
    if (ok && !_ib_is_interpolant(P, V, pts))
        ok = 0;
    if (ok)
    {
        nmod_poly_mat_row_degree(rd, P, shift);
        for (slong i = 0; i < m; i++)
            if (rd[i] != rdeg[i] || rdeg[i] != rdeg_ref[i])
                ok = 0;
    }
    flint_free(rd);
    return ok;
}

TEST_FUNCTION_START(nmod_poly_mat_interpolant_basis, state)
{
    for (slong iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong m, D, pool;
        ulong p;
        nmod_mat_t V, Vs;
        nmod_poly_mat_t P, Q;
        nn_ptr pts, spts;
        slong * shift, * rdeg, * rdeg_ref, * sh, * perm;
        int shift_kind, ok = 1;
        _nmod_poly_mat_basis_dac_param_t prm;

        if (n_randint(state, 2))
        {
            m = 1 + n_randint(state, 12);
            D = n_randint(state, 120);
        }
        else
        {
            m = 1 + n_randint(state, 4);
            D = n_randint(state, 250);
        }

        if (n_randint(state, 4) == 0)
            p = n_nth_prime(1 + n_randint(state, 5));
        else
            p = n_randtest_prime(state, 0);

        /* points from a pool: all distinct, few distinct, or one */
        if (n_randint(state, 3) == 0)
            pool = 1 + n_randint(state, 4);
        else
            pool = 1 + n_randint(state, D + 1);
        pool = FLINT_MIN(pool, (slong) FLINT_MIN(p, (ulong) WORD_MAX));

        pts = _nmod_vec_init(D + 1);
        spts = _nmod_vec_init(D + 1);
        perm = flint_malloc((D + 1) * sizeof(slong));
        {
            nn_ptr pl = _nmod_vec_init(pool);
            for (slong k = 0; k < pool; k++)
                pl[k] = n_randint(state, p);
            for (slong b = 0; b < D; b++)
                pts[b] = pl[n_randint(state, pool)];
            /* some inputs with as many columns at each point (pmintbasis) */
            if (n_randint(state, 4) == 0 && pool > 0)
                for (slong b = 0; b < D; b++)
                    pts[b] = pl[b % pool];
            _nmod_vec_clear(pl);
        }

        shift = flint_malloc(m * sizeof(slong));
        rdeg = flint_malloc(m * sizeof(slong));
        rdeg_ref = flint_malloc(m * sizeof(slong));
        sh = flint_malloc(m * sizeof(slong));
        shift_kind = n_randint(state, 3);
        for (slong i = 0; i < m; i++)
        {
            if (shift_kind == 0)
                shift[i] = 0;
            else if (shift_kind == 1)
                shift[i] = (slong) n_randint(state, 41) - 20;
            else
                shift[i] = i * (D + 1);
        }

        nmod_mat_init(V, m, D, p);
        nmod_mat_randtest(V, state);
        if (m > 1 && n_randint(state, 4) == 0)
            for (slong b = 0; b < D; b++)  /* dependent rows */
                nmod_mat_entry(V, m - 1, b) = nmod_mat_entry(V, 0, b);

        nmod_poly_mat_init(P, m, m, p);
        nmod_poly_mat_init(Q, m, m, p);

        /* reference: the s-Popov basis by one Krylov call */
        nmod_poly_mat_interpolant_basis_krylov(Q, V, pts, shift);
        for (slong i = 0; i < m; i++)
            rdeg_ref[i] = shift[i] + nmod_poly_degree(nmod_poly_mat_entry(Q, i, i));

        /* general interface */
        for (int form = 0; form < 2 && ok; form++)
        {
            const int mode = n_randint(state, 3);
            const int zero = (shift_kind == 0) && n_randint(state, 2);

            for (slong i = 0; i < m; i++)
                sh[i] = shift[i];

            if (mode == 0)
            {
                nmod_poly_mat_interpolant_basis(P, NULL, V, pts, zero ? NULL : sh,
                                                form ? POPOV : ORD_WEAK_POPOV);
                for (slong i = 0; i < m; i++)
                    rdeg[i] = rdeg_ref[i];
            }
            else if (mode == 1)
                nmod_poly_mat_interpolant_basis(P, rdeg, V, pts, zero ? NULL : sh,
                                                form ? POPOV : ORD_WEAK_POPOV);
            else
            {
                nmod_poly_mat_interpolant_basis(P, sh, V, pts, sh,
                                                form ? POPOV : ORD_WEAK_POPOV);
                for (slong i = 0; i < m; i++)
                    rdeg[i] = sh[i];
            }

            if (!_ib_check(P, rdeg, V, pts, shift, rdeg_ref))
                ok = 0;
            if (ok && form && !nmod_poly_mat_equal(P, Q))
                ok = 0;
            if (!ok)
                TEST_FUNCTION_FAIL("general interface, form %d, mode %d: m = %wd, D = %wd, pool = %wd, p = %wu\n",
                                   form, mode, m, D, pool, p);
        }

        /* the divide and conquer with other parameters, on the columns
         * sorted by point */
        for (slong b = 0; b < D; b++)
            perm[b] = b;
        for (slong b = 1; b < D; b++)  /* insertion sort by point */
        {
            slong c = perm[b], l = b;
            while (l > 0 && pts[perm[l - 1]] > pts[c])
            {
                perm[l] = perm[l - 1];
                l--;
            }
            perm[l] = c;
        }
        nmod_mat_init(Vs, m, D, p);
        for (slong b = 0; b < D; b++)
        {
            spts[b] = pts[perm[b]];
            for (slong i = 0; i < m; i++)
                nmod_mat_entry(Vs, i, b) = nmod_mat_entry(V, i, perm[b]);
        }

        for (int variant = 0; variant < 2 && ok; variant++)
        {
            _nmod_poly_mat_basis_dac_param_init(prm);
            prm->leaf_cutoff = n_randint(state, 2) ? 0.0 : 6.0;
            prm->direct_cutoff = variant ? 0.0 : 0.5 + n_randint(state, 64);

            for (slong i = 0; i < m; i++)
                sh[i] = shift[i];
            _nmod_poly_mat_interpolant_basis(P, sh, Vs, spts, prm);

            if (!_ib_check(P, sh, Vs, spts, shift, rdeg_ref))
                TEST_FUNCTION_FAIL("divide and conquer, variant %d: m = %wd, D = %wd, pool = %wd, p = %wu\n",
                                   variant, m, D, pool, p);
        }

        nmod_mat_clear(Vs);
        nmod_mat_clear(V);
        nmod_poly_mat_clear(P);
        nmod_poly_mat_clear(Q);
        _nmod_vec_clear(pts);
        _nmod_vec_clear(spts);
        flint_free(perm);
        flint_free(shift);
        flint_free(rdeg);
        flint_free(rdeg_ref);
        flint_free(sh);
    }

    TEST_FUNCTION_END(state);
}
