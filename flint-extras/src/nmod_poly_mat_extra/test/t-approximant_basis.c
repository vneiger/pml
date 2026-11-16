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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_extra/impl.h"

/* p F[:,b] = 0 mod x^orders[b] for every row p of P */
/* TODO make general available interface for testing approximants with different orders? */
static int
_ab_is_approximant(const nmod_poly_mat_t P, const nmod_poly_mat_t F, const slong * orders)
{
    nmod_poly_mat_t PF;
    int ok = 1;

    nmod_poly_mat_init(PF, P->r, F->c, F->modulus);
    nmod_poly_mat_mul(PF, P, F);
    for (slong i = 0; i < PF->r && ok; i++)
        for (slong b = 0; b < PF->c && ok; b++)
        {
            const nmod_poly_struct * e = nmod_poly_mat_entry(PF, i, b);
            for (slong k = 0; k < FLINT_MIN(orders[b], e->length); k++)
                if (e->coeffs[k] != 0)
                    ok = 0;
        }
    nmod_poly_mat_clear(PF);
    return ok;
}

/* P is an approximant basis, shift-ordered weak Popov, with shifted row
 * degree rdeg equal to rdeg_ref: since P is shift-reduced, rows in the module
 * and the same shifted row degree as a basis imply that P is a basis */
static int
_ab_check(const nmod_poly_mat_t P, const slong * rdeg, const nmod_poly_mat_t F,
          const slong * orders, const slong * shift, const slong * rdeg_ref)
{
    const slong m = F->r;
    slong * rd = flint_malloc((m + 1) * sizeof(slong));
    int ok = 1;

    if (!nmod_poly_mat_is_ordered_weak_popov(P, shift, ROW_LOWER))
        ok = 0;
    if (ok && !_ab_is_approximant(P, F, orders))
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

/* orders: random, uniform, close (d or d - 1), or geometric (d, d/2, ...),
 * with some zeros */
static void
_ab_rand_orders(slong * orders, slong n, slong dmax, flint_rand_t state)
{
    const int kind = n_randint(state, 4);
    const slong d = 1 + n_randint(state, dmax);

    for (slong b = 0; b < n; b++)
    {
        if (kind == 0)
            orders[b] = n_randint(state, dmax + 1);
        else if (kind == 1)
            orders[b] = d;
        else if (kind == 2)
            orders[b] = FLINT_MAX(0, d - (slong) n_randint(state, 2));
        else
            orders[b] = d >> FLINT_MIN(b, FLINT_BITS - 2);
        if (n_randint(state, 8) == 0)
            orders[b] = 0;
    }
}

TEST_FUNCTION_START(nmod_poly_mat_approximant_basis, state)
{
    for (slong iter = 0; iter < 400 * flint_test_multiplier(); iter++)
    {
        slong m, n, dmax, D;
        ulong p;
        nmod_poly_mat_t F, P, Q, Pe, Fe;
        slong * orders, * shift, * rdeg, * rdeg_ref, * sh;
        int shift_kind, use_krylov_ref, ok = 1;
        _nmod_poly_mat_basis_dac_param_t prm;

        /* small m with large orders makes the divide and conquer split */
        if (n_randint(state, 2))
        {
            m = 1 + n_randint(state, 10);
            n = 1 + n_randint(state, 6);
            dmax = 1 + n_randint(state, 60);
        }
        else
        {
            m = 1 + n_randint(state, 4);
            n = 1 + n_randint(state, 4);
            dmax = 1 + n_randint(state, 200);
        }

        if (n_randint(state, 4) == 0)
            p = n_nth_prime(1 + n_randint(state, 5));
        else
            p = n_randtest_prime(state, 0);

        orders = flint_malloc(n * sizeof(slong));
        shift = flint_malloc(m * sizeof(slong));
        rdeg = flint_malloc(m * sizeof(slong));
        rdeg_ref = flint_malloc(m * sizeof(slong));
        sh = flint_malloc(m * sizeof(slong));

        _ab_rand_orders(orders, n, dmax, state);
        D = 0;
        dmax = 0;
        for (slong b = 0; b < n; b++)
        {
            D += orders[b];
            dmax = FLINT_MAX(dmax, orders[b]);
        }

        shift_kind = n_randint(state, 4);
        for (slong i = 0; i < m; i++)
        {
            if (shift_kind == 0)
                shift[i] = 0;
            else if (shift_kind == 1)
                shift[i] = (slong) n_randint(state, 41) - 20;
            else if (shift_kind == 2)
                shift[i] = i * (D + 1);           /* Hermite */
            else
                shift[i] = (m - 1 - i) * (D + 1); /* reverse Hermite */
        }

        nmod_poly_mat_init(F, m, n, p);
        nmod_poly_mat_init(Fe, m, n, p);
        nmod_poly_mat_init(P, m, m, p);
        nmod_poly_mat_init(Q, m, m, p);
        nmod_poly_mat_init(Pe, m, m, p);
        nmod_poly_mat_randtest(F, state, dmax + 1 + n_randint(state, 5));
        if (m > 1 && n_randint(state, 4) == 0)
            for (slong b = 0; b < n; b++)  /* dependent rows */
                nmod_poly_set(nmod_poly_mat_entry(F, m - 1, b), nmod_poly_mat_entry(F, 0, b));

        /* reference shifted row degree: pmbasis on the input raised to the
         * largest order */
        for (slong b = 0; b < n; b++)
            for (slong i = 0; i < m; i++)
            {
                nmod_poly_struct * g = nmod_poly_mat_entry(Fe, i, b);
                nmod_poly_set_trunc(g, nmod_poly_mat_entry(F, i, b), orders[b]);
                nmod_poly_shift_left(g, g, dmax - orders[b]);
            }
        for (slong i = 0; i < m; i++)
            rdeg_ref[i] = shift[i];
        if (dmax > 0)
            nmod_poly_mat_pmbasis(Pe, rdeg_ref, Fe, dmax);

        /* the s-Popov basis, by one Krylov call, when it is not too large */
        use_krylov_ref = (D <= 300);
        if (use_krylov_ref)
        {
            nmod_poly_mat_approximant_basis_krylov(Q, F, orders, shift);
            for (slong i = 0; i < m; i++)
                if (rdeg_ref[i] != shift[i] + nmod_poly_degree(nmod_poly_mat_entry(Q, i, i)))
                    ok = 0;
            if (!ok)
                TEST_FUNCTION_FAIL("reference: pmbasis and Krylov disagree, m = %wd, n = %wd, D = %wd\n", m, n, D);
        }

        /* general interface: ordered weak Popov, then Popov; shift NULL or
         * not, rdeg NULL, separate, or aliasing shift */
        for (int form = 0; form < 2 && ok; form++)
        {
            const int mode = n_randint(state, 3);
            const int zero = (shift_kind == 0) && n_randint(state, 2);

            for (slong i = 0; i < m; i++)
                sh[i] = shift[i];

            if (mode == 0)
                nmod_poly_mat_approximant_basis(P, NULL, F, orders, zero ? NULL : sh,
                                                form ? POPOV : ORD_WEAK_POPOV);
            else if (mode == 1)
                nmod_poly_mat_approximant_basis(P, rdeg, F, orders, zero ? NULL : sh,
                                                form ? POPOV : ORD_WEAK_POPOV);
            else
                nmod_poly_mat_approximant_basis(P, sh, F, orders, sh,
                                                form ? POPOV : ORD_WEAK_POPOV);

            for (slong i = 0; i < m && mode != 1; i++)
                rdeg[i] = rdeg_ref[i];  /* not returned (mode 0) or in sh */
            if (mode == 2)
                for (slong i = 0; i < m; i++)
                    rdeg[i] = sh[i];
            else if (mode == 0)
                for (slong i = 0; i < m; i++)
                    sh[i] = shift[i];

            if (!_ab_check(P, rdeg, F, orders, shift, rdeg_ref))
                ok = 0;
            if (ok && form && !nmod_poly_mat_is_popov(P, shift, ROW_LOWER))
                ok = 0;
            if (ok && form && use_krylov_ref && !nmod_poly_mat_equal(P, Q))
                ok = 0;
            if (!ok)
                TEST_FUNCTION_FAIL("general interface, form %d, mode %d: m = %wd, n = %wd, D = %wd, p = %wu\n",
                                   form, mode, m, n, D, p);
        }

        /* the divide and conquer with other parameters, so that all its
         * branches are used on small inputs: no Krylov call (cutoff 0) or
         * many, both splits, close orders or not */
        for (int variant = 0; variant < 4 && ok; variant++)
        {
            _nmod_poly_mat_basis_dac_param_init(prm);
            prm->leaf_cutoff = n_randint(state, 2) ? 0.0 : 8.0;
            prm->direct_cutoff = (variant & 1) ? 0.0 : 0.5 + n_randint(state, 64);
            prm->split = (variant >> 1) & 1;
            prm->raise_slack = n_randint(state, 2) ? -1.0 : 0.25 * n_randint(state, 8);

            for (slong i = 0; i < m; i++)
                sh[i] = shift[i];
            _nmod_poly_mat_approximant_basis(P, sh, F, orders, prm);

            if (!_ab_check(P, sh, F, orders, shift, rdeg_ref))
                TEST_FUNCTION_FAIL("divide and conquer, variant %d: m = %wd, n = %wd, D = %wd, p = %wu\n",
                                   variant, m, n, D, p);
        }

        nmod_poly_mat_clear(F);
        nmod_poly_mat_clear(Fe);
        nmod_poly_mat_clear(P);
        nmod_poly_mat_clear(Q);
        nmod_poly_mat_clear(Pe);
        flint_free(orders);
        flint_free(shift);
        flint_free(rdeg);
        flint_free(rdeg_ref);
        flint_free(sh);
    }

    TEST_FUNCTION_END(state);
}
