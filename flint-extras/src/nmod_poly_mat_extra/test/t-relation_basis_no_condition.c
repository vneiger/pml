/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Approximant and interpolant bases without condition: order 0 (or all
 * orders 0, or no column) and no point (or matrices E_k without column).
 * Every entry point returns the identity and leaves the shift unchanged,
 * without reading pts and E when there is no point (they are NULL here); see
 * "Order zero" in nmod_poly_mat_approximant.h and "No point" in
 * nmod_poly_mat_interpolant.h. Also: for the general interface in Popov form,
 * the columns of order 0 can be removed without changing the output. */

#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_poly_mat_extra/impl.h"  /* for _nmod_poly_mat_pmbasis, _pmintbasis */

/* P is the identity and sh = sh0 */
static int
_nc_is_trivial(const nmod_poly_mat_t P, const slong * sh, const slong * sh0)
{
    if (!nmod_poly_mat_is_one(P))
        return 0;
    for (slong i = 0; i < P->r; i++)
        if (sh[i] != sh0[i])
            return 0;
    return 1;
}

#define _NC_CHECK(name, call)                                                 \
do {                                                                          \
    for (slong i = 0; i < m; i++)                                             \
        sh[i] = sh0[i];                                                       \
    nmod_poly_mat_zero(P);                                                    \
    call;                                                                     \
    if (!_nc_is_trivial(P, sh, sh0))                                          \
        TEST_FUNCTION_FAIL("%s: m = %wd, n = %wd, order = %wd, p = %wu\n",    \
                           name, m, n, order, p);                             \
} while (0)

TEST_FUNCTION_START(nmod_poly_mat_relation_basis_no_condition, state)
{
    for (slong iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        const slong m = 1 + n_randint(state, 8);
        /* odd iterations: order 0 (and no point); even ones: no column */
        const slong n = (iter % 2) ? (slong) n_randint(state, 6) : 0;
        const slong order = (iter % 2) ? 0 : 1 + (slong) n_randint(state, 300);
        const ulong p = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong thres = 1 + n_randint(state, 16);

        nmod_poly_mat_t F, P;
        nmod_poly_mat_init(F, m, n, p);
        nmod_poly_mat_randtest(F, state, 1 + n_randint(state, 20));
        nmod_poly_mat_init(P, m, m, p);

        slong * sh0 = flint_malloc(m * sizeof(slong));
        slong * sh = flint_malloc(m * sizeof(slong));
        for (slong i = 0; i < m; i++)
            sh0[i] = (slong) n_randint(state, 20) - 10;

        slong * orders = flint_malloc((n + 1) * sizeof(slong));
        ulong * pts = flint_malloc((n + 1) * sizeof(ulong));
        for (slong b = 0; b < n; b++)
        {
            orders[b] = order;
            pts[b] = n_randint(state, p);
        }

        /* approximant bases */
        _NC_CHECK("nmod_poly_mat_mbasis", nmod_poly_mat_mbasis(P, sh, F, order));
        _NC_CHECK("nmod_poly_mat_pmbasis", nmod_poly_mat_pmbasis(P, sh, F, order));
        _NC_CHECK("_nmod_poly_mat_pmbasis",
                  _nmod_poly_mat_pmbasis(P, sh, F, order, 6.0 * n_randint(state, 2), thres));
        {
            nmod_mat_poly_t A, Fm;
            nmod_mat_poly_init(A, m, m, p);
            nmod_mat_poly_init(Fm, m, n, p);
            nmod_mat_poly_set_from_poly_mat(Fm, F);
            _NC_CHECK("nmod_mat_poly_mbasis",
                      nmod_mat_poly_mbasis(A, sh, Fm, order);
                      nmod_poly_mat_set_from_mat_poly(P, A));
            nmod_mat_poly_clear(A);
            nmod_mat_poly_clear(Fm);
        }
        for (int form = NONE; form <= POPOV; form++)
        {
            /* rdeg aliasing shift; and the zero shift (rdeg must be zero,
               sh0 is added back for the check) */
            _NC_CHECK("nmod_poly_mat_approximant_basis",
                      nmod_poly_mat_approximant_basis(P, sh, F, orders, sh, form));
            _NC_CHECK("nmod_poly_mat_approximant_basis (zero shift)",
                      nmod_poly_mat_approximant_basis(P, sh, F, orders, NULL, form);
                      for (slong i = 0; i < m; i++) sh[i] += sh0[i]);
        }
        _NC_CHECK("nmod_poly_mat_approximant_basis_krylov",
                  nmod_poly_mat_approximant_basis_krylov(P, F, orders, sh));
        _NC_CHECK("nmod_poly_mat_mpade_basis_krylov",
                  nmod_poly_mat_mpade_basis_krylov(P, F, pts, orders, sh));

        nmod_poly_mat_one(P);
        if (!nmod_poly_mat_is_approximant_basis(P, F, order, sh0, ROW_LOWER))
            TEST_FUNCTION_FAIL("nmod_poly_mat_is_approximant_basis: m = %wd, n = %wd, order = %wd\n",
                               m, n, order);

        /* interpolant bases: no point (pts and E are NULL), or d points
           with matrices E_k of size m x 0 */
        if (order == 0)
        {
            nmod_mat_t V;
            nmod_mat_init(V, m, 0, p);

            _NC_CHECK("nmod_poly_mat_mintbasis",
                      nmod_poly_mat_mintbasis(P, sh, NULL, NULL, 0));
            _NC_CHECK("nmod_poly_mat_pmintbasis",
                      nmod_poly_mat_pmintbasis(P, sh, NULL, NULL, 0));
            _NC_CHECK("_nmod_poly_mat_pmintbasis",
                      _nmod_poly_mat_pmintbasis(P, sh, NULL, NULL, 0, 6.0, thres));
            /* r = 0 is not read either */
            _NC_CHECK("nmod_poly_mat_pmintbasis_geometric",
                      nmod_poly_mat_pmintbasis_geometric(P, sh, NULL, NULL, 0, 0));
            _NC_CHECK("nmod_poly_mat_pmintbasis_geometric_auto",
                      nmod_poly_mat_pmintbasis_geometric_auto(P, sh, NULL, NULL, 0));
            for (int form = NONE; form <= POPOV; form++)
                _NC_CHECK("nmod_poly_mat_interpolant_basis",
                          nmod_poly_mat_interpolant_basis(P, sh, V, NULL, sh, form));
            _NC_CHECK("nmod_poly_mat_interpolant_basis_krylov",
                      nmod_poly_mat_interpolant_basis_krylov(P, V, NULL, sh));

            nmod_poly_mat_one(P);
            if (!nmod_poly_mat_is_interpolant_basis(P, NULL, NULL, 0, sh0, ROW_LOWER))
                TEST_FUNCTION_FAIL("nmod_poly_mat_is_interpolant_basis: m = %wd\n", m);

            nmod_mat_clear(V);
        }
        else if (p > 2 * (ulong) order + 1)
        {
            const slong d = order;
            ulong * ipts = flint_malloc(d * sizeof(ulong));
            nmod_mat_struct * E = flint_malloc(d * sizeof(nmod_mat_struct));
            for (slong k = 0; k < d; k++)
            {
                ipts[k] = n_randint(state, p);
                nmod_mat_init(E + k, m, 0, p);
            }

            _NC_CHECK("nmod_poly_mat_mintbasis (no column)",
                      nmod_poly_mat_mintbasis(P, sh, ipts, E, d));
            _NC_CHECK("nmod_poly_mat_pmintbasis (no column)",
                      nmod_poly_mat_pmintbasis(P, sh, ipts, E, d));
            _NC_CHECK("_nmod_poly_mat_pmintbasis (no column)",
                      _nmod_poly_mat_pmintbasis(P, sh, ipts, E, d, 6.0, thres));
            _NC_CHECK("nmod_poly_mat_pmintbasis_geometric_auto (no column)",
                      nmod_poly_mat_pmintbasis_geometric_auto(P, sh, NULL, E, d));

            for (slong k = 0; k < d; k++)
                nmod_mat_clear(E + k);
            flint_free(E);
            flint_free(ipts);
        }

        flint_free(pts);
        flint_free(orders);
        flint_free(sh);
        flint_free(sh0);
        nmod_poly_mat_clear(P);
        nmod_poly_mat_clear(F);
    }

    /* general interface in Popov form: columns of order 0 put no condition,
       removing them does not change the (unique) basis */
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        const slong m = 1 + n_randint(state, 8);
        const slong n = 1 + n_randint(state, 6);
        const ulong p = n_randprime(state, 2 + n_randint(state, 63), 1);

        nmod_poly_mat_t F, F1, P, P1;
        nmod_poly_mat_init(F, m, n, p);
        nmod_poly_mat_randtest(F, state, 1 + n_randint(state, 20));
        nmod_poly_mat_init(P, m, m, p);
        nmod_poly_mat_init(P1, m, m, p);

        slong * orders = flint_malloc(n * sizeof(slong));
        slong * orders1 = flint_malloc(n * sizeof(slong));
        slong * cols = flint_malloc(n * sizeof(slong));
        slong * shift = flint_malloc(m * sizeof(slong));
        slong * rdeg = flint_malloc(m * sizeof(slong));
        slong * rdeg1 = flint_malloc(m * sizeof(slong));
        slong n1 = 0;
        for (slong b = 0; b < n; b++)
        {
            orders[b] = n_randint(state, 2) ? 0 : (slong) n_randint(state, 40);
            if (orders[b] > 0)
            {
                cols[n1] = b;
                orders1[n1++] = orders[b];
            }
        }
        for (slong i = 0; i < m; i++)
            shift[i] = (slong) n_randint(state, 20) - 10;

        nmod_poly_mat_init(F1, m, n1, p);
        for (slong i = 0; i < m; i++)
            for (slong k = 0; k < n1; k++)
                nmod_poly_set(nmod_poly_mat_entry(F1, i, k),
                              nmod_poly_mat_entry(F, i, cols[k]));

        nmod_poly_mat_approximant_basis(P, rdeg, F, orders, shift, POPOV);
        nmod_poly_mat_approximant_basis(P1, rdeg1, F1, orders1, shift, POPOV);

        int ok = nmod_poly_mat_equal(P, P1);
        for (slong i = 0; i < m; i++)
            if (rdeg[i] != rdeg1[i])
                ok = 0;
        if (!ok)
            TEST_FUNCTION_FAIL("columns of order 0: m = %wd, n = %wd (%wd of positive order), p = %wu\n",
                               m, n, n1, p);

        flint_free(rdeg1);
        flint_free(rdeg);
        flint_free(shift);
        flint_free(cols);
        flint_free(orders1);
        flint_free(orders);
        nmod_poly_mat_clear(P1);
        nmod_poly_mat_clear(P);
        nmod_poly_mat_clear(F1);
        nmod_poly_mat_clear(F);
    }

    TEST_FUNCTION_END(state);
}

#undef _NC_CHECK
