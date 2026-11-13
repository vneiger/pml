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
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_mat_extra/krylov_impl.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_interpolant.h"

static slong
_apbk_sum_diag_degrees(const nmod_poly_mat_t P)
{
    slong i, s = 0;
    for (i = 0; i < P->r; i++)
        s += nmod_poly_degree(nmod_poly_mat_entry(P, i, i));
    return s;
}

/* rank of the Krylov matrix of E with respect to the Jordan matrix */
static slong
_apbk_krylov_rank(const nmod_mat_t E, nn_srcptr eig, const slong * sizes, slong nb)
{
    _nmod_mat_krylov_op_t op;
    nmod_mat_t B;
    slong * rows, * exps, r;

    nmod_mat_init(B, E->c, E->c, E->mod.n);
    rows = flint_malloc((E->c + 1) * sizeof(slong));
    exps = flint_malloc((E->c + 1) * sizeof(slong));
    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nb, E->mod);
    r = _nmod_mat_krylov_basis_naive(B, rows, exps, E, op, NULL, NULL);
    _nmod_mat_krylov_op_clear(op);
    nmod_mat_clear(B);
    flint_free(rows);
    flint_free(exps);

    return r;
}

TEST_FUNCTION_START(nmod_poly_mat_approximant_basis_krylov, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong m, n, i, j, k, sigma, nb, off;
        ulong p;
        nmod_poly_mat_t F, P, PF;
        nmod_mat_t E;
        slong * orders, * shifts, * sizes;
        nn_ptr pts, eig;
        int use_shifts;

        m = n_randint(state, 6);
        n = n_randint(state, 4);

        if (n_randint(state, 3) == 0)
            p = n_nth_prime(1 + n_randint(state, 4));
        else
            p = n_randtest_prime(state, 0);

        nmod_poly_mat_init(F, m, n, p);
        nmod_poly_mat_init(P, m, m, p);
        nmod_poly_mat_init(PF, m, n, p);
        orders = flint_malloc((n + 1) * sizeof(slong));
        sizes = flint_malloc((n + 1) * sizeof(slong));
        shifts = flint_malloc((m + 1) * sizeof(slong));
        pts = _nmod_vec_init(n + 1);
        eig = _nmod_vec_init(n + 1);

        sigma = 0;
        for (j = 0; j < n; j++)
        {
            orders[j] = n_randint(state, 7);
            sigma += orders[j];
            pts[j] = n_randint(state, p);
        }
        if (n > 1 && n_randint(state, 2))
            pts[n - 1] = pts[0];   /* repeated point */

        for (i = 0; i < m; i++)
            shifts[i] = (slong) n_randint(state, 11) - 5;
        use_shifts = n_randint(state, 3) != 0;

        nmod_poly_mat_randtest(F, state, 1 + n_randint(state, 8));
        if (m > 1 && n_randint(state, 3) == 0)
            for (j = 0; j < n; j++)  /* dependent rows */
                nmod_poly_set(nmod_poly_mat_entry(F, m - 1, j), nmod_poly_mat_entry(F, 0, j));

        /* ------------------------------------------------------------ */
        /* Hermite-Pade approximation                                   */
        /* ------------------------------------------------------------ */
        nmod_poly_mat_approximant_basis_krylov(P, F, orders, use_shifts ? shifts : NULL);

        if (!nmod_poly_mat_is_popov(P, use_shifts ? shifts : NULL, ROW_LOWER))
            TEST_FUNCTION_FAIL("approximant: not Popov\nF = %{nmod_poly_mat}\n", F);

        nmod_poly_mat_mul(PF, P, F);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                for (k = 0; k < orders[j]; k++)
                    if (nmod_poly_get_coeff_ui(nmod_poly_mat_entry(PF, i, j), k) != 0)
                        TEST_FUNCTION_FAIL("approximant: P F != 0 mod x^order\n");

        nmod_mat_init(E, m, sigma, p);
        nb = 0;
        off = 0;
        for (j = 0; j < n; j++)
        {
            if (orders[j] == 0)
                continue;
            for (i = 0; i < m; i++)
                for (k = 0; k < orders[j]; k++)
                    nmod_mat_entry(E, i, off + k) =
                        nmod_poly_get_coeff_ui(nmod_poly_mat_entry(F, i, j), k);
            sizes[nb++] = orders[j];
            off += orders[j];
        }

        if (_apbk_sum_diag_degrees(P) != _apbk_krylov_rank(E, NULL, sizes, nb))
            TEST_FUNCTION_FAIL("approximant: not a basis (degree)\n");
        nmod_mat_clear(E);

        /* ------------------------------------------------------------ */
        /* M-Pade approximation                                         */
        /* ------------------------------------------------------------ */
        nmod_poly_mat_mpade_basis_krylov(P, F, pts, orders, use_shifts ? shifts : NULL);

        if (!nmod_poly_mat_is_popov(P, use_shifts ? shifts : NULL, ROW_LOWER))
            TEST_FUNCTION_FAIL("mpade: not Popov\n");

        nmod_poly_mat_mul(PF, P, F);
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                nmod_poly_t g;
                nmod_poly_init(g, p);
                nmod_poly_taylor_shift(g, nmod_poly_mat_entry(PF, i, j), pts[j]);
                for (k = 0; k < orders[j]; k++)
                    if (nmod_poly_get_coeff_ui(g, k) != 0)
                        TEST_FUNCTION_FAIL("mpade: P F != 0 mod (x - pt)^order\n");
                nmod_poly_clear(g);
            }
        }

        nmod_mat_init(E, m, sigma, p);
        nb = 0;
        off = 0;
        for (j = 0; j < n; j++)
        {
            if (orders[j] == 0)
                continue;
            for (i = 0; i < m; i++)
            {
                nmod_poly_t g;
                nmod_poly_init(g, p);
                nmod_poly_taylor_shift(g, nmod_poly_mat_entry(F, i, j), pts[j]);
                for (k = 0; k < orders[j]; k++)
                    nmod_mat_entry(E, i, off + k) = nmod_poly_get_coeff_ui(g, k);
                nmod_poly_clear(g);
            }
            eig[nb] = pts[j];
            sizes[nb++] = orders[j];
            off += orders[j];
        }

        if (_apbk_sum_diag_degrees(P) != _apbk_krylov_rank(E, eig, sizes, nb))
            TEST_FUNCTION_FAIL("mpade: not a basis (degree)\n");
        nmod_mat_clear(E);

        /* ------------------------------------------------------------ */
        /* interpolation at (possibly repeated) points                  */
        /* ------------------------------------------------------------ */
        {
            slong N = n_randint(state, 12);
            nmod_mat_t Vals;
            nn_ptr x = _nmod_vec_init(N + 1);

            nmod_mat_init(Vals, m, N, p);
            nmod_mat_randtest(Vals, state);
            for (j = 0; j < N; j++)
                x[j] = (j > 0 && n_randint(state, 3) == 0) ? x[j - 1] : n_randint(state, p);

            nmod_poly_mat_interpolant_basis_krylov(P, Vals, x, use_shifts ? shifts : NULL);

            if (!nmod_poly_mat_is_popov(P, use_shifts ? shifts : NULL, ROW_LOWER))
                TEST_FUNCTION_FAIL("interpolant: not Popov\n");

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < N; j++)
                {
                    ulong acc = 0;
                    for (k = 0; k < m; k++)
                        acc = nmod_add(acc, nmod_mul(nmod_poly_evaluate_nmod(
                            nmod_poly_mat_entry(P, i, k), x[j]),
                            nmod_mat_entry(Vals, k, j), Vals->mod), Vals->mod);
                    if (acc != 0)
                        TEST_FUNCTION_FAIL("interpolant: not an interpolant\n");
                }
            }

            if (_apbk_sum_diag_degrees(P) != _apbk_krylov_rank(Vals, x, NULL, N))
                TEST_FUNCTION_FAIL("interpolant: not a basis (degree)\n");

            nmod_mat_clear(Vals);
            _nmod_vec_clear(x);
        }

        nmod_poly_mat_clear(F);
        nmod_poly_mat_clear(P);
        nmod_poly_mat_clear(PF);
        flint_free(orders);
        flint_free(sizes);
        flint_free(shifts);
        _nmod_vec_clear(pts);
        _nmod_vec_clear(eig);
    }

    TEST_FUNCTION_END(state);
}
