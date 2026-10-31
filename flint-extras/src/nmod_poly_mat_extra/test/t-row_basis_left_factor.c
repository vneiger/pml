/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_row_basis.h"

/* random m x n input: random, or product of random m x k and k x n
 * matrices (rank <= k), with possibly some zero rows */
static void _lftest_rand_input(nmod_poly_mat_t A, slong len, flint_rand_t state)
{
    const slong m = A->r;
    const slong n = A->c;
    if (n_randint(state, 2) == 0 || m == 0 || n == 0)
        nmod_poly_mat_randtest(A, state, len);
    else
    {
        const slong k = n_randint(state, FLINT_MIN(m, n) + 1);
        nmod_poly_mat_t W, B;
        nmod_poly_mat_init(W, m, k, A->modulus);
        nmod_poly_mat_init(B, k, n, A->modulus);
        nmod_poly_mat_randtest(W, state, 1 + n_randint(state, len));
        nmod_poly_mat_randtest(B, state, len);
        nmod_poly_mat_mul(A, W, B);
        nmod_poly_mat_clear(W);
        nmod_poly_mat_clear(B);
    }
    for (slong i = 0; i < m; i++)
        if (n_randint(state, 6) == 0)
            for (slong j = 0; j < n; j++)
                nmod_poly_zero(nmod_poly_mat_entry(A, i, j));
}

/* checks lfac (first r columns) times the first r rows of rbas is A, and
 * other columns of lfac are zero */
static int _lftest_verify(const nmod_poly_mat_t A, const nmod_poly_mat_t lfac,
                          const nmod_poly_mat_t rbas, slong r)
{
    for (slong i = 0; i < lfac->r; i++)
        for (slong j = r; j < lfac->c; j++)
            if (!nmod_poly_is_zero(nmod_poly_mat_entry(lfac, i, j)))
                return 0;
    if (r == 0)
        return nmod_poly_mat_is_zero(A);
    nmod_poly_mat_t C, R, prod;
    nmod_poly_mat_window_init(C, lfac, 0, 0, A->r, r);
    nmod_poly_mat_window_init(R, rbas, 0, 0, r, A->c);
    nmod_poly_mat_init(prod, A->r, A->c, A->modulus);
    nmod_poly_mat_mul(prod, C, R);
    const int res = nmod_poly_mat_equal(prod, A);
    nmod_poly_mat_clear(prod);
    nmod_poly_mat_window_clear(R);
    nmod_poly_mat_window_clear(C);
    return res;
}

/* row basis of A in shift-weak Popov form (shift may be NULL), in the first
 * rows of rbas (n x n, extra rows zero); returns the rank */
static slong _lftest_weak_popov_row_basis(nmod_poly_mat_t rbas, const nmod_poly_mat_t A,
                                          const slong * shift)
{
    nmod_poly_mat_t W;
    nmod_poly_mat_init_set(W, A);
    slong * pivind = flint_malloc(FLINT_MAX(1, A->r) * sizeof(slong));
    const slong r = (A->r == 0 || A->c == 0) ? 0 :
                    nmod_poly_mat_weak_popov_iter(W, shift, NULL, pivind, NULL, ROW_LOWER);
    nmod_poly_mat_zero(rbas);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < A->c; j++)
            nmod_poly_swap(nmod_poly_mat_entry(rbas, i, j), nmod_poly_mat_entry(W, i, j));
    flint_free(pivind);
    nmod_poly_mat_clear(W);
    return r;
}

TEST_FUNCTION_START(nmod_poly_mat_row_basis_left_factor, state)
{
    int i, result;

    for (i = 0; i < 400 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 61), 1);
        const slong m = n_randint(state, 9);
        const slong n = n_randint(state, 7);
        const slong len = 1 + n_randint(state, 8);

        nmod_poly_mat_t A, rbas, lfac, lfac2;
        nmod_poly_mat_init(A, m, n, prime);
        nmod_poly_mat_init(rbas, n, n, prime);
        nmod_poly_mat_init(lfac, m, n, prime);
        nmod_poly_mat_init(lfac2, m, n, prime);
        _lftest_rand_input(A, len, state);

        slong * shift = flint_malloc(FLINT_MAX(1, n) * sizeof(slong));
        for (slong k = 0; k < n; k++)
            shift[k] = (slong) n_randint(state, 3 * len) - len;
        const slong * sh = (n_randint(state, 2) == 0) ? NULL : shift;

        /* row basis in sh-weak Popov form, then left factor with and
         * without verification */
        const slong r = _lftest_weak_popov_row_basis(rbas, A, sh);
        result = nmod_poly_mat_row_basis_left_factor(lfac, A, rbas, r, sh, 1)
                 && _lftest_verify(A, lfac, rbas, r);
        if (result)
            result = nmod_poly_mat_row_basis_left_factor(lfac2, A, rbas, r, sh, 0)
                     && nmod_poly_mat_equal(lfac, lfac2);

        /* the basis is reduced for the shift used, but generally not for
         * another one: then either 0 is returned, or the result is correct */
        if (result && sh != NULL)
        {
            int res2 = nmod_poly_mat_row_basis_left_factor(lfac2, A, rbas, r, NULL, 1);
            if (res2 && !_lftest_verify(A, lfac2, rbas, r))
                result = 0;
        }

        /* rows outside the row span are detected when verifying: append to
         * A a row which is not in the row span of x * R */
        if (result && r > 0)
        {
            nmod_poly_mat_t xR, B;
            nmod_poly_mat_init(xR, n, n, prime);
            nmod_poly_mat_init(B, m + 1, n, prime);
            nmod_poly_mat_shift_left(xR, rbas, 1);
            for (slong ii = 0; ii < m; ii++)
                for (slong j = 0; j < n; j++)
                    nmod_poly_shift_left(nmod_poly_mat_entry(B, ii, j),
                                         nmod_poly_mat_entry(A, ii, j), 1);
            /* last row: first row of R (not divisible by x in the row span
             * of x * R, since R has no zero row) */
            for (slong j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(B, m, j), nmod_poly_mat_entry(rbas, 0, j));
            nmod_poly_mat_t lf3;
            nmod_poly_mat_init(lf3, m + 1, n, prime);
            if (nmod_poly_mat_row_basis_left_factor(lf3, B, xR, r, sh, 1))
            {
                flint_printf("FAIL: row outside the row span not detected\n");
                result = 0;
            }
            /* without the last row, it is in the row span of x * R */
            nmod_poly_mat_t Bw;
            nmod_poly_mat_window_init(Bw, B, 0, 0, m, n);
            nmod_poly_mat_t lf4;
            nmod_poly_mat_init(lf4, m, n, prime);
            if (!nmod_poly_mat_row_basis_left_factor(lf4, Bw, xR, r, sh, 1)
                || !nmod_poly_mat_equal(lf4, lfac))
            {
                flint_printf("FAIL: left factor for x * A and x * R\n");
                result = 0;
            }
            nmod_poly_mat_clear(lf4);
            nmod_poly_mat_window_clear(Bw);
            nmod_poly_mat_clear(lf3);
            nmod_poly_mat_clear(B);
            nmod_poly_mat_clear(xR);
        }

        if (!result)
        {
            flint_printf("input (%wd x %wd, rank %wd):\n", m, n, r);
            nmod_poly_mat_print_pretty(A, "x");
            flint_printf("\n");
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, n = %wd, len = %wd, shifted = %d\n",
                               prime, m, n, len, sh != NULL);
        }

        flint_free(shift);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(rbas);
        nmod_poly_mat_clear(lfac);
        nmod_poly_mat_clear(lfac2);
    }

    /* non-reduced basis: R = [[1, x], [1, x + 1]] has singular leading
     * matrix, 0 must be returned */
    {
        nmod_poly_mat_t R, A, C;
        nmod_poly_mat_init(R, 2, 2, 7);
        nmod_poly_mat_init(A, 1, 2, 7);
        nmod_poly_mat_init(C, 1, 2, 7);
        nmod_poly_one(nmod_poly_mat_entry(R, 0, 0));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(R, 0, 1), 1, 1);
        nmod_poly_one(nmod_poly_mat_entry(R, 1, 0));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(R, 1, 1), 1, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(R, 1, 1), 0, 1);
        nmod_poly_one(nmod_poly_mat_entry(A, 0, 0));
        if (nmod_poly_mat_row_basis_left_factor(C, A, R, 2, NULL, 0))
            TEST_FUNCTION_FAIL("non-reduced basis not detected\n");
        nmod_poly_mat_clear(R);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
