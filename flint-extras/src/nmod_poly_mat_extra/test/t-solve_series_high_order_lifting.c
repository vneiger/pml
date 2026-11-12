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
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_linsolve.h"

/* random m x m matrix with entries of length < len, whose constant
 * coefficient is replaced by a random matrix of rank `rk0` */
static void _holtest_rand(nmod_poly_mat_t A, slong len, slong rk0, flint_rand_t state)
{
    const slong m = A->r;
    nmod_poly_mat_randtest(A, state, len);
    if (m > 0)
    {
        nmod_mat_t C;
        nmod_mat_init(C, m, m, A->modulus);
        nmod_mat_randrank(C, state, rk0);
        if (n_randint(state, 2))
            nmod_mat_randops(C, state, m);
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, i, j), 0, nmod_mat_entry(C, i, j));
        nmod_mat_clear(C);
    }
}

TEST_FUNCTION_START(nmod_poly_mat_solve_series_high_order_lifting, state)
{
    int i, result;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        /* moduli from 2 bits (2 or 3) to 64 bits */
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        /* dimensions: mostly small, sometimes larger; n may be 0, or larger
         * than m */
        const int big = (n_randint(state, 10) == 0);
        const slong m = n_randint(state, big ? 16 : 7);
        const slong n = n_randint(state, big ? 16 : 5);
        const slong lenA = 1 + n_randint(state, big ? 12 : 10);
        /* precision: sometimes nonpositive, sometimes less than lenA */
        const slong prec = (n_randint(state, 8) == 0) ? (slong) n_randint(state, 3) - 2
                                                      : (slong) n_randint(state, big ? 200 : 100);
        /* length of B: shorter than A, up to twice the degree of A (one
         * chunk, block size length(B)), or longer (several chunks) */
        slong lenB;
        switch (n_randint(state, 3))
        {
            case 0: lenB = n_randint(state, lenA + 1); break;
            case 1: lenB = lenA + n_randint(state, lenA); break;
            default: lenB = 2 * lenA + n_randint(state, 10 * lenA);
        }
        /* rank of A(0): mostly full */
        const slong rk0 = (n_randint(state, 4) == 0) ? (slong) n_randint(state, m + 1) : m;

        nmod_poly_mat_t A, B, X, Y, Z;
        nmod_poly_mat_init(A, m, m, prime);
        nmod_poly_mat_init(B, n, m, prime);
        nmod_poly_mat_init(X, n, m, prime);
        nmod_poly_mat_init(Y, n, m, prime);
        nmod_poly_mat_init(Z, n, m, prime);
        _holtest_rand(A, lenA, rk0, state);
        nmod_poly_mat_randtest(B, state, lenB);

        const int expected = (rk0 == m || prec <= 0);

        /* return value, correctness (same as Dixon's lifting), X not
         * modified in case of failure */
        nmod_poly_mat_randtest(X, state, 3);
        nmod_poly_mat_set(Z, X);
        result = (nmod_poly_mat_solve_series_high_order_lifting(X, A, B, prec) == expected);
        if (result && expected)
            result = nmod_poly_mat_solve_series_dixon(Y, A, B, 0, prec)
                     && nmod_poly_mat_equal(X, Y);
        if (result && !expected)
            result = nmod_poly_mat_equal(X, Z);

        /* aliasing X and B */
        if (result)
        {
            nmod_poly_mat_set(Y, B);
            result = (nmod_poly_mat_solve_series_high_order_lifting(Y, A, Y, prec) == expected)
                     && nmod_poly_mat_equal(Y, expected ? X : B);
        }

        /* aliasing X and A */
        if (result && n == m)
        {
            nmod_poly_mat_set(Y, A);
            result = (nmod_poly_mat_solve_series_high_order_lifting(Y, Y, B, prec) == expected)
                     && nmod_poly_mat_equal(Y, expected ? X : A);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, n = %wd, lenA = %wd, lenB = %wd, prec = %wd, rk0 = %wd\n",
                               prime, m, n, lenA, lenB, prec, rk0);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(X);
        nmod_poly_mat_clear(Y);
        nmod_poly_mat_clear(Z);
    }

    TEST_FUNCTION_END(state);
}
