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

#include "nmod_poly_mat_inverse.h"

/* checks that X has length at most order and that X * A and A * X are the
 * identity modulo x^order (zero if order == 0) */
static int _invtest_check(const nmod_poly_mat_t X, const nmod_poly_mat_t A, slong order)
{
    const slong m = A->r;
    if (m == 0)
        return 1;
    if (nmod_poly_mat_max_length(X) > FLINT_MAX(order, 0))
        return 0;
    nmod_poly_mat_t P;
    nmod_poly_mat_init(P, m, m, A->modulus);
    int ok = 1;
    for (int side = 0; ok && side < 2; side++)
    {
        if (side == 0)
            nmod_poly_mat_mul(P, X, A);
        else
            nmod_poly_mat_mul(P, A, X);
        nmod_poly_mat_truncate(P, FLINT_MAX(order, 0));
        ok = (order <= 0) ? nmod_poly_mat_is_zero(P) : nmod_poly_mat_is_one(P);
    }
    nmod_poly_mat_clear(P);
    return ok;
}

/* random m x m matrix with entries of length < len; with some probability,
 * the constant coefficient is replaced by a random matrix of rank `rk0` */
static void _invtest_rand(nmod_poly_mat_t A, slong len, slong rk0, flint_rand_t state)
{
    const slong m = A->r;
    nmod_poly_mat_randtest(A, state, len);
    if (m > 0 && rk0 >= 0)
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

TEST_FUNCTION_START(nmod_poly_mat_inv_series, state)
{
    int i, result;

    for (i = 0; i < 400 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 62), 1);
        const slong m = n_randint(state, 7);
        /* input length: small, or large compared to the order */
        const slong len = 1 + (slong) n_randint(state, (n_randint(state, 3) == 0) ? 4 : 100);
        /* order: sometimes large, so that Newton iteration is used by the
         * general interface */
        const slong order = (n_randint(state, 4) == 0) ? (slong) (64 + n_randint(state, 100)) : (slong) n_randint(state, 40) - 2;
        /* constant coefficient: random, full rank, or rank deficient */
        const slong rk0 = (n_randint(state, 3) == 0) ? -1 : (n_randint(state, 4) == 0 ? (slong) n_randint(state, m + 1) : m);
        const slong cutoff = (slong) n_randint(state, 12) - 1;

        nmod_poly_mat_t A, X, Y, Z;
        nmod_poly_mat_init(A, m, m, prime);
        nmod_poly_mat_init(X, m, m, prime);
        nmod_poly_mat_init(Y, m, m, prime);
        nmod_poly_mat_init(Z, m, m, prime);
        _invtest_rand(A, len, rk0, state);

        nmod_mat_t inv0;
        nmod_mat_init(inv0, m, m, prime);
        nmod_poly_mat_get_coeff_mat(inv0, A, 0);
        const int invertible = nmod_mat_inv(inv0, inv0);
        const int expected = (invertible || order <= 0 || m == 0);

        /* general interface, basecase, Newton: return value, correctness,
         * same results */
        nmod_poly_mat_randtest(X, state, 3);
        nmod_poly_mat_set(Z, X);
        result = (nmod_poly_mat_inv_series(X, A, order) == expected);
        if (result && expected)
            result = _invtest_check(X, A, order);
        if (result && !expected)  /* imat not modified */
            result = nmod_poly_mat_equal(X, Z);
        if (result)
        {
            result = (nmod_poly_mat_inv_series_basecase(Y, A, order) == expected)
                     && (!expected || nmod_poly_mat_equal(X, Y));
        }
        if (result)
        {
            result = (nmod_poly_mat_inv_series_newton(Y, A, order, cutoff) == expected)
                     && (!expected || nmod_poly_mat_equal(X, Y));
        }

        /* aliasing */
        if (result)
        {
            nmod_poly_mat_set(Y, A);
            result = (nmod_poly_mat_inv_series(Y, Y, order) == expected)
                     && nmod_poly_mat_equal(Y, expected ? X : A);
        }
        if (result)
        {
            nmod_poly_mat_set(Y, A);
            result = (nmod_poly_mat_inv_series_newton(Y, Y, order, cutoff) == expected)
                     && nmod_poly_mat_equal(Y, expected ? X : A);
        }

        /* variants with the inverse of the constant coefficient */
        if (result && invertible)
        {
            _nmod_poly_mat_inv_series(Y, A, inv0, order);
            result = nmod_poly_mat_equal(X, Y);
            _nmod_poly_mat_inv_series_basecase(Y, A, inv0, order);
            result = result && nmod_poly_mat_equal(X, Y);
            _nmod_poly_mat_inv_series_newton(Y, A, inv0, order, cutoff);
            result = result && nmod_poly_mat_equal(X, Y);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, len = %wd, order = %wd, rk0 = %wd, cutoff = %wd\n",
                               prime, m, len, order, rk0, cutoff);

        nmod_mat_clear(inv0);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(X);
        nmod_poly_mat_clear(Y);
        nmod_poly_mat_clear(Z);
    }

    TEST_FUNCTION_END(state);
}
