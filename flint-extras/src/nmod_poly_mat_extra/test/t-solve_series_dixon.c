/*
    Copyright (C) 2025 Gilles Villard
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
#include "nmod_poly_mat_linsolve.h"
#include "nmod_poly_mat_multiply.h"

/* checks that X has length at most prec and that X A = B mod x^prec (X zero
 * if prec <= 0) */
static int _dixontest_check(const nmod_poly_mat_t X, const nmod_poly_mat_t A,
                            const nmod_poly_mat_t B, slong prec)
{
    if (prec <= 0)
        return nmod_poly_mat_is_zero(X);
    if (nmod_poly_mat_max_length(X) > prec)
        return 0;
    nmod_poly_mat_t P;
    nmod_poly_mat_init(P, B->r, B->c, B->modulus);
    nmod_poly_mat_mul(P, X, A);
    nmod_poly_mat_sub(P, P, B);
    nmod_poly_mat_truncate(P, prec);
    const int ok = nmod_poly_mat_is_zero(P);
    nmod_poly_mat_clear(P);
    return ok;
}

/* random m x m matrix with entries of length < len, whose constant
 * coefficient is replaced by a random matrix of rank `rk0` */
static void _dixontest_rand(nmod_poly_mat_t A, slong len, slong rk0, flint_rand_t state)
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

TEST_FUNCTION_START(nmod_poly_mat_solve_series_dixon, state)
{
    int i, result;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        /* moduli from 2 bits (2 or 3) to 64 bits */
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        /* dimensions: mostly small, sometimes larger; n may be 0, or larger
         * than m */
        const int big = (n_randint(state, 10) == 0);
        const slong m = n_randint(state, big ? 20 : 7);
        const slong n = n_randint(state, big ? 20 : 7);
        const slong lenA = 1 + n_randint(state, big ? 12 : 10);
        /* precision: sometimes nonpositive */
        const slong prec = (n_randint(state, 8) == 0) ? (slong) n_randint(state, 3) - 2
                                                      : (slong) n_randint(state, big ? 200 : 100);
        /* length of B: short, or at least prec */
        const slong lenB = (n_randint(state, 4) == 0) ? FLINT_MAX(prec, 0) + (slong) n_randint(state, 30)
                                                      : (slong) n_randint(state, 30);
        /* order: default (nonpositive), 1, random, at least prec */
        slong order;
        switch (n_randint(state, 4))
        {
            case 0: order = - (slong) n_randint(state, 2); break;
            case 1: order = 1; break;
            case 2: order = 1 + n_randint(state, 20); break;
            default: order = FLINT_MAX(prec, 0) + n_randint(state, 5);
        }
        /* rank of A(0): mostly full */
        const slong rk0 = (n_randint(state, 4) == 0) ? (slong) n_randint(state, m + 1) : m;

        nmod_poly_mat_t A, B, X, Y, Z;
        nmod_poly_mat_init(A, m, m, prime);
        nmod_poly_mat_init(B, n, m, prime);
        nmod_poly_mat_init(X, n, m, prime);
        nmod_poly_mat_init(Y, n, m, prime);
        nmod_poly_mat_init(Z, n, m, prime);
        _dixontest_rand(A, lenA, rk0, state);
        nmod_poly_mat_randtest(B, state, lenB);

        const int expected = (rk0 == m || prec <= 0);

        /* return value, correctness, X not modified in case of failure */
        nmod_poly_mat_randtest(X, state, 3);
        nmod_poly_mat_set(Z, X);
        result = (nmod_poly_mat_solve_series_dixon(X, A, B, order, prec) == expected);
        if (result && expected)
            result = _dixontest_check(X, A, B, prec);
        if (result && !expected)
            result = nmod_poly_mat_equal(X, Z);

        /* same as B times truncated inverse at order prec */
        if (result && expected && prec > 0 && m > 0)
        {
            nmod_poly_mat_t Ainv;
            nmod_poly_mat_init(Ainv, m, m, prime);
            nmod_poly_mat_inv_series(Ainv, A, prec);
            nmod_poly_mat_mulmid(Y, B, Ainv, 0, prec);
            result = nmod_poly_mat_equal(X, Y);
            nmod_poly_mat_clear(Ainv);
        }

        /* aliasing X and B */
        if (result)
        {
            nmod_poly_mat_set(Y, B);
            result = (nmod_poly_mat_solve_series_dixon(Y, A, Y, order, prec) == expected)
                     && nmod_poly_mat_equal(Y, expected ? X : B);
        }

        /* aliasing X and A */
        if (result && n == m)
        {
            nmod_poly_mat_set(Y, A);
            result = (nmod_poly_mat_solve_series_dixon(Y, Y, B, order, prec) == expected)
                     && nmod_poly_mat_equal(Y, expected ? X : A);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, n = %wd, lenA = %wd, lenB = %wd, prec = %wd, order = %wd, rk0 = %wd\n",
                               prime, m, n, lenA, lenB, prec, order, rk0);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(X);
        nmod_poly_mat_clear(Y);
        nmod_poly_mat_clear(Z);
    }

    TEST_FUNCTION_END(state);
}
