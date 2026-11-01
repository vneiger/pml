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

#include "nmod_poly_mat_utils.h"

/* rank of pmat evaluated at pt, computed naively */
static slong _rktest_rank_at(const nmod_poly_mat_t pmat, ulong pt)
{
    nmod_mat_t M;
    nmod_mat_init(M, pmat->r, pmat->c, pmat->modulus);
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_mat_entry(M, i, j) = nmod_poly_evaluate_nmod(nmod_poly_mat_entry(pmat, i, j), pt);
    const slong rk = nmod_mat_rank(M);
    nmod_mat_clear(M);
    return rk;
}

/* random input: random, product of random m x k and k x n matrices (rank
 * <= k), with possibly zero rows, or (when square) x*I - M with M singular */
static void _rktest_rand_input(nmod_poly_mat_t A, slong len, flint_rand_t state)
{
    const slong m = A->r;
    const slong n = A->c;
    const ulong type = n_randint(state, 3);

    if (type == 0 || m == 0 || n == 0)
        nmod_poly_mat_randtest(A, state, len);
    else if (type == 1 || m != n)
    {
        const slong k = n_randint(state, FLINT_MIN(m, n) + 1);
        nmod_poly_mat_t W, B;
        nmod_poly_mat_init(W, m, k, A->modulus);
        nmod_poly_mat_init(B, k, n, A->modulus);
        nmod_poly_mat_randtest(W, state, len);
        nmod_poly_mat_randtest(B, state, len);
        nmod_poly_mat_mul(A, W, B);
        nmod_poly_mat_clear(W);
        nmod_poly_mat_clear(B);
    }
    else
    {
        /* x*I - M, M random with a zero column: rank deficient at 0 */
        nmod_poly_mat_zero(A);
        for (slong i = 0; i < n; i++)
        {
            for (slong j = 1; j < n; j++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, i, j), 0,
                                       n_randint(state, A->modulus));
            nmod_poly_neg(nmod_poly_mat_entry(A, i, i), nmod_poly_mat_entry(A, i, i));
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, i, i), 1, 1);
        }
        return;
    }

    for (slong i = 0; i < m; i++)
        if (n_randint(state, 6) == 0)
            for (slong j = 0; j < n; j++)
                nmod_poly_zero(nmod_poly_mat_entry(A, i, j));
}

TEST_FUNCTION_START(nmod_poly_mat_rank_at_point, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        const int large = n_randint(state, 2);
        const ulong prime = large ? n_randprime(state, 60, 1)
                                  : n_randprime(state, 2 + n_randint(state, 4), 1);
        const slong m = n_randint(state, 8);
        const slong n = n_randint(state, 8);
        const slong len = 1 + n_randint(state, 6);
        const ulong pt = n_randint(state, prime);

        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, m, n, prime);
        _rktest_rand_input(A, len, state);

        const slong rk = (m == 0 || n == 0) ? 0 : nmod_poly_mat_rank(A);
        const slong r0 = nmod_poly_mat_rank_at_zero(A);
        const slong r1 = nmod_poly_mat_rank_at_point(A, 1);
        const slong rs = nmod_poly_mat_rank_at_some_point(A);
        const slong rp = nmod_poly_mat_rank_at_point(A, pt);

        /* consistency with the naive evaluation */
        result = (r0 == _rktest_rank_at(A, 0))
                 && (r1 == _rktest_rank_at(A, 1))
                 && (rs == _rktest_rank_at(A, NMOD_POLY_MAT_SOME_POINT % prime))
                 && (rp == _rktest_rank_at(A, pt));

        /* lower bounds on the rank */
        result = result && r0 <= rk && r1 <= rk && rs <= rk && rp <= rk;

        /* for a large prime, the fixed point gives the rank (the failure
         * probability is at most about min(m,n) * len / 2^59, negligible) */
        if (result && large)
            result = (rs == rk);

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, n = %wd, len = %wd, pt = %wu, "
                               "rank = %wd, at 0: %wd, at 1: %wd, at some point: %wd, at pt: %wd\n",
                               prime, m, n, len, pt, rk, r0, r1, rs, rp);

        nmod_poly_mat_clear(A);
    }

    /* x*I - M with M singular, large prime: deficient at 0, full rank at the
     * fixed point */
    {
        const ulong prime = UWORD(1152921504606846883);
        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, 3, 3, prime);
        for (slong k = 0; k < 3; k++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, k, k), 1, 1);
        /* M = [[0,1,0],[0,0,1],[0,0,0]] */
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 1), 0, prime - 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 1, 2), 0, prime - 1);
        if (nmod_poly_mat_rank_at_zero(A) != 2 || nmod_poly_mat_rank_at_some_point(A) != 3
            || nmod_poly_mat_rank_at_point(A, 1) != 3)
            TEST_FUNCTION_FAIL("x*I - M with M nilpotent\n");
        nmod_poly_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
