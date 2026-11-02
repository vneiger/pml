/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_linsolve.h"

/* checks that den is monic, that u A = den b, and that den and the entries
 * of u have no nontrivial common divisor */
static int _linsolvetest_check(const nmod_poly_mat_t u, const nmod_poly_t den,
                               const nmod_poly_mat_t A, const nmod_poly_mat_t b)
{
    if (nmod_poly_is_zero(den) || den->coeffs[den->length - 1] != 1)
        return 0;

    nmod_poly_mat_t P, Q;
    nmod_poly_mat_init(P, 1, A->c, A->modulus);
    nmod_poly_mat_init(Q, 1, A->c, A->modulus);
    nmod_poly_mat_mul(P, u, A);
    nmod_poly_mat_scalar_mul_nmod_poly(Q, b, den);
    int ok = nmod_poly_mat_equal(P, Q);
    nmod_poly_mat_clear(P);
    nmod_poly_mat_clear(Q);

    nmod_poly_t g;
    nmod_poly_init_mod(g, den->mod);
    nmod_poly_set(g, den);
    for (slong i = 0; i < u->c; i++)
        nmod_poly_gcd(g, g, nmod_poly_mat_entry(u, 0, i));
    ok = ok && nmod_poly_is_one(g);
    nmod_poly_clear(g);

    return ok;
}

/* shifted degree of the solution [u, den] for the shift s = row degree of
 * [A ; b] (zero rows counting as degree 0), that is,
 * max(max_i (deg(u_i) + s_i), deg(den) + s_m) */
static slong _linsolvetest_sdeg(const nmod_poly_mat_t u, const nmod_poly_t den,
                                const nmod_poly_mat_t A, const nmod_poly_mat_t b)
{
    const slong m = A->r;
    slong d = -1;
    for (slong i = 0; i <= m; i++)
    {
        const nmod_poly_struct * c = (i < m) ? nmod_poly_mat_entry(u, 0, i) : den;
        if (nmod_poly_is_zero(c))
            continue;
        slong s = 0;
        for (slong j = 0; j < A->c; j++)
        {
            const nmod_poly_struct * e = (i < m) ? nmod_poly_mat_entry(A, i, j) : nmod_poly_mat_entry(b, 0, j);
            s = FLINT_MAX(s, nmod_poly_degree(e));
        }
        d = FLINT_MAX(d, nmod_poly_degree(c) + s);
    }
    return d;
}

/* example where the kernel basis row giving the solution does not have its
 * shifted pivot on the last index: A = [[x, 1], [x^2 + 1, x]], b = [1, 0],
 * whose solution is u = [-x, 1], den = 1 */
static int _linsolvetest_example(void)
{
    const ulong prime = 101;
    nmod_poly_mat_t A, b, u, v;
    nmod_poly_t den;
    nmod_poly_mat_init(A, 2, 2, prime);
    nmod_poly_mat_init(b, 1, 2, prime);
    nmod_poly_mat_init(u, 1, 2, prime);
    nmod_poly_mat_init(v, 1, 2, prime);
    nmod_poly_init(den, prime);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 0), 1, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 1, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 1, 0), 2, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 1, 1), 1, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(b, 0, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(v, 0, 0), 1, prime - 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(v, 0, 1), 0, 1);

    int ok = nmod_poly_mat_linsolve_via_kernel(u, den, A, b)
             && nmod_poly_mat_equal(u, v) && nmod_poly_is_one(den);

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(u);
    nmod_poly_mat_clear(v);
    nmod_poly_clear(den);
    return ok;
}

/* example over GF(2) where A does not have full row rank and the kernel
 * basis computed for [A ; b] has two rows with nonzero last entry, of
 * shifted degrees 4 and 3: the solution must have shifted degree 3 (the
 * minimum over all solutions); entries are given as bit masks (bit k is the
 * coefficient of degree k) */
static int _linsolvetest_example_mindeg(void)
{
    const ulong mA[4][5] = {{1, 7, 2, 4, 1},
                            {4, 2, 4, 2, 0},
                            {5, 4, 7, 6, 3},
                            {5, 7, 7, 5, 1}};
    const ulong mb[5] = {1, 0, 0, 1, 3};
    nmod_poly_mat_t A, b, u;
    nmod_poly_t den;
    nmod_poly_mat_init(A, 4, 5, 2);
    nmod_poly_mat_init(b, 1, 5, 2);
    nmod_poly_mat_init(u, 1, 4, 2);
    nmod_poly_init(den, 2);
    for (slong j = 0; j < 5; j++)
        for (slong k = 0; k < 3; k++)
        {
            for (slong i = 0; i < 4; i++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, i, j), k, (mA[i][j] >> k) & 1);
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(b, 0, j), k, (mb[j] >> k) & 1);
        }

    int ok = nmod_poly_mat_linsolve_via_kernel(u, den, A, b)
             && _linsolvetest_check(u, den, A, b)
             && (_linsolvetest_sdeg(u, den, A, b) == 3);

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(u);
    nmod_poly_clear(den);
    return ok;
}

TEST_FUNCTION_START(nmod_poly_mat_linsolve_via_kernel, state)
{
    int i, result;

    if (!_linsolvetest_example())
        TEST_FUNCTION_FAIL("fixed example (pivot)\n");
    if (!_linsolvetest_example_mindeg())
        TEST_FUNCTION_FAIL("fixed example (minimal degree)\n");

    for (i = 0; i < 400 * flint_test_multiplier(); i++)
    {
        /* moduli from 2 bits (2 or 3) to 64 bits */
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong kind = n_randint(state, 4);
        /* kind 3: square A */
        const slong n = n_randint(state, 7);
        const slong m = (kind == 3) ? n : (slong) n_randint(state, 7);
        const slong len = 1 + n_randint(state, 6);

        nmod_poly_mat_t A, b, Ab, u, u0, v;
        nmod_poly_t den, den0;
        nmod_poly_mat_init(A, m, n, prime);
        nmod_poly_mat_init(b, 1, n, prime);
        nmod_poly_mat_init(Ab, m + 1, n, prime);
        nmod_poly_mat_init(u, 1, m, prime);
        nmod_poly_mat_init(u0, 1, m, prime);
        nmod_poly_mat_init(v, 1, m, prime);
        nmod_poly_init(den, prime);
        nmod_poly_init(den0, prime);

        /* A: random, or of rank at most r (kinds 1 and 2, sometimes) */
        if (kind != 3 && n_randint(state, 2))
        {
            const slong r = n_randint(state, FLINT_MIN(m, n) + 1);
            nmod_poly_mat_t L, R;
            nmod_poly_mat_init(L, m, r, prime);
            nmod_poly_mat_init(R, r, n, prime);
            nmod_poly_mat_randtest(L, state, len);
            nmod_poly_mat_randtest(R, state, len);
            nmod_poly_mat_mul(A, L, R);
            nmod_poly_mat_clear(L);
            nmod_poly_mat_clear(R);
        }
        else
            nmod_poly_mat_randtest(A, state, len);

        /* b: random, or v A for some random v (kind 2) */
        if (kind == 2)
        {
            nmod_poly_mat_randtest(v, state, len);
            nmod_poly_mat_mul(b, v, A);
        }
        else
            nmod_poly_mat_randtest(b, state, len);

        nmod_poly_mat_concat_vertical(Ab, A, b);
        const slong rkA = nmod_poly_mat_rank(A);
        const int expected = (nmod_poly_mat_rank(Ab) == rkA);

        /* return value, correctness, u and den not modified in case of
         * failure */
        nmod_poly_mat_randtest(u, state, 3);
        nmod_poly_randtest(den, state, 3);
        nmod_poly_mat_set(u0, u);
        nmod_poly_set(den0, den);
        result = (nmod_poly_mat_linsolve_via_kernel(u, den, A, b) == expected);
        if (result && expected)
            result = _linsolvetest_check(u, den, A, b);
        if (result && !expected)
            result = nmod_poly_mat_equal(u, u0) && nmod_poly_equal(den, den0);

        if (result && expected && kind == 2)
        {
            /* minimal shifted degree among all solutions, compared with the
             * known solution [v, 1] */
            nmod_poly_t one;
            nmod_poly_init(one, prime);
            nmod_poly_one(one);
            result = (_linsolvetest_sdeg(u, den, A, b) <= _linsolvetest_sdeg(v, one, A, b));
            /* full row rank: unique solution, which is v (then den = 1) */
            if (result && rkA == m)
                result = nmod_poly_mat_equal(u, v) && nmod_poly_is_one(den);
            nmod_poly_clear(one);
        }

        /* square nonsingular A: compare with FLINT's solver, applied to the
         * transposes (A^T u^T = den b^T) */
        if (result && expected && kind != 2 && m == n && m > 0 && rkA == m)
        {
            nmod_poly_mat_t At, bt, w, P, Q;
            nmod_poly_t d;
            nmod_poly_mat_init(At, n, m, prime);
            nmod_poly_mat_init(bt, n, 1, prime);
            nmod_poly_mat_init(w, m, 1, prime);
            nmod_poly_mat_init(P, m, 1, prime);
            nmod_poly_mat_init(Q, m, 1, prime);
            nmod_poly_init(d, prime);
            nmod_poly_mat_transpose(At, A);
            nmod_poly_mat_transpose(bt, b);
            result = nmod_poly_mat_solve(w, d, At, bt);
            /* w / d == u^T / den */
            nmod_poly_mat_scalar_mul_nmod_poly(P, w, den);
            nmod_poly_mat_transpose(Q, u);
            nmod_poly_mat_scalar_mul_nmod_poly(Q, Q, d);
            result = result && nmod_poly_mat_equal(P, Q);
            nmod_poly_mat_clear(At);
            nmod_poly_mat_clear(bt);
            nmod_poly_mat_clear(w);
            nmod_poly_mat_clear(P);
            nmod_poly_mat_clear(Q);
            nmod_poly_clear(d);
        }

        /* aliasing u and b */
        if (result && n == m)
        {
            nmod_poly_mat_set(u0, b);
            nmod_poly_set(den0, den);
            result = (nmod_poly_mat_linsolve_via_kernel(u0, den0, A, u0) == expected)
                     && nmod_poly_mat_equal(u0, expected ? u : b)
                     && nmod_poly_equal(den0, den);
        }

        /* aliasing u and A */
        if (result && m == 1 && n == 1)
        {
            nmod_poly_mat_set(u0, A);
            nmod_poly_set(den0, den);
            result = (nmod_poly_mat_linsolve_via_kernel(u0, den0, u0, b) == expected)
                     && nmod_poly_mat_equal(u0, expected ? u : A)
                     && nmod_poly_equal(den0, den);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, kind = %wd, m = %wd, n = %wd, len = %wd, rank = %wd\n",
                               prime, kind, m, n, len, rkA);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(b);
        nmod_poly_mat_clear(Ab);
        nmod_poly_mat_clear(u);
        nmod_poly_mat_clear(u0);
        nmod_poly_mat_clear(v);
        nmod_poly_clear(den);
        nmod_poly_clear(den0);
    }

    TEST_FUNCTION_END(state);
}
