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
static int _lsgeomtest_check(const nmod_poly_mat_t u, const nmod_poly_t den,
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

TEST_FUNCTION_START(nmod_poly_mat_linsolve_geometric, state)
{
    int i, result;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        /* moduli from 2 bits to 64 bits; the small ones do not allow the
         * geometric progression (then via_kernel is used) */
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        /* kind of A: 0 random, 1 one row multiplied by x - 1 (the first
         * point of the progression is 1: A(1) singular, generically of rank
         * n - 1), 2 two rows multiplied by x - 1 (rank of A(1) at most n - 2),
         * 3 singular, 4 a few rows of large degree; kind of b: 0 random, 1
         * zero, 2 v A for some random v */
        const slong kindA = n_randint(state, 5);
        const slong kindb = (n_randint(state, 8) == 0) ? 1 : (n_randint(state, 4) == 0 ? 2 : 0);
        const slong n = (n_randint(state, 8) == 0) ? 8 + n_randint(state, 5) : n_randint(state, 7);
        const slong len = 1 + n_randint(state, (n <= 6) ? 6 : 3);

        nmod_poly_mat_t A, b, u, u0, v, w;
        nmod_poly_t den, den0, f;
        nmod_poly_mat_init(A, n, n, prime);
        nmod_poly_mat_init(b, 1, n, prime);
        nmod_poly_mat_init(u, 1, n, prime);
        nmod_poly_mat_init(u0, 1, n, prime);
        nmod_poly_mat_init(v, 1, n, prime);
        nmod_poly_mat_init(w, 1, n, prime);
        nmod_poly_init(den, prime);
        nmod_poly_init(den0, prime);
        nmod_poly_init(f, prime);

        if (kindA == 3 && n > 0)
        {
            const slong r = n_randint(state, n);
            nmod_poly_mat_t L, R;
            nmod_poly_mat_init(L, n, r, prime);
            nmod_poly_mat_init(R, r, n, prime);
            nmod_poly_mat_randtest(L, state, len);
            nmod_poly_mat_randtest(R, state, len);
            nmod_poly_mat_mul(A, L, R);
            nmod_poly_mat_clear(L);
            nmod_poly_mat_clear(R);
        }
        else
            nmod_poly_mat_randtest(A, state, len);

        if ((kindA == 1 && n > 0) || (kindA == 2 && n > 1))
        {
            /* f = x - 1 */
            nmod_poly_set_coeff_ui(f, 1, 1);
            nmod_poly_set_coeff_ui(f, 0, prime - 1);
            const slong i1 = n_randint(state, n);
            const slong i2 = (kindA == 2) ? (i1 + 1 + (slong) n_randint(state, n - 1)) % n : i1;
            for (slong j = 0; j < n; j++)
            {
                nmod_poly_mul(nmod_poly_mat_entry(A, i1, j), nmod_poly_mat_entry(A, i1, j), f);
                if (kindA == 2)
                    nmod_poly_mul(nmod_poly_mat_entry(A, i2, j), nmod_poly_mat_entry(A, i2, j), f);
            }
        }

        if (kindA == 4 && n > 0)
        {
            const slong nlong = 1 + n_randint(state, FLINT_MAX(n / 3, 1));
            for (slong k = 0; k < nlong; k++)
            {
                const slong ii = n_randint(state, n);
                const slong llen = 1 + n_randint(state, (n + 1) * len + 1);
                for (slong j = 0; j < n; j++)
                    nmod_poly_randtest(nmod_poly_mat_entry(A, ii, j), state, llen);
            }
        }

        /* b: random of length len, or shorter, or much longer, sometimes
         * with zero entries */
        slong lenb = len;
        if (kindb == 0)
        {
            const slong r = n_randint(state, 3);
            lenb = (r == 0) ? len : (r == 1 ? 1 + (slong) n_randint(state, len)
                                            : len * (n + 1) * (1 + (slong) n_randint(state, 4)));
            nmod_poly_mat_randtest(b, state, lenb);
            if (n > 1 && n_randint(state, 4) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(b, 0, n_randint(state, n)));
        }
        else if (kindb == 2)
        {
            nmod_poly_mat_randtest(v, state, len);
            nmod_poly_mat_mul(b, v, A);
        }

        const int nonsingular = (nmod_poly_mat_rank(A) == n);

        /* reference: via_kernel */
        const int expected = nmod_poly_mat_linsolve_via_kernel(w, den0, A, b);

        /* return value, correctness, and u and den not modified in case of
         * failure; if A is nonsingular, the solution is unique */
        nmod_poly_mat_randtest(u, state, 3);
        nmod_poly_randtest(den, state, 3);
        nmod_poly_mat_set(u0, u);
        nmod_poly_set(f, den);
        result = (nmod_poly_mat_linsolve_geometric(u, den, A, b) == expected);
        if (result && expected)
            result = _lsgeomtest_check(u, den, A, b);
        if (result && expected && nonsingular)
            result = nmod_poly_mat_equal(u, w) && nmod_poly_equal(den, den0);
        if (result && expected && nonsingular && kindb == 2)
            result = nmod_poly_mat_equal(u, v) && nmod_poly_is_one(den);
        if (result && !expected)
            result = nmod_poly_mat_equal(u, u0) && nmod_poly_equal(den, f);

        /* aliasing u and b; and chunks of a few points for the evaluation of
         * A (by default, there is only one chunk at these sizes) */
        if (result)
        {
            const slong chunk = n_randint(state, 2) ? 0 : 1 + n_randint(state, 4 * len);
            nmod_poly_mat_set(u0, b);
            result = (_nmod_poly_mat_linsolve_geometric(u0, f, A, u0, chunk) == expected)
                     && nmod_poly_mat_equal(u0, expected ? u : b)
                     && nmod_poly_equal(f, den);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, kindA = %wd, kindb = %wd, n = %wd, len = %wd, lenb = %wd, nonsingular = %d\n",
                               prime, kindA, kindb, n, len, lenb, nonsingular);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(b);
        nmod_poly_mat_clear(u);
        nmod_poly_mat_clear(u0);
        nmod_poly_mat_clear(v);
        nmod_poly_mat_clear(w);
        nmod_poly_clear(den);
        nmod_poly_clear(den0);
        nmod_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
