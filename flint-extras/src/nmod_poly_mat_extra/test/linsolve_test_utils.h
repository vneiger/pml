/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Shared test-only helpers for the solution of u A = den b over the
 * polynomials (t-linsolve_*.c): check of the output, and generators of the
 * random inputs */

#ifndef LINSOLVE_TEST_UTILS_H
#define LINSOLVE_TEST_UTILS_H

#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

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

/* dimension n of a square system (mostly at most 6, sometimes 8 to 12),
 * and length bound len for the entries of A (shorter for larger n) */
static void _linsolvetest_rand_dims(slong * n, slong * len, flint_rand_t state)
{
    *n = (n_randint(state, 8) == 0) ? 8 + n_randint(state, 5) : n_randint(state, 7);
    *len = 1 + n_randint(state, (*n <= 6) ? 6 : 3);
}

/* A: random square matrix with entries of length less than len; if singular
 * (and A is not empty), of rank less than its dimension n, as a product
 * L R of random matrices of sizes n x r and r x n for a random r < n */
static void _linsolvetest_rand_A(nmod_poly_mat_t A, flint_rand_t state,
                                 slong len, int singular)
{
    const slong n = A->r;
    if (singular && n > 0)
    {
        const slong r = n_randint(state, n);
        nmod_poly_mat_t L, R;
        nmod_poly_mat_init(L, n, r, A->modulus);
        nmod_poly_mat_init(R, r, n, A->modulus);
        nmod_poly_mat_randtest(L, state, len);
        nmod_poly_mat_randtest(R, state, len);
        nmod_poly_mat_mul(A, L, R);
        nmod_poly_mat_clear(L);
        nmod_poly_mat_clear(R);
    }
    else
        nmod_poly_mat_randtest(A, state, len);
}

/* replaces a few rows of A (at most max(n/3, 1), for A of dimension n >= 1)
 * by random rows of length up to (n+1) len + 1: a few rows of much larger
 * degree than the others */
static void _linsolvetest_long_rows(nmod_poly_mat_t A, flint_rand_t state, slong len)
{
    const slong n = A->r;
    const slong nlong = 1 + n_randint(state, FLINT_MAX(n / 3, 1));
    for (slong k = 0; k < nlong; k++)
    {
        const slong ii = n_randint(state, n);
        const slong llen = 1 + n_randint(state, (n + 1) * len + 1);
        for (slong j = 0; j < A->c; j++)
            nmod_poly_randtest(nmod_poly_mat_entry(A, ii, j), state, llen);
    }
}

/* kind of the right-hand side b: 0 random (most often), 1 zero, 2 v A for
 * some random v */
static slong _linsolvetest_rand_kindb(flint_rand_t state)
{
    return (n_randint(state, 8) == 0) ? 1 : (n_randint(state, 4) == 0 ? 2 : 0);
}

/* b of the given kind (see _linsolvetest_rand_kindb), for A of dimension n:
 * for kind 0, random of length len, or shorter, or much longer (up to about
 * 4 (n+1) len, the numerators of the solution having degree about n len),
 * sometimes with a zero entry; for kind 1, b is left as it is (b must be
 * zero); for kind 2, v is set to a random vector of length less than len,
 * and b = v A. Returns the length bound used for b (for messages). */
static slong _linsolvetest_rand_b(nmod_poly_mat_t b, nmod_poly_mat_t v,
                                  const nmod_poly_mat_t A, flint_rand_t state,
                                  slong kindb, slong len)
{
    const slong n = A->r;
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
    return lenb;
}

#endif  /* LINSOLVE_TEST_UTILS_H */
