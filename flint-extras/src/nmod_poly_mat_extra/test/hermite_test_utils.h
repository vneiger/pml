/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Shared test-only helpers for row bases, Hermite diagonals and
 * determinants: reference computations based on the (upper, row-wise)
 * Hermite normal form algorithms, and generators of structured random
 * matrices (nontrivial Smith forms, triangular, characteristic matrices of
 * derogatory constant matrices, unimodular, ...) */

#ifndef HERMITE_TEST_UTILS_H
#define HERMITE_TEST_UTILS_H

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_utils.h"

/*------------------------------------------------------------*/
/* reference computations                                     */
/*------------------------------------------------------------*/

/* hnf <- upper row-wise Hermite normal form of mat (same dimensions,
 * zero rows at the bottom), returns rank of mat */
static slong _htest_upper_hnf(nmod_poly_mat_t hnf, const nmod_poly_mat_t mat)
{
    nmod_poly_mat_set(hnf, mat);
    if (mat->r == 0 || mat->c == 0)
        return 0;
    slong * pivind = flint_malloc(FLINT_MAX(1, FLINT_MIN(mat->r, mat->c)) * sizeof(slong));
    slong * mrp = flint_malloc(FLINT_MAX(1, mat->r) * sizeof(slong));
    slong rk = nmod_poly_mat_hnf_ur_revlex_xgcd(hnf, NULL, pivind, mrp);
    flint_free(pivind);
    flint_free(mrp);
    return rk;
}

/* returns 1 if the first r1 rows of A1 and the first r2 rows of A2 generate
 * the same K[x]-module (A1 and A2 must have the same number of columns) */
static int _htest_same_row_space(const nmod_poly_mat_t A1, slong r1,
                                 const nmod_poly_mat_t A2, slong r2)
{
    if (A1->c != A2->c)
        return 0;
    nmod_poly_mat_t B1, B2, H1, H2;
    nmod_poly_mat_window_init(B1, A1, 0, 0, r1, A1->c);
    nmod_poly_mat_window_init(B2, A2, 0, 0, r2, A2->c);
    nmod_poly_mat_init(H1, r1, A1->c, A1->modulus);
    nmod_poly_mat_init(H2, r2, A2->c, A2->modulus);
    slong rk1 = _htest_upper_hnf(H1, B1);
    slong rk2 = _htest_upper_hnf(H2, B2);
    int res = (rk1 == rk2);
    for (slong i = 0; res && i < rk1; i++)
        for (slong j = 0; res && j < A1->c; j++)
            if (!nmod_poly_equal(nmod_poly_mat_entry(H1, i, j), nmod_poly_mat_entry(H2, i, j)))
                res = 0;
    nmod_poly_mat_clear(H1);
    nmod_poly_mat_clear(H2);
    nmod_poly_mat_window_clear(B1);
    nmod_poly_mat_window_clear(B2);
    return res;
}

/* for square mat, computes the diagonal of its upper row-wise (ROW_UPPER)
 * Hermite normal form; returns 0 if mat is singular (then diag is
 * undefined), 1 otherwise */
static int _htest_upper_hnf_diagonal(nmod_poly_struct * diag, const nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    if (n == 0)
        return 1;
    nmod_poly_mat_t hnf;
    nmod_poly_mat_init(hnf, n, n, mat->modulus);
    slong rk = _htest_upper_hnf(hnf, mat);
    if (rk == n)
        for (slong i = 0; i < n; i++)
            nmod_poly_set(diag + i, nmod_poly_mat_entry(hnf, i, i));
    nmod_poly_mat_clear(hnf);
    return (rk == n);
}

/* for square mat, computes the diagonal of its lower row-wise (ROW_LOWER)
 * Hermite normal form, using the upper row-wise algorithm on the matrix with
 * reversed rows and columns; returns 0 if mat is singular (then diag is
 * undefined), 1 otherwise */
static int _htest_lower_hnf_diagonal(nmod_poly_struct * diag, const nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    if (n == 0)
        return 1;
    nmod_poly_mat_t rev;
    nmod_poly_mat_init(rev, n, n, mat->modulus);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(rev, n-1-i, n-1-j), nmod_poly_mat_entry(mat, i, j));
    int res = _htest_upper_hnf_diagonal(diag, rev);
    if (res)
        for (slong i = 0; i < n/2; i++)
            nmod_poly_swap(diag + i, diag + n - 1 - i);
    nmod_poly_mat_clear(rev);
    return res;
}

/*------------------------------------------------------------*/
/* generators                                                 */
/*------------------------------------------------------------*/

/* random small prime, with high probability of 2, 3, 5 */
static ulong _htest_rand_prime(flint_rand_t state)
{
    const ulong c = n_randint(state, 10);
    if (c == 0) return 2;
    if (c == 1) return 3;
    if (c == 2) return 5;
    return n_randprime(state, 2 + n_randint(state, 63), 1);
}

/* random unimodular matrix (lower unitriangular times upper unitriangular,
 * random off-diagonal entries of length < len, rows then shuffled) */
static void _htest_rand_unimodular(nmod_poly_mat_t U, slong len, flint_rand_t state)
{
    const slong n = U->r;
    nmod_poly_mat_t L, V;
    nmod_poly_mat_init(L, n, n, U->modulus);
    nmod_poly_mat_init(V, n, n, U->modulus);
    for (slong i = 0; i < n; i++)
    {
        nmod_poly_one(nmod_poly_mat_entry(L, i, i));
        nmod_poly_one(nmod_poly_mat_entry(V, i, i));
        for (slong j = 0; j < i; j++)
        {
            nmod_poly_randtest(nmod_poly_mat_entry(L, i, j), state, len);
            nmod_poly_randtest(nmod_poly_mat_entry(V, j, i), state, len);
        }
    }
    nmod_poly_mat_mul(U, L, V);
    slong * perm = flint_malloc(FLINT_MAX(1, n) * sizeof(slong));
    _perm_randtest(perm, n, state);
    nmod_poly_mat_permute_rows(U, perm, NULL);
    flint_free(perm);
    nmod_poly_mat_clear(L);
    nmod_poly_mat_clear(V);
}

/* random monic polynomial of degree < len (degree 0 allowed) */
static void _htest_rand_monic(nmod_poly_t f, slong len, flint_rand_t state)
{
    slong d = (len <= 1) ? 0 : (slong) n_randint(state, len);
    nmod_poly_randtest_monic(f, state, d + 1);
}

/* random n x n matrix with a nontrivial Smith form: U * D * V with U, V
 * random (unimodular or plain random) and D diagonal with a divisibility
 * chain of random monic polynomials */
static void _htest_rand_nontrivial_smith(nmod_poly_mat_t F, slong len, flint_rand_t state)
{
    const slong n = F->r;
    nmod_poly_mat_t U, D, V;
    nmod_poly_mat_init(U, n, n, F->modulus);
    nmod_poly_mat_init(D, n, n, F->modulus);
    nmod_poly_mat_init(V, n, n, F->modulus);
    const slong lenuv = FLINT_MAX(1, len/2);
    if (n_randint(state, 2))
        _htest_rand_unimodular(U, lenuv, state);
    else
        nmod_poly_mat_rand(U, state, lenuv);
    if (n_randint(state, 2))
        _htest_rand_unimodular(V, lenuv, state);
    else
        nmod_poly_mat_rand(V, state, lenuv);
    nmod_poly_t f;
    nmod_poly_init(f, F->modulus);
    nmod_poly_one(f);
    for (slong i = 0; i < n; i++)
    {
        nmod_poly_t g;
        nmod_poly_init(g, F->modulus);
        _htest_rand_monic(g, (n_randint(state, 3) == 0) ? 3 : 1, state);
        nmod_poly_mul(f, f, g);
        nmod_poly_set(nmod_poly_mat_entry(D, i, i), f);
        nmod_poly_clear(g);
    }
    nmod_poly_clear(f);
    nmod_poly_mat_mul(F, U, D);
    nmod_poly_mat_mul(F, F, V);
    nmod_poly_mat_clear(U);
    nmod_poly_mat_clear(D);
    nmod_poly_mat_clear(V);
}

/* random lower (lower != 0) or upper triangular matrix with nonzero diagonal */
static void _htest_rand_triangular(nmod_poly_mat_t F, slong len, int lower, flint_rand_t state)
{
    const slong n = F->r;
    nmod_poly_mat_zero(F);
    for (slong i = 0; i < n; i++)
    {
        _htest_rand_monic(nmod_poly_mat_entry(F, i, i), FLINT_MAX(1, len), state);
        for (slong j = 0; j < i; j++)
        {
            if (lower)
                nmod_poly_randtest(nmod_poly_mat_entry(F, i, j), state, len);
            else
                nmod_poly_randtest(nmod_poly_mat_entry(F, j, i), state, len);
        }
    }
}

/* x*I - M where M is similar to a block diagonal matrix with several copies
 * of a same random block (derogatory: nontrivial invariant factors) */
static void _htest_rand_charmat_derogatory(nmod_poly_mat_t F, flint_rand_t state)
{
    const slong n = F->r;
    const ulong p = F->modulus;
    nmod_mat_t M, B, P, Pinv, T;
    nmod_mat_init(M, n, n, p);
    nmod_mat_zero(M);
    slong b = 1 + n_randint(state, FLINT_MAX(1, n/2));
    nmod_mat_init(B, b, b, p);
    nmod_mat_randtest(B, state);
    for (slong i0 = 0; i0 + b <= n; i0 += b)
        for (slong i = 0; i < b; i++)
            for (slong j = 0; j < b; j++)
                nmod_mat_entry(M, i0 + i, i0 + j) = nmod_mat_entry(B, i, j);
    for (slong i = (n/b)*b; i < n; i++)
        nmod_mat_entry(M, i, i) = n_randint(state, p);
    /* random similarity transformation */
    nmod_mat_init(P, n, n, p);
    nmod_mat_init(Pinv, n, n, p);
    nmod_mat_init(T, n, n, p);
    do
        nmod_mat_randtest(P, state);
    while (!nmod_mat_inv(Pinv, P));
    nmod_mat_mul(T, P, M);
    nmod_mat_mul(M, T, Pinv);
    /* F = x*I - M */
    nmod_poly_mat_zero(F);
    for (slong i = 0; i < n; i++)
    {
        for (slong j = 0; j < n; j++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, j), 0,
                                   nmod_neg(nmod_mat_entry(M, i, j), M->mod));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, i), 1, 1);
    }
    nmod_mat_clear(M);
    nmod_mat_clear(B);
    nmod_mat_clear(P);
    nmod_mat_clear(Pinv);
    nmod_mat_clear(T);
}

/* multiplies column j of F by x^k */
static void _htest_shift_column(nmod_poly_mat_t F, slong j, slong k)
{
    for (slong i = 0; i < F->r; i++)
        nmod_poly_shift_left(nmod_poly_mat_entry(F, i, j), nmod_poly_mat_entry(F, i, j), k);
}

/* random square matrix, picked among many structures; most of them are
 * nonsingular, a few may be singular */
static void _htest_rand_square(nmod_poly_mat_t F, slong len, flint_rand_t state)
{
    const slong n = F->r;
    const ulong c = n_randint(state, 16);
    if (n == 0)
        return;
    switch (c)
    {
        case 0: case 1:
            nmod_poly_mat_rand(F, state, len);
            break;
        case 2:
            nmod_poly_mat_randtest(F, state, len);
            break;
        case 3:
        {
            slong * rlen = flint_malloc(n * sizeof(slong));
            for (slong i = 0; i < n; i++)
                rlen[i] = n_randtest(state) % (len + 1);
            nmod_poly_mat_rand_row_degree(F, state, rlen);
            flint_free(rlen);
            break;
        }
        case 4:
        {
            slong * clen = flint_malloc(n * sizeof(slong));
            for (slong j = 0; j < n; j++)
                clen[j] = n_randtest(state) % (len + 1);
            nmod_poly_mat_rand_column_degree(F, state, clen);
            flint_free(clen);
            break;
        }
        case 5: case 6: case 7:
            _htest_rand_nontrivial_smith(F, len, state);
            break;
        case 8:
            _htest_rand_triangular(F, len, 1, state);
            break;
        case 9:
            _htest_rand_triangular(F, len, 0, state);
            break;
        case 10:
            _htest_rand_charmat_derogatory(F, state);
            break;
        case 11:
            _htest_rand_unimodular(F, FLINT_MAX(1, len/2), state);
            break;
        case 12:
            nmod_poly_mat_rand(F, state, 1);
            break;
        case 13:
            nmod_poly_mat_rand(F, state, len);
            _htest_shift_column(F, n_randint(state, n), 1 + n_randint(state, 3));
            break;
        case 14:
            nmod_poly_mat_randtest_sparse(F, state, len, 0.3);
            break;
        default:
            /* product of a random matrix and a matrix with nontrivial Smith form */
        {
            nmod_poly_mat_t G;
            nmod_poly_mat_init(G, n, n, F->modulus);
            _htest_rand_nontrivial_smith(G, FLINT_MAX(1, len/2), state);
            nmod_poly_mat_rand(F, state, FLINT_MAX(1, len/2));
            nmod_poly_mat_mul(F, F, G);
            nmod_poly_mat_clear(G);
        }
    }
}

/* random m x n matrix (any shape), picked among several structures */
static void _htest_rand_rect(nmod_poly_mat_t A, slong len, flint_rand_t state)
{
    const slong m = A->r;
    const slong n = A->c;
    if (m == 0 || n == 0)
        return;
    const ulong c = n_randint(state, 8);
    if (c <= 1)
        nmod_poly_mat_rand(A, state, len);
    else if (c == 2)
        nmod_poly_mat_randtest(A, state, len);
    else if (c == 3)
        nmod_poly_mat_randtest_sparse(A, state, len, 0.2);
    else
    {
        /* A = W * B, with W random m x k and B random k x n structured, for
         * some k <= min(m, n) (rank <= k) */
        slong k = (c == 4) ? FLINT_MIN(m, n) : 1 + n_randint(state, FLINT_MIN(m, n));
        nmod_poly_mat_t W, B, Bsq;
        nmod_poly_mat_init(W, m, k, A->modulus);
        nmod_poly_mat_init(B, k, n, A->modulus);
        nmod_poly_mat_init(Bsq, n, n, A->modulus);
        nmod_poly_mat_rand(W, state, FLINT_MAX(1, len/2));
        if (c == 7)  /* constant left factor */
            nmod_poly_mat_rand(W, state, 1);
        _htest_rand_square(Bsq, FLINT_MAX(1, len/2), state);
        for (slong i = 0; i < k; i++)
            for (slong j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(B, i, j), nmod_poly_mat_entry(Bsq, i, j));
        nmod_poly_mat_mul(A, W, B);
        nmod_poly_mat_clear(W);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(Bsq);
    }
}

#endif /* HERMITE_TEST_UTILS_H */
