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

#include "nmod_poly_mat_forms.h"

/* shifted degree of entry (i,j), or WORD_MIN if zero */
static slong _ddbtest_weight(const nmod_poly_mat_t A, const slong * shift, int rowwise,
                             slong i, slong j)
{
    const slong d = nmod_poly_degree(nmod_poly_mat_entry(A, i, j));
    if (d < 0)
        return WORD_MIN;
    return d + (shift ? (rowwise ? shift[j] : shift[i]) : 0);
}

/* minimum of the shift entries, 0 if shift is NULL */
static slong _ddbtest_min_shift(const slong * shift, slong len)
{
    if (shift == NULL)
        return 0;
    slong mn = shift[0];
    for (slong i = 1; i < len; i++)
        mn = FLINT_MIN(mn, shift[i]);
    return mn;
}

/* brute force: maximum over injective maps from the k = min(m,n) indices of
 * the smaller dimension to the larger one, by recursion on a */
static slong _ddbtest_brute(const nmod_poly_mat_t A, const slong * shift, int rowwise,
                            slong a, char * used)
{
    const slong m = A->r, n = A->c;
    const int tall = (m >= n);
    const slong k = tall ? n : m, N = tall ? m : n;
    if (a == k)
        return 0;
    slong best = WORD_MIN;
    for (slong b = 0; b < N; b++)
    {
        if (used[b])
            continue;
        const slong wab = tall ? _ddbtest_weight(A, shift, rowwise, b, a)
                               : _ddbtest_weight(A, shift, rowwise, a, b);
        if (wab == WORD_MIN)
            continue;
        used[b] = 1;
        const slong rest = _ddbtest_brute(A, shift, rowwise, a + 1, used);
        used[b] = 0;
        if (rest != WORD_MIN && wab + rest > best)
            best = wab + rest;
    }
    return best;
}

TEST_FUNCTION_START(nmod_poly_mat_degree_det_bound, state)
{
    int i, result;

    /* comparison with brute force, rectangular and square, both
     * orientations, signed shifts, zero entries */
    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        const slong m = n_randint(state, 7);
        const slong n = n_randint(state, 7);
        const int rowwise = n_randint(state, 2);
        const orientation_t orient = rowwise ? ROW_LOWER : COL_UPPER;
        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, m, n, 17);
        for (slong ii = 0; ii < m; ii++)
            for (slong jj = 0; jj < n; jj++)
                if (n_randint(state, 3))
                    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, ii, jj), n_randint(state, 10), 1);
        const slong slen = rowwise ? n : m;
        slong * shift = flint_malloc(FLINT_MAX(1, slen) * sizeof(slong));
        for (slong k = 0; k < slen; k++)
            shift[k] = (slong) n_randint(state, 21) - 10;
        const slong * sh = n_randint(state, 4) ? shift : NULL;

        const slong k = FLINT_MIN(m, n), N = FLINT_MAX(m, n);
        char * used = flint_calloc(FLINT_MAX(1, N), 1);
        slong expected = (k == 0) ? 0 : _ddbtest_brute(A, sh, rowwise, 0, used);
        const int has_map = (expected != WORD_MIN);
        if (!has_map)
            expected = k * _ddbtest_min_shift(sh, slen) - 1;
        slong * assign = flint_malloc(FLINT_MAX(1, k) * sizeof(slong));
        for (slong a = 0; a < k; a++)
            assign[a] = -2;
        const slong res = nmod_poly_mat_degree_det_bound(assign, A, sh, orient);
        result = (res == expected) && (nmod_poly_mat_degree_det_bound(NULL, A, sh, orient) == res);

        /* no map: assign is not modified */
        for (slong a = 0; result && !has_map && a < k; a++)
            result = (assign[a] == -2);

        /* the returned map is injective, avoids zero entries, and reaches
         * the bound */
        if (result && has_map)
        {
            for (slong b = 0; b < N; b++)
                used[b] = 0;
            slong sum = 0;
            for (slong a = 0; result && a < k; a++)
            {
                const slong b = assign[a];
                if (b < 0 || b >= N || used[b])
                {
                    result = 0;
                    break;
                }
                used[b] = 1;
                const slong wab = (m >= n) ? _ddbtest_weight(A, sh, rowwise, b, a)
                                           : _ddbtest_weight(A, sh, rowwise, a, b);
                if (wab == WORD_MIN)
                    result = 0;
                else
                    sum += wab;
            }
            result = result && (sum == res);
        }

        if (!result)
            TEST_FUNCTION_FAIL("m = %wd, n = %wd, rowwise = %d, result %wd, expected %wd\n",
                               m, n, rowwise, res, expected);

        flint_free(assign);
        flint_free(used);
        flint_free(shift);
        nmod_poly_mat_clear(A);
    }

    /* shift entries close to the largest allowed absolute value
     * WORD_MAX / (16 k), comparison with brute force */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        const slong m = 1 + n_randint(state, 6);
        const slong n = 1 + n_randint(state, 6);
        const int rowwise = n_randint(state, 2);
        const orientation_t orient = rowwise ? ROW_UPPER : COL_LOWER;
        const slong k = FLINT_MIN(m, n), N = FLINT_MAX(m, n);
        const slong lim = WORD_MAX / (16 * k);
        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, m, n, 17);
        for (slong ii = 0; ii < m; ii++)
            for (slong jj = 0; jj < n; jj++)
                if (n_randint(state, 3))
                    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, ii, jj), n_randint(state, 10), 1);
        const slong slen = rowwise ? n : m;
        slong * shift = flint_malloc(slen * sizeof(slong));
        for (slong kk = 0; kk < slen; kk++)
        {
            const slong r = n_randint(state, 4);
            shift[kk] = (r == 0) ? lim : (r == 1) ? -lim
                        : (n_randint(state, 2) ? 1 : -1) * (lim - (slong) n_randint(state, 1000));
        }

        char * used = flint_calloc(N, 1);
        slong expected = _ddbtest_brute(A, shift, rowwise, 0, used);
        if (expected == WORD_MIN)
            expected = k * _ddbtest_min_shift(shift, slen) - 1;
        const slong res = nmod_poly_mat_degree_det_bound(NULL, A, shift, orient);
        result = (res == expected);
        if (!result)
            TEST_FUNCTION_FAIL("large shifts, m = %wd, n = %wd, rowwise = %d, result %wd, expected %wd\n",
                               m, n, rowwise, res, expected);

        flint_free(used);
        flint_free(shift);
        nmod_poly_mat_clear(A);
    }

    /* square matrices with random coefficients over a large field: the bound
     * is reached by the determinant (failure probability negligible); if
     * there is no assignment, the matrix is singular and the bound is -1 */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        const slong n = 1 + n_randint(state, 6);
        const ulong prime = n_randprime(state, 60, 1);
        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, n, n, prime);
        for (slong ii = 0; ii < n; ii++)
            for (slong jj = 0; jj < n; jj++)
                if (n_randint(state, 3))
                {
                    const slong len = 1 + n_randint(state, 8);
                    nmod_poly_randtest(nmod_poly_mat_entry(A, ii, jj), state, len);
                    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, ii, jj), len - 1,
                                           1 + n_randint(state, prime - 1));
                }
        nmod_poly_t det;
        nmod_poly_init(det, prime);
        nmod_poly_mat_det(det, A);
        const slong bound = nmod_poly_mat_degree_det_bound(NULL, A, NULL, ROW_LOWER);
        result = (nmod_poly_degree(det) == bound);
        nmod_poly_clear(det);
        nmod_poly_mat_clear(A);
        if (!result)
            TEST_FUNCTION_FAIL("determinant degree, n = %wd, bound = %wd\n", n, bound);
    }

    TEST_FUNCTION_END(state);
}
