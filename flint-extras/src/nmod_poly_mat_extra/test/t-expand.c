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

#include "nmod_poly_mat_utils.h"

/* x^d-adic expansion by blocks of columns (rows == 0) or of rows (rows == 1),
 * with random orders and d:
 *   - expanding then compressing a polynomial matrix gives it back, with
 *     column b (resp. row a) truncated at order orders[b] (resp. orders[a]);
 *   - compressing then expanding a matrix whose pieces have the degrees of an
 *     expansion gives it back;
 *   - compressing A times the expansion (resp. the expansion times B) gives
 *     A times the matrix (resp. the matrix times B): the overlapping parts of
 *     the pieces are added.
 * The outputs hold random data beforehand. */
static int
_expand_round_trips(slong m, slong n, ulong prime, int rows, flint_rand_t state)
{
    const slong norders = rows ? m : n;
    const slong d = 1 + n_randint(state, 10);
    const slong l = n_randint(state, 6);
    slong * orders = flint_malloc(FLINT_MAX(1, norders) * sizeof(slong));
    slong E = 0, off;
    int result;

    for (slong b = 0; b < norders; b++)
    {
        orders[b] = n_randint(state, 4) ? n_randint(state, 30) : 0;
        E += (orders[b] + d - 1) / d;
    }

    nmod_poly_mat_t F, T, G, X, Y, A, P, C, R;
    nmod_poly_mat_init(F, m, n, prime);
    nmod_poly_mat_init(T, m, n, prime);
    nmod_poly_mat_init(G, m, n, prime);
    nmod_poly_mat_init(X, rows ? E : m, rows ? n : E, prime);
    nmod_poly_mat_init(Y, rows ? E : m, rows ? n : E, prime);

    /* polynomial -> expanded -> polynomial; T is F truncated at the orders */
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        {
            const slong o = orders[rows ? i : j];
            nmod_poly_randtest(nmod_poly_mat_entry(F, i, j), state, n_randint(state, o + 5));
            nmod_poly_set_trunc(nmod_poly_mat_entry(T, i, j), nmod_poly_mat_entry(F, i, j), o);
        }
    nmod_poly_mat_randtest(X, state, n_randint(state, 12));
    nmod_poly_mat_randtest(G, state, n_randint(state, 40));

    if (rows)
    {
        nmod_poly_mat_expand_rows(X, F, orders, d);
        nmod_poly_mat_compress_rows(G, X, orders, d);
    }
    else
    {
        nmod_poly_mat_expand_columns(X, F, orders, d);
        nmod_poly_mat_compress_columns(G, X, orders, d);
    }
    result = nmod_poly_mat_equal(G, T);

    /* expanded -> polynomial -> expanded: piece k of block b has degree less
     * than min(d, orders[b] - k d) */
    off = 0;
    for (slong b = 0; b < norders; b++)
        for (slong k = 0; k * d < orders[b]; k++, off++)
            for (slong t = 0; t < (rows ? n : m); t++)
                nmod_poly_randtest(rows ? nmod_poly_mat_entry(X, off, t)
                                        : nmod_poly_mat_entry(X, t, off),
                                   state, n_randint(state, FLINT_MIN(d, orders[b] - k * d) + 1));
    nmod_poly_mat_randtest(Y, state, n_randint(state, 12));
    nmod_poly_mat_randtest(G, state, n_randint(state, 40));

    if (rows)
    {
        nmod_poly_mat_compress_rows(G, X, orders, d);
        nmod_poly_mat_expand_rows(Y, G, orders, d);
    }
    else
    {
        nmod_poly_mat_compress_columns(G, X, orders, d);
        nmod_poly_mat_expand_columns(Y, G, orders, d);
    }
    result = result && nmod_poly_mat_equal(Y, X);

    /* products: A (l x m) times the expansion of T by columns, or the
     * expansion of T by rows times A (n x l) */
    nmod_poly_mat_init(A, rows ? n : l, rows ? l : m, prime);
    nmod_poly_mat_randtest(A, state, n_randint(state, 7));
    nmod_poly_mat_init(P, rows ? E : l, rows ? l : E, prime);
    nmod_poly_mat_init(C, rows ? m : l, rows ? l : n, prime);
    nmod_poly_mat_init(R, rows ? m : l, rows ? l : n, prime);
    nmod_poly_mat_randtest(C, state, n_randint(state, 40));

    if (rows)
    {
        nmod_poly_mat_expand_rows(X, T, orders, d);
        nmod_poly_mat_mul(P, X, A);
        nmod_poly_mat_compress_rows(C, P, orders, d);
        nmod_poly_mat_mul(R, T, A);
    }
    else
    {
        nmod_poly_mat_expand_columns(X, T, orders, d);
        nmod_poly_mat_mul(P, A, X);
        nmod_poly_mat_compress_columns(C, P, orders, d);
        nmod_poly_mat_mul(R, A, T);
    }
    result = result && nmod_poly_mat_equal(C, R);

    nmod_poly_mat_clear(R);
    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(P);
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(Y);
    nmod_poly_mat_clear(X);
    nmod_poly_mat_clear(G);
    nmod_poly_mat_clear(T);
    nmod_poly_mat_clear(F);
    flint_free(orders);

    return result;
}

TEST_FUNCTION_START(nmod_poly_mat_expand, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong m = n_randint(state, 8);
        const slong n = n_randint(state, 8);

        for (int rows = 0; rows < 2; rows++)
            if (!_expand_round_trips(m, n, prime, rows, state))
                TEST_FUNCTION_FAIL("%s, prime = %wu, m = %wd, n = %wd\n",
                                   rows ? "rows" : "columns", prime, m, n);
    }

    /* cancellation in the compression: x^d - x^d is the zero polynomial */
    for (int rows = 0; rows < 2; rows++)
    {
        const ulong prime = 7;
        const slong d = 3, orders[1] = {2 * d};
        nmod_poly_mat_t X, G;
        nmod_poly_mat_init(X, rows ? 2 : 1, rows ? 1 : 2, prime);
        nmod_poly_mat_init(G, 1, 1, prime);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(X, 0, 0), d, 1);
        nmod_poly_set_coeff_ui(rows ? nmod_poly_mat_entry(X, 1, 0) : nmod_poly_mat_entry(X, 0, 1), 0, prime - 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(G, 0, 0), 5, 1);
        if (rows)
            nmod_poly_mat_compress_rows(G, X, orders, d);
        else
            nmod_poly_mat_compress_columns(G, X, orders, d);
        if (!nmod_poly_mat_is_zero(G) || nmod_poly_mat_entry(G, 0, 0)->length != 0)
            TEST_FUNCTION_FAIL("%s: cancellation, the output is not normalised\n",
                               rows ? "rows" : "columns");
        nmod_poly_mat_clear(G);
        nmod_poly_mat_clear(X);
    }

    TEST_FUNCTION_END(state);
}
