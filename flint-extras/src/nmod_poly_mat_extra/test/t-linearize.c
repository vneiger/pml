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

/* Round trips through the linearization by blocks of columns (rows == 0) or
 * of rows (rows == 1): linearizing then polynomializing a polynomial matrix
 * whose column b (resp. row a) has degree less than orders[b] (resp.
 * orders[a]) gives it back, and polynomializing then linearizing a constant
 * matrix gives it back. The outputs hold random data beforehand. */
static int
_linearize_round_trips(slong m, slong n, ulong prime, int rows, flint_rand_t state)
{
    const slong norders = rows ? m : n;
    slong * orders = flint_malloc(FLINT_MAX(1, norders) * sizeof(slong));
    slong D = 0;
    int result;

    for (slong b = 0; b < norders; b++)
    {
        orders[b] = n_randint(state, 4) ? n_randint(state, 20) : 0;
        D += orders[b];
    }

    nmod_poly_mat_t F, G;
    nmod_mat_t C, C2;
    nmod_poly_mat_init(F, m, n, prime);
    nmod_poly_mat_init(G, m, n, prime);
    nmod_mat_init(C, rows ? D : m, rows ? n : D, prime);
    nmod_mat_init(C2, rows ? D : m, rows ? n : D, prime);

    /* polynomial -> constant -> polynomial */
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
            nmod_poly_randtest(nmod_poly_mat_entry(F, i, j), state,
                               n_randint(state, orders[rows ? i : j] + 1));
    nmod_mat_randtest(C, state);
    nmod_poly_mat_randtest(G, state, n_randint(state, 24));

    if (rows)
    {
        nmod_poly_mat_linearize_rows(C, F, orders);
        nmod_poly_mat_polynomialize_rows(G, C, orders);
    }
    else
    {
        nmod_poly_mat_linearize_columns(C, F, orders);
        nmod_poly_mat_polynomialize_columns(G, C, orders);
    }
    result = nmod_poly_mat_equal(G, F);

    /* constant -> polynomial -> constant */
    nmod_mat_randtest(C, state);
    nmod_mat_randtest(C2, state);
    nmod_poly_mat_randtest(G, state, n_randint(state, 24));

    if (rows)
    {
        nmod_poly_mat_polynomialize_rows(G, C, orders);
        nmod_poly_mat_linearize_rows(C2, G, orders);
    }
    else
    {
        nmod_poly_mat_polynomialize_columns(G, C, orders);
        nmod_poly_mat_linearize_columns(C2, G, orders);
    }
    result = result && nmod_mat_equal(C2, C);

    nmod_mat_clear(C2);
    nmod_mat_clear(C);
    nmod_poly_mat_clear(G);
    nmod_poly_mat_clear(F);
    flint_free(orders);

    return result;
}

TEST_FUNCTION_START(nmod_poly_mat_linearize, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong m = n_randint(state, 12);
        const slong n = n_randint(state, 12);

        for (int rows = 0; rows < 2; rows++)
            if (!_linearize_round_trips(m, n, prime, rows, state))
                TEST_FUNCTION_FAIL("%s, prime = %wu, m = %wd, n = %wd\n",
                                   rows ? "rows" : "columns", prime, m, n);
    }

    TEST_FUNCTION_END(state);
}
