/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_extra.h"   /* for nmod_find_root */
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_extra/impl.h"

/* random entries of unequal lengths, with zero entries, rows and columns
   sprinkled in: this is what the approximant and interpolant basis
   algorithms feed to the multiplication, and it is what exercises the
   per-entry output lengths of the implementation */
static void _rand_mat(nmod_poly_mat_t M, flint_rand_t state, slong deg,
                      int ragged, int zeros)
{
    slong i, j;

    for (i = 0; i < M->r; i++)
        for (j = 0; j < M->c; j++)
        {
            const slong d = ragged ? (slong) n_randint(state, deg + 1) : deg;
            nmod_poly_randtest(nmod_poly_mat_entry(M, i, j), state, d + 1);
        }

    if (!zeros || M->r == 0 || M->c == 0)
        return;

    for (i = 0; i < M->r; i++)
        for (j = 0; j < M->c; j++)
            if (n_randint(state, 4) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(M, i, j));

    if (n_randint(state, 3) == 0)   /* a zero row */
    {
        i = (slong) n_randint(state, M->r);
        for (j = 0; j < M->c; j++)
            nmod_poly_zero(nmod_poly_mat_entry(M, i, j));
    }
    if (n_randint(state, 3) == 0)   /* a zero column */
    {
        j = (slong) n_randint(state, M->c);
        for (i = 0; i < M->r; i++)
            nmod_poly_zero(nmod_poly_mat_entry(M, i, j));
    }
}

/* membytes == 0 goes through the public entry point; otherwise the
   progression is built here and the bound is forced, which is how the
   row and column grouping gets exercised at testable sizes */
static int _check(flint_rand_t state, ulong prime, slong m, slong k, slong n,
                  slong deg, int ragged, int zeros, ulong membytes, int alias)
{
    nmod_poly_mat_t A, B, C1, C2;
    int res;

    nmod_poly_mat_init(A, m, k, prime);
    nmod_poly_mat_init(B, k, n, prime);
    nmod_poly_mat_init(C1, m, n, prime);
    nmod_poly_mat_init(C2, m, n, prime);

    _rand_mat(A, state, deg, ragged, zeros);
    _rand_mat(B, state, deg, ragged, zeros);
    nmod_poly_mat_randtest(C2, state, 5);   /* stale data in the output */

    nmod_poly_mat_mul(C1, A, B);

    if (membytes == 0)
    {
        if (alias == 1 && n == k)   /* C2 and A have the same shape */
        {
            nmod_poly_mat_set(C2, A);
            nmod_poly_mat_mul_geometric(C2, C2, B);
        }
        else if (alias == 2 && m == k)   /* C2 and B have the same shape */
        {
            nmod_poly_mat_set(C2, B);
            nmod_poly_mat_mul_geometric(C2, A, C2);
        }
        else
            nmod_poly_mat_mul_geometric(C2, A, B);
    }
    else
    {
        const slong len1 = nmod_poly_mat_max_length(A);
        const slong len2 = nmod_poly_mat_max_length(B);

        if (len1 == 0 || len2 == 0)
            nmod_poly_mat_zero(C2);
        else
        {
            const slong len = len1 + len2 - 1;
            nmod_t mod;
            nmod_geometric_progression_t G;
            ulong w;

            nmod_init(&mod, prime);
            w = nmod_find_root(2 * len, mod);
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
            _nmod_geometric_progression_init_function(G, w, len, mod, UWORD(3));
#else
            nmod_geometric_progression_init(G, w, len, mod);
#endif
            _nmod_poly_mat_mul_geometric_precomp_bounded(C2, A, len1, B, len2,
                                                         G, membytes);
            nmod_geometric_progression_clear(G);
        }
    }

    res = nmod_poly_mat_equal(C1, C2);

    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);

    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_mul_geometric, state)
{
    slong i;
    const slong nthreads_save = flint_get_num_threads();

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        /* the field must be large enough for a geometric progression of
           2*len points; the bits below leave ample room */
        const ulong bits = 16 + n_randint(state, 49);
        const ulong prime = n_randprime(state, bits, 1);
        const slong m = 1 + n_randint(state, 20);
        const slong k = 1 + n_randint(state, 20);
        const slong n = 1 + n_randint(state, 20);
        const slong deg = n_randint(state, 40);
        const int ragged = n_randint(state, 2);
        const int zeros = n_randint(state, 2);
        const int alias = (int) n_randint(state, 3);
        ulong membytes = 0;

        flint_set_num_threads(1 + n_randint(state, 4));

        /* one run in three forces the grouping, with a bound small
           enough to cut both the rows of A and the columns of B */
        if (n_randint(state, 3) == 0)
            membytes = 8 * (1 + n_randint(state, 4000));

        if (!_check(state, prime, m, k, n, deg, ragged, zeros, membytes, alias))
            TEST_FUNCTION_FAIL(
                    "m = %wd, k = %wd, n = %wd, deg = %wd, bits = %wu\n"
                    "ragged = %d, zeros = %d, alias = %d, membytes = %wu\n"
                    "threads = %wd\n",
                    m, k, n, deg, bits, ragged, zeros, alias, membytes,
                    flint_get_num_threads());
    }

    /* a few larger, more lopsided shapes: long and thin, wide and short,
       one long entry against many short ones */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 50 + n_randint(state, 14), 1);
        static const slong shapes[6][3] = { {1, 40, 1}, {40, 1, 40}, {2, 2, 2},
                                            {30, 5, 3}, {3, 5, 30}, {17, 17, 17} };
        const slong * s = shapes[i % 6];
        const slong deg = 1 + n_randint(state, 300);
        const ulong membytes = (n_randint(state, 2) == 0)
                               ? 0 : 8 * (1 + n_randint(state, 20000));

        flint_set_num_threads(1 + n_randint(state, 4));

        if (!_check(state, prime, s[0], s[1], s[2], deg, 1, 1, membytes, 0))
            TEST_FUNCTION_FAIL("m = %wd, k = %wd, n = %wd, deg = %wd, "
                               "membytes = %wu, threads = %wd\n",
                               s[0], s[1], s[2], deg, membytes,
                               flint_get_num_threads());
    }

    flint_set_num_threads(nthreads_save);

    TEST_FUNCTION_END(state);
}
