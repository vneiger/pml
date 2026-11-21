/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra.h"

/* Checks that the rows of A are in row echelon form with `rank` nonzero rows
 * on top: for the "upper" forms (upper = 1) the pivot of a row is the index of
 * its leftmost nonzero entry, for the "lower" forms (upper = 0) that of its
 * rightmost one, and the pivots increase strictly from row to row. If unit is
 * set the pivot entries are 1, and if reduced is set (which implies unit) the
 * pivot columns are zero outside the pivot. Records in *nonunit whether some
 * pivot is not 1 and in *nonreduced whether some pivot column has a nonzero
 * entry outside its pivot. Returns 1 if all holds, 0 otherwise. */
static int
_check_row_echelon(const nmod_mat_t A, slong rank, int upper, int unit,
                   int reduced, int * nonunit, int * nonreduced)
{
    slong i, j, k, prev = -1;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (nmod_mat_entry(A, i, j) >= A->mod.n)
                return 0;

    for (i = rank; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (nmod_mat_entry(A, i, j) != 0)
                return 0;

    for (i = 0; i < rank; i++)
    {
        slong piv = -1;

        if (upper)
        {
            for (j = 0; j < A->c && piv == -1; j++)
                if (nmod_mat_entry(A, i, j) != 0)
                    piv = j;
        }
        else
        {
            for (j = A->c - 1; j >= 0 && piv == -1; j--)
                if (nmod_mat_entry(A, i, j) != 0)
                    piv = j;
        }

        if (piv <= prev)  /* zero row (piv == -1) or pivots not increasing */
            return 0;
        prev = piv;

        if (nmod_mat_entry(A, i, piv) != 1)
        {
            *nonunit = 1;
            if (unit || reduced)
                return 0;
        }

        for (k = 0; k < rank; k++)
            if (k != i && nmod_mat_entry(A, k, piv) != 0)
            {
                *nonreduced = 1;
                if (reduced)
                    return 0;
            }
    }

    return 1;
}

/* Same as _check_row_echelon for the columns of A (pivot of a column: index
 * of its uppermost nonzero entry for the "lower" forms, of its bottommost one
 * for the "upper" forms), through the transpose. */
static int
_check_column_echelon(const nmod_mat_t A, slong rank, int upper, int unit,
                      int reduced, int * nonunit, int * nonreduced)
{
    int res;
    nmod_mat_t T;
    nmod_mat_init(T, A->c, A->r, A->mod.n);
    nmod_mat_transpose(T, A);
    res = _check_row_echelon(T, rank, !upper, unit, reduced, nonunit, nonreduced);
    nmod_mat_clear(T);
    return res;
}

TEST_FUNCTION_START(nmod_mat_rand, state)
{
    slong iter;
    /* per function, counts of the cases where a pivot different from 1, or
     * a nonzero entry in a pivot column outside the pivot, can be expected,
     * and of the cases where one was seen; used to check, over all
     * iterations, that the non-unit and non-reduced forms are not unit or
     * reduced by mistake */
    slong nonunit_cases[10] = {0}, nonunit_seen[10] = {0};
    slong nonreduced_cases[10] = {0}, nonreduced_seen[10] = {0};

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong m, n, rank, fn;
        int unit, full, ok, nonunit = 0, nonreduced = 0;
        ulong p;
        nmod_mat_t A;

        m = n_randint(state, 13);
        n = n_randint(state, 13);
        full = (n_randint(state, 4) == 0);
        rank = full ? FLINT_MIN(m, n) : (slong) n_randint(state, FLINT_MIN(m, n) + 1);
        unit = n_randint(state, 2);

        switch (n_randint(state, 4))
        {
            case 0: p = 2; break;
            case 1: p = n_randprime(state, 2 + n_randint(state, 20), 1); break;
            case 2: p = n_randprime(state, 22 + n_randint(state, 40), 1); break;
            default: p = n_randtest_prime(state, 0); break;
        }

        nmod_mat_init(A, m, n, p);
        nmod_mat_randtest(A, state);  /* entries to overwrite */

        fn = n_randint(state, 10);
        switch (fn)
        {
            case 0:  /* dense, prescribed rank */
                if (full)
                    nmod_mat_rand_fullrank(A, state);
                else
                    nmod_mat_randrank_dense(A, state, rank);
                ok = 1;
                break;
            case 1:
                if (full)
                    nmod_mat_rand_fullrank_lref(A, state, unit);
                else
                    nmod_mat_rand_lref(A, state, rank, unit);
                ok = _check_row_echelon(A, rank, 0, unit, 0, &nonunit, &nonreduced);
                break;
            case 2:
                if (full)
                    nmod_mat_rand_fullrank_lrref(A, state);
                else
                    nmod_mat_rand_lrref(A, state, rank);
                ok = _check_row_echelon(A, rank, 0, 1, 1, &nonunit, &nonreduced);
                break;
            case 3:
                if (full)
                    nmod_mat_rand_fullrank_uref(A, state, unit);
                else
                    nmod_mat_rand_uref(A, state, rank, unit);
                ok = _check_row_echelon(A, rank, 1, unit, 0, &nonunit, &nonreduced);
                break;
            case 4:
                if (full)
                    nmod_mat_rand_fullrank_urref(A, state);
                else
                    nmod_mat_rand_urref(A, state, rank);
                ok = _check_row_echelon(A, rank, 1, 1, 1, &nonunit, &nonreduced);
                break;
            case 5:
                if (full)
                    nmod_mat_rand_fullrank_lcef(A, state, unit);
                else
                    nmod_mat_rand_lcef(A, state, rank, unit);
                ok = _check_column_echelon(A, rank, 0, unit, 0, &nonunit, &nonreduced);
                break;
            case 6:
                if (full)
                    nmod_mat_rand_fullrank_lrcef(A, state);
                else
                    nmod_mat_rand_lrcef(A, state, rank);
                ok = _check_column_echelon(A, rank, 0, 1, 1, &nonunit, &nonreduced);
                break;
            case 7:
                if (full)
                    nmod_mat_rand_fullrank_ucef(A, state, unit);
                else
                    nmod_mat_rand_ucef(A, state, rank, unit);
                ok = _check_column_echelon(A, rank, 1, unit, 0, &nonunit, &nonreduced);
                break;
            case 8:
                if (full)
                    nmod_mat_rand_fullrank_urcef(A, state);
                else
                    nmod_mat_rand_urcef(A, state, rank);
                ok = _check_column_echelon(A, rank, 1, 1, 1, &nonunit, &nonreduced);
                break;
            default:  /* dense, prescribed rank, on a window */
            {
                nmod_mat_t B, W;
                nmod_mat_init(B, m + 2, n + 3, p);
                nmod_mat_randtest(B, state);
                nmod_mat_window_init(W, B, 1, 2, m + 1, n + 2);
                nmod_mat_randrank_dense(W, state, rank);
                nmod_mat_set(A, W);
                nmod_mat_window_clear(W);
                nmod_mat_clear(B);
                ok = 1;
                break;
            }
        }

        if (!ok)
            TEST_FUNCTION_FAIL("function %wd, m = %wd, n = %wd, rank = %wd, unit = %d, "
                               "p = %wu: not of the expected form\n%{nmod_mat}\n",
                               fn, m, n, rank, unit, p, A);

        if (nmod_mat_rank(A) != rank)
            TEST_FUNCTION_FAIL("function %wd, m = %wd, n = %wd, p = %wu: rank %wd "
                               "instead of %wd\n%{nmod_mat}\n",
                               fn, m, n, p, nmod_mat_rank(A), rank, A);

        /* non-unit forms: pivots are random nonzero entries, so with p > 2
           some must differ from 1; non-reduced forms: entries above or below
           the pivots are random, so for rank >= 2 some must be nonzero */
        if ((fn == 1 || fn == 3 || fn == 5 || fn == 7) && !unit && p > 2 && rank > 0)
        {
            nonunit_cases[fn] += 1;
            nonunit_seen[fn] += nonunit;
        }
        if ((fn == 1 || fn == 3 || fn == 5 || fn == 7) && rank >= 2 && p > 2)
        {
            nonreduced_cases[fn] += 1;
            nonreduced_seen[fn] += nonreduced;
        }

        nmod_mat_clear(A);
    }

    /* the probability that 20 such cases show nothing is negligible */
    for (slong fn = 1; fn < 9; fn += 2)
    {
        if (nonunit_cases[fn] >= 20 && nonunit_seen[fn] == 0)
            TEST_FUNCTION_FAIL("function %wd, non-unit echelon forms: all pivots "
                               "equal to 1 in %wd cases\n", fn, nonunit_cases[fn]);
        if (nonreduced_cases[fn] >= 20 && nonreduced_seen[fn] == 0)
            TEST_FUNCTION_FAIL("function %wd, non-reduced echelon forms: all "
                               "reduced in %wd cases\n", fn, nonreduced_cases[fn]);
    }

    TEST_FUNCTION_END(state);
}
