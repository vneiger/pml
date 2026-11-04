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

/* sum of row degrees and of column degrees, compared with the entries of the
 * row degree and column degree computed entry by entry */
TEST_FUNCTION_START(nmod_poly_mat_degree_sum, state)
{
    for (slong iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong r = n_randint(state, 7);
        const slong c = n_randint(state, 7);
        const slong len = n_randint(state, 6);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, r, c, prime);
        /* sparse: some zero entries, rows and columns */
        nmod_poly_mat_randtest_sparse(mat, state, len, 0.3);

        /* shifts: NULL or random, for rows (length c) and columns (length r) */
        const int null_shift = n_randint(state, 2);
        slong * rshift = flint_malloc((c + 1) * sizeof(slong));
        slong * cshift = flint_malloc((r + 1) * sizeof(slong));
        for (slong j = 0; j < c; j++)
            rshift[j] = (slong) n_randint(state, 11) - 5;
        for (slong i = 0; i < r; i++)
            cshift[i] = (slong) n_randint(state, 11) - 5;

        slong * degs = flint_malloc((FLINT_MAX(r, c) + 1) * sizeof(slong));
        slong sum, ref;

        const ulong nzr = nmod_poly_mat_row_degree_zero(degs, mat, null_shift ? NULL : rshift);
        ref = 0;
        for (slong i = 0; i < r; i++)
            ref += degs[i];
        if (nmod_poly_mat_row_degree_sum(&sum, mat, null_shift ? NULL : rshift) != nzr || sum != ref)
            TEST_FUNCTION_FAIL("row_degree_sum: r = %wd, c = %wd, sum = %wd, expected %wd\n",
                               r, c, sum, ref);

        const ulong nzc = nmod_poly_mat_column_degree_zero(degs, mat, null_shift ? NULL : cshift);
        ref = 0;
        for (slong j = 0; j < c; j++)
            ref += degs[j];
        if (nmod_poly_mat_column_degree_sum(&sum, mat, null_shift ? NULL : cshift) != nzc || sum != ref)
            TEST_FUNCTION_FAIL("column_degree_sum: r = %wd, c = %wd, sum = %wd, expected %wd\n",
                               r, c, sum, ref);

        flint_free(degs);
        flint_free(rshift);
        flint_free(cshift);
        nmod_poly_mat_clear(mat);
    }

    TEST_FUNCTION_END(state);
}
