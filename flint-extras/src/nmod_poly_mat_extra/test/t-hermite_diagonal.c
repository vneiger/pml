/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_row_basis.h"
#include "hermite_test_utils.h"

/* reference diagonal for orientation orient; returns nonsingularity */
static int _hdtest_reference(nmod_poly_struct * diag, const nmod_poly_mat_t F, orientation_t orient)
{
    if (orient == ROW_LOWER)
        return _htest_lower_hnf_diagonal(diag, F);
    if (orient == ROW_UPPER)
        return _htest_upper_hnf_diagonal(diag, F);

    const slong n = F->r;
    nmod_poly_mat_t T;
    nmod_poly_mat_init(T, n, n, F->modulus);
    nmod_poly_mat_transpose(T, F);
    int res = (orient == COL_UPPER) ? _htest_lower_hnf_diagonal(diag, T)
                                    : _htest_upper_hnf_diagonal(diag, T);
    nmod_poly_mat_clear(T);
    return res;
}

static int _one_test_hermite_diagonal(const nmod_poly_mat_t F, orientation_t orient)
{
    const slong n = F->r;
    nmod_poly_struct * diag = flint_malloc(FLINT_MAX(1, n) * sizeof(nmod_poly_struct));
    nmod_poly_struct * diag_ref = flint_malloc(FLINT_MAX(1, n) * sizeof(nmod_poly_struct));
    for (slong i = 0; i < n; i++)
    {
        nmod_poly_init(diag + i, F->modulus);
        nmod_poly_init(diag_ref + i, F->modulus);
    }

    int nonsing;
    if (orient == ROW_LOWER)
    {
        /* default thresholds; recursion down to dimension 2 with weak Popov
         * row bases; recursion down to dimension 2 with Zhou-Labahn row
         * bases */
        nonsing = nmod_poly_mat_hnf_diagonal_lnz(diag, F);
        for (slong variant = 0; variant < 2 && nonsing >= 0; variant++)
        {
            int nonsing2 = _nmod_poly_mat_hnf_diagonal_lnz(diag_ref, F,
                               (variant == 0) ? NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES : -1, -1);
            for (slong i = 0; nonsing == nonsing2 && nonsing && i < n; i++)
                if (!nmod_poly_equal(diag + i, diag_ref + i))
                    nonsing2 = -1;
            if (nonsing != nonsing2)
            {
                flint_printf("FAIL: different results for variant %wd of the recursion\n", variant);
                nonsing = -1;
            }
        }
    }
    else
        nonsing = nmod_poly_mat_hnf_diagonal(diag, F, orient);
    int nonsing_ref = _hdtest_reference(diag_ref, F, orient);

    int res = (nonsing == nonsing_ref);
    if (!res)
        flint_printf("FAIL: singularity detection (%d instead of %d), orientation %d\n",
                     nonsing, nonsing_ref, (int) orient);

    for (slong i = 0; res && nonsing && i < n; i++)
        if (!nmod_poly_equal(diag + i, diag_ref + i))
        {
            flint_printf("FAIL: wrong diagonal entry %wd, orientation %d\n", i, (int) orient);
            res = 0;
        }

    if (!res)
    {
        flint_printf("input matrix (%wd x %wd, modulus %wu):\n", n, n, F->modulus);
        nmod_poly_mat_print_pretty(F, "x");
        flint_printf("\ncomputed diagonal:\n");
        for (slong i = 0; i < n; i++)
            nmod_poly_print_pretty(diag + i, "x"), flint_printf("\n");
        flint_printf("expected diagonal:\n");
        for (slong i = 0; i < n; i++)
            nmod_poly_print_pretty(diag_ref + i, "x"), flint_printf("\n");
    }

    for (slong i = 0; i < n; i++)
    {
        nmod_poly_clear(diag + i);
        nmod_poly_clear(diag_ref + i);
    }
    flint_free(diag);
    flint_free(diag_ref);
    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_hnf_diagonal, state)
{
    int i, result;

    for (i = 0; i < 600 * flint_test_multiplier(); i++)
    {
        const ulong prime = _htest_rand_prime(state);
        const slong n = n_randint(state, 13);
        const slong len = 1 + n_randint(state, 7);

        nmod_poly_mat_t F;
        nmod_poly_mat_init(F, n, n, prime);
        _htest_rand_square(F, len, state);

        /* main orientation, and one other at random */
        const orientation_t orient = (orientation_t) n_randint(state, 4);
        result = _one_test_hermite_diagonal(F, ROW_LOWER)
                 && _one_test_hermite_diagonal(F, orient);

        nmod_poly_mat_clear(F);

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, n = %wd, len = %wd\n", prime, n, len);
    }

    TEST_FUNCTION_END(state);
}
