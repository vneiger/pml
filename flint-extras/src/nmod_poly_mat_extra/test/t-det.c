/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_determinant.h"
#include "hermite_test_utils.h"

/* test one given input, comparing with FLINT's nmod_poly_mat_det */
static int core_test_determinant(const nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    int res = 1;

    nmod_poly_t det, det_correct;
    nmod_poly_init(det, mat->modulus);
    nmod_poly_init(det_correct, mat->modulus);
    nmod_poly_mat_det(det_correct, mat);

    /* Mulders and Storjohann's algorithm, row by row variant */
    {
        nmod_poly_mat_t copy_mat;
        nmod_poly_mat_init_set(copy_mat, mat);
        nmod_poly_mat_det_iter(det, copy_mat);
        nmod_poly_mat_clear(copy_mat);
        if (!nmod_poly_equal(det_correct, det))
        {
            flint_printf("FAIL: det_iter\n");
            res = 0;
        }
    }

    /* expansion by minors */
    if (res && n <= 9)
    {
        nmod_poly_mat_det_minors(det, mat);
        if (!nmod_poly_equal(det_correct, det))
        {
            flint_printf("FAIL: det_minors\n");
            res = 0;
        }
    }

    /* Labahn-Neiger-Zhou */
    if (res)
    {
        nmod_poly_mat_det_lnz(det, mat);
        if (!nmod_poly_equal(det_correct, det))
        {
            flint_printf("FAIL: det_lnz\n");
            res = 0;
        }
    }

    /* general interface */
    if (res)
    {
        nmod_poly_mat_determinant(det, mat);
        if (!nmod_poly_equal(det_correct, det))
        {
            flint_printf("FAIL: determinant (general interface)\n");
            res = 0;
        }
    }

    /* evaluation at geometric progression */
    if (res)
    {
        nmod_poly_one(det);
        const int ok = nmod_poly_mat_det_geometric(det, mat);
        if (ok && !nmod_poly_equal(det_correct, det))
        {
            flint_printf("FAIL: det_geometric\n");
            res = 0;
        }
        if (!ok && (!nmod_poly_is_one(det)
                    || NMOD_POLY_CAN_USE_GEOMETRIC(mat->modulus, n * (nmod_poly_mat_degree(mat) + 1) + 1)))
        {
            flint_printf("FAIL: det_geometric refused a large enough modulus, or modified det\n");
            res = 0;
        }
    }

    if (!res)
    {
        flint_printf("input matrix (%wd x %wd, modulus %wu):\n", n, n, mat->modulus);
        nmod_poly_mat_print_pretty(mat, "x");
        flint_printf("\n");
    }

    nmod_poly_clear(det_correct);
    nmod_poly_clear(det);
    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_det, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 62);
        ulong rdim = 1 + n_randint(state, 20);
        ulong order = n_randint(state, 40);

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, rdim, rdim, prime);

        if (i < 30)
            nmod_poly_mat_randtest_sparse(mat, state, order, 0.2);
        else if (i < 60)
            nmod_poly_mat_rand(mat, state, order);
        else
            nmod_poly_mat_randtest(mat, state, order);

        result = core_test_determinant(mat);

        nmod_poly_mat_clear(mat);

        if (!result)
            TEST_FUNCTION_FAIL("dim = %wu, degree = %wu, p = %wu\n",
                    rdim, order, prime);
    }

    /* structured matrices, small dimensions, small primes */
    for (i = 0; i < 400 * flint_test_multiplier(); i++)
    {
        const ulong prime = _htest_rand_prime(state);
        const slong n = n_randint(state, 13);
        const slong len = 1 + n_randint(state, 7);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, n, n, prime);
        _htest_rand_square(mat, len, state);

        result = core_test_determinant(mat);

        nmod_poly_mat_clear(mat);

        if (!result)
            TEST_FUNCTION_FAIL("(structured) dim = %wd, len = %wd, p = %wu\n",
                    n, len, prime);
    }

    /* unbalanced degrees: one row (or one column) of degree < len, the
     * others of degree < 2 */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 62), 1);
        const slong n = 2 + n_randint(state, 30);
        const slong len = 2 + n_randint(state, 60);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, n, n, prime);
        slong * lens = flint_malloc(n * sizeof(slong));
        for (slong j = 0; j < n; j++)
            lens[j] = 2;
        lens[n_randint(state, n)] = len;
        if (i % 2 == 0)
            nmod_poly_mat_rand_row_degree(mat, state, lens);
        else
            nmod_poly_mat_rand_column_degree(mat, state, lens);
        flint_free(lens);

        result = core_test_determinant(mat);

        nmod_poly_mat_clear(mat);

        if (!result)
            TEST_FUNCTION_FAIL("(unbalanced) dim = %wd, len = %wd, p = %wu\n",
                    n, len, prime);
    }

    /* a few larger matrices, for the recursive algorithms */
    for (i = 0; i < 4 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 20 + n_randint(state, 44), 1);
        const slong n = 24 + n_randint(state, 25);
        const slong len = 1 + n_randint(state, 6);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, n, n, prime);
        _htest_rand_square(mat, len, state);

        result = core_test_determinant(mat);

        nmod_poly_mat_clear(mat);

        if (!result)
            TEST_FUNCTION_FAIL("(larger) dim = %wd, len = %wd, p = %wu\n",
                    n, len, prime);
    }

    TEST_FUNCTION_END(state);
}
