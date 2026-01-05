/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"

// test one given input
int core_test_kernel_zls(const nmod_poly_mat_t mat, flint_rand_t state)
{
    slong n = mat->c;

    /* pick shift with entries >= cdeg */
    slong * cdeg = FLINT_ARRAY_ALLOC(n, slong);
    nmod_poly_mat_column_degree(cdeg, mat, NULL);
    slong * shift = FLINT_ARRAY_ALLOC(n, slong);
    for (slong j = 0; j < n; j++)
        shift[j] = FLINT_MAX(0, cdeg[j]) + n_randint(state, 30);

    nmod_poly_mat_t N;
    nmod_poly_mat_init(N, n, n, mat->modulus);
    slong degN[n];
    slong nz = nmod_poly_mat_kernel_zls(N, degN, mat, shift, 2.);

    nmod_poly_mat_t Nnz;
    nmod_poly_mat_window_init(Nnz, N, 0, 0, n, nz);
    int verif = nmod_poly_mat_is_kernel(Nnz, shift, mat, REDUCED, COL_UPPER);

    nmod_poly_mat_clear(N);
    nmod_poly_mat_window_clear(Nnz);
    flint_free(cdeg);
    flint_free(shift);

    return verif;
}

TEST_FUNCTION_START(nmod_poly_mat_kernel_zls, state)
{
    int i, result;

    for (i = 0; i < 16 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 63);
        ulong rdim = 1 + n_randint(state, 20);
        ulong cdim = rdim + 1 + n_randint(state, 20);
        ulong deg = n_randint(state, 100);

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t A;

        if (i < 4)
        {
            cdim=rdim;
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 1.0);
        }
        else if (i < 8)
        {
            cdim=rdim+1;
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.8);
        }
        else if (i < 12)
        {
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.2);
        }
        else
        {
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.84);
        }

        result = core_test_kernel_zls(A, state);

        nmod_poly_mat_clear(A);

        if (!result)
        {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

    // Degree zero
    // -----------

    for (i = 0; i < 2; i++)
    {
        ulong rdim = 40 + n_randint(state, 8);
        ulong cdim = rdim + n_randint(state, 8);
        ulong deg = 0;

        ulong prime = 7;

        nmod_poly_mat_t A;

        nmod_poly_mat_init(A, rdim, cdim, prime);
        nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.8);

        result = core_test_kernel_zls(A, state);

        nmod_poly_mat_clear(A);

        if (!result)
        {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

    // Column dimension smaller that row dimension
    // -------------------------------------------

    for (i = 0; i < 2; i++)
    {
        ulong rdim = 20 + n_randint(state, 8);
        ulong cdim = rdim -8;
        ulong deg = 1;

        ulong prime = 2;

        nmod_poly_mat_t A;

        nmod_poly_mat_init(A, rdim, cdim, prime);
        nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.2);

        result = core_test_kernel_zls(A, state);

        nmod_poly_mat_clear(A);

        if (!result)
        {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

    TEST_FUNCTION_END(state);
}
