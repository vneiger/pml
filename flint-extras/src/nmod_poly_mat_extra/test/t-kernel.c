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
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h"

// test one given input
/* TODO does not test generation */
/* TODO does not test reducedness */
/* TODO does not seem tested with shifts */
int core_test_kernel_zls(const nmod_poly_mat_t mat)
{
    slong m = mat->r;
    slong n = mat->c;

    /* careful: this is currently the column degree due to working on the right */
    slong * rdeg = FLINT_ARRAY_ALLOC(n, slong);
    nmod_poly_mat_column_degree(rdeg, mat, NULL);

    nmod_poly_mat_t N;
    nmod_poly_mat_init(N, n, n, mat->modulus);

    slong degN[n];

    slong nz = nmod_poly_mat_kernel_zls(N, degN, mat, NULL, 2.);

    nmod_poly_mat_t Nt;
    nmod_poly_mat_init(Nt, nz, n, mat->modulus);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < nz; j++)
            nmod_poly_set(nmod_poly_mat_entry(Nt, j, i), nmod_poly_mat_entry(N, i, j));        

    nmod_poly_mat_t Mt;
    nmod_poly_mat_init(Mt, n, m, mat->modulus);
    nmod_poly_mat_transpose(Mt, mat);

    int verif = nmod_poly_mat_is_kernel(Nt, Mt, rdeg, ROW_LOWER);

    nmod_poly_mat_clear(N);
    nmod_poly_mat_clear(Nt);
    nmod_poly_mat_clear(Mt);
    flint_free(rdeg);

    return verif;
}

TEST_FUNCTION_START(nmod_poly_mat_kernel_zls, state)
{
    int i,result;

    for (i = 0; i < 16 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 4);
        ulong rdim = 1 + n_randint(state, 4);
        ulong cdim = rdim + 1 + n_randint(state, 2);
        ulong deg = n_randint(state, 5);
        /* ulong nbits = 2 + n_randint(state, 30); */
        /* ulong rdim = 1 + n_randint(state, 60); */
        /* ulong cdim = rdim + 1 + n_randint(state, 20); */
        /* ulong deg = n_randint(state, 20); */

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

        result = core_test_kernel_zls(A);

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

        result = core_test_kernel_zls(A);

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

        result = core_test_kernel_zls(A);

        nmod_poly_mat_clear(A);

        if (!result)
        {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

    TEST_FUNCTION_END(state);
}
