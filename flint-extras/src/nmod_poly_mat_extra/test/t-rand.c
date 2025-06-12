/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

int check_rand_form(slong field_prime, slong iterations, flint_rand_t state, slong nrows, slong ncols, slong len)
{
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, nrows, ncols, field_prime);

    slong * rdeg = flint_malloc(nrows * sizeof(slong));
    slong * cdeg = flint_malloc(ncols * sizeof(slong));
    for (slong i = 0; i < nrows; i++)
        rdeg[i] = n_randint(state, 100);
    for (slong i = 0; i < ncols; i++)
        cdeg[i] = n_randint(state, 100);

    slong * urshift = flint_calloc(ncols, sizeof(slong));
    slong * ucshift = flint_calloc(nrows, sizeof(slong));
    slong * rshift = flint_malloc(ncols * sizeof(slong));
    slong * cshift = flint_malloc(nrows * sizeof(slong));
    for (slong i = 0; i < ncols; i++)
        rshift[i] = (slong)n_randint(state, 100) - 50;
    for (slong i = 0; i < nrows; i++)
        cshift[i] = (slong)n_randint(state, 100) - 50;

    slong * rpivind = flint_malloc(nrows * sizeof(slong));
    slong * rpivdeg = flint_malloc(nrows * sizeof(slong));
    slong * cpivind = flint_malloc(ncols * sizeof(slong));
    slong * cpivdeg = flint_malloc(ncols * sizeof(slong));
    for (slong i = 0; i < nrows/2; i++)
    {
        rpivind[i] = i;
        rpivdeg[i] = n_randint(state, 100);
    }
    for (slong i = nrows/2; i < nrows; i++)
    {
        rpivind[i] = ncols-nrows+i;
        rpivdeg[i] = n_randint(state, 100);
    }
    for (slong i = 0; i < ncols/2; i++)
    {
        cpivind[i] = i;
        cpivdeg[i] = n_randint(state, 100);
    }
    for (slong i = ncols/2; i < ncols; i++)
    {
        cpivind[i] = nrows-ncols+i;
        cpivdeg[i] = n_randint(state, 100);
    }

    fmpz_mat_t dmat;
    fmpz_mat_init(dmat, nrows, ncols);
    for (slong i = 0; i < nrows; i++)
        for (slong j = 0; j < ncols; j++)
            fmpz_randtest_unsigned(fmpz_mat_entry(dmat,i,j), state, 5);

    int verbose = 0;
    if (verbose)
    {
        flint_printf("nrows\tncols\tlen\n");
        flint_printf("%ld\t%ld\t%ld\n", nrows, ncols, len);
        flint_printf("rdeg:\n");
        slongvec_print_sagemath(rdeg, nrows);
        flint_printf("cdeg:\n");
        slongvec_print_sagemath(cdeg, ncols);
        flint_printf("dmat:\n");
        fmpz_mat_print_pretty(dmat);
        flint_printf("\nrshift:\n");
        slongvec_print_sagemath(rshift, ncols);
        flint_printf("cshift:\n");
        slongvec_print_sagemath(cshift, nrows);
        flint_printf("rpivind\n");
        slongvec_print_sagemath(rpivind, nrows);
        flint_printf("rpivdeg\n");
        slongvec_print_sagemath(rpivdeg, nrows);
        flint_printf("cpivind\n");
        slongvec_print_sagemath(cpivind, ncols);
        flint_printf("cpivdeg\n");
        slongvec_print_sagemath(cpivdeg, ncols);
        printf("\n");
    }

    for (slong k = 0; k < iterations; ++k)
    {
        if (verbose) flint_printf("\n(Showing degree matrices)\n");

        // random uniform
        if (verbose) flint_printf("\nUniformly random matrix, with prescribed degree:\n");
        nmod_poly_mat_rand(mat, state, len);
        //nmod_poly_mat_print_pretty(mat, "x");
        if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random uniform, rdeg
        if (verbose) flint_printf("\nUniformly random matrix, with prescribed row degree:\n");
        nmod_poly_mat_rand_row_degree(mat, state, rdeg);
        //nmod_poly_mat_print_pretty(mat, "x");
        if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random uniform, cdeg
        if (verbose) flint_printf("\nUniformly random matrix, with prescribed column degree:\n");
        nmod_poly_mat_rand_column_degree(mat, state, cdeg);
        //nmod_poly_mat_print_pretty(mat, "x");
        if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random uniform, dmat
        if (verbose) flint_printf("\nUniformly random matrix, with prescribed degree matrix:\n");
        nmod_poly_mat_rand_degree_matrix(mat, state, dmat);
        //nmod_poly_mat_print_pretty(mat, "x");
        if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random Popov, row-wise
        if (nrows <= ncols)
        {
            if (verbose) flint_printf("\nUniformly random matrix, row-wise 0-Popov:\n");
            _nmod_poly_mat_rand_popov_row_lower(mat, state, rpivind, rpivdeg, urshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, urshift, ROW_LOWER))
                return 1;

            if (verbose) flint_printf("\nUniformly random matrix, row-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, rpivind, rpivdeg, NULL, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER))
                return 1;

            if (verbose) flint_printf("\nUniformly random matrix, row-wise shifted-Popov:\n");
            _nmod_poly_mat_rand_popov_row_lower(mat, state, rpivind, rpivdeg, rshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, rshift, ROW_LOWER))
                return 1;

            if (verbose) flint_printf("\nUniformly random matrix, row-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, rpivind, rpivdeg, rshift, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, rshift, ROW_LOWER))
                return 1;
        }

        // random Popov, column-wise
        if (ncols <= nrows)
        {
            if (verbose) flint_printf("\nUniformly random matrix, column-wise 0-Popov:\n");
            _nmod_poly_mat_rand_popov_col_upper(mat, state, cpivind, cpivdeg, ucshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, ucshift, COL_UPPER))
                return 1;

            if (verbose) flint_printf("\nUniformly random matrix, column-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, cpivind, cpivdeg, NULL, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, COL_UPPER))
                return 1;

            if (verbose) flint_printf("\nUniformly random matrix, column-wise shifted-Popov:\n");
            _nmod_poly_mat_rand_popov_col_upper(mat, state, cpivind, cpivdeg, cshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, cshift, COL_UPPER))
                return 1;

            if (verbose) flint_printf("\nUniformly random matrix, column-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, cpivind, cpivdeg, cshift, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, cshift, COL_UPPER))
                return 1;
        }

        // general interface with square case, to test NULL arguments
        if (ncols == nrows)
        {
            if (verbose) flint_printf("\nUniformly random square matrix, row-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, rpivdeg, NULL, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER))
                return 1;

            if (verbose) flint_printf("\nUniformly random square matrix, row-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, rpivdeg, rshift, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, rshift, ROW_LOWER))
                return 1;

            if (verbose) flint_printf("\nUniformly random square matrix, column-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, cpivdeg, NULL, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, COL_UPPER))
                return 1;

            if (verbose) flint_printf("\nUniformly random square matrix, column-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, cpivdeg, cshift, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            if (verbose) nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, cshift, COL_UPPER))
                return 1;
        }
    }
    nmod_poly_mat_clear(mat);
    flint_free(rdeg);
    flint_free(cdeg);
    return 0;
}

TEST_FUNCTION_START(nmod_poly_mat_rand, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 63);
        ulong nrows = 1 + n_randint(state, 30);
        ulong ncols = 1 + n_randint(state, 30);
        ulong len = 1 + n_randint(state, 100);
        ulong prime = n_randprime(state, bits, 1);

        result = check_rand_form(prime, 10, state, nrows, ncols, len);

        if (result)
            TEST_FUNCTION_FAIL("m = %wu, n = %wu, prime = %wu\n", nrows, ncols, prime);
    }

    TEST_FUNCTION_END(state);
}
