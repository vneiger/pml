#include <time.h>
#include <stdlib.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

int check(slong field_prime, slong iterations, flint_rand_t state, slong nrows, slong ncols, slong len)
{
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, nrows, ncols, field_prime);

    srand(time(NULL));

    slong * rdeg = flint_malloc(nrows * sizeof(slong));
    slong * cdeg = flint_malloc(ncols * sizeof(slong));
    for (slong i = 0; i < nrows; i++)
        rdeg[i] = rand() % 100;
    for (slong i = 0; i < ncols; i++)
        cdeg[i] = rand() % 100;

    slong * urshift = flint_calloc(ncols, sizeof(slong));
    slong * ucshift = flint_calloc(nrows, sizeof(slong));
    slong * rshift = flint_malloc(ncols * sizeof(slong));
    slong * cshift = flint_malloc(nrows * sizeof(slong));
    for (slong i = 0; i < ncols; i++)
        rshift[i] = (rand() % 100) - 50;
    for (slong i = 0; i < nrows; i++)
        cshift[i] = (rand() % 100) - 50;

    slong * rpivind = flint_malloc(nrows * sizeof(slong));
    slong * rpivdeg = flint_malloc(nrows * sizeof(slong));
    slong * cpivind = flint_malloc(ncols * sizeof(slong));
    slong * cpivdeg = flint_malloc(ncols * sizeof(slong));
    for (slong i = 0; i < nrows/2; i++)
    {
        rpivind[i] = i;
        rpivdeg[i] = rand() % 100;
    }
    for (slong i = nrows/2; i < nrows; i++)
    {
        rpivind[i] = ncols-nrows+i;
        rpivdeg[i] = rand() % 100;
    }
    for (slong i = 0; i < ncols/2; i++)
    {
        cpivind[i] = i;
        cpivdeg[i] = rand() % 100;
    }
    for (slong i = ncols/2; i < ncols; i++)
    {
        cpivind[i] = nrows-ncols+i;
        cpivdeg[i] = rand() % 100;
    }



    fmpz_mat_t dmat;
    fmpz_mat_init(dmat, nrows, ncols);
    for (slong i = 0; i < nrows; i++)
        for (slong j = 0; j < ncols; j++)
            fmpz_randtest_unsigned(fmpz_mat_entry(dmat,i,j), state, 5);

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

    for (slong k = 0; k < iterations; ++k)
    {
        flint_printf("\n(Showing degree matrices)\n");

        // random uniform
        flint_printf("\nUniformly random matrix, with prescribed degree:\n");
        nmod_poly_mat_rand(mat, state, len);
        //nmod_poly_mat_print_pretty(mat, "x");
        nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random uniform, rdeg
        flint_printf("\nUniformly random matrix, with prescribed row degree:\n");
        nmod_poly_mat_rand_row_degree(mat, state, rdeg);
        //nmod_poly_mat_print_pretty(mat, "x");
        nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random uniform, cdeg
        flint_printf("\nUniformly random matrix, with prescribed column degree:\n");
        nmod_poly_mat_rand_column_degree(mat, state, cdeg);
        //nmod_poly_mat_print_pretty(mat, "x");
        nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random uniform, dmat
        flint_printf("\nUniformly random matrix, with prescribed degree matrix:\n");
        nmod_poly_mat_rand_degree_matrix(mat, state, dmat);
        //nmod_poly_mat_print_pretty(mat, "x");
        nmod_poly_mat_degree_matrix_print_pretty(mat);

        // random Popov, row-wise
        if (nrows <= ncols)
        {
            flint_printf("\nUniformly random matrix, row-wise 0-Popov:\n");
            _nmod_poly_mat_rand_popov_row_lower(mat, state, rpivind, rpivdeg, urshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, urshift, ROW_LOWER))
                return 1;

            flint_printf("\nUniformly random matrix, row-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, rpivind, rpivdeg, NULL, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER))
                return 1;

            flint_printf("\nUniformly random matrix, row-wise shifted-Popov:\n");
            _nmod_poly_mat_rand_popov_row_lower(mat, state, rpivind, rpivdeg, rshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, rshift, ROW_LOWER))
                return 1;

            flint_printf("\nUniformly random matrix, row-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, rpivind, rpivdeg, rshift, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, rshift, ROW_LOWER))
                return 1;
        }

        // random Popov, column-wise
        if (ncols <= nrows)
        {
            flint_printf("\nUniformly random matrix, column-wise 0-Popov:\n");
            _nmod_poly_mat_rand_popov_col_upper(mat, state, cpivind, cpivdeg, ucshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, ucshift, COL_UPPER))
                return 1;

            flint_printf("\nUniformly random matrix, column-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, cpivind, cpivdeg, NULL, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, COL_UPPER))
                return 1;

            flint_printf("\nUniformly random matrix, column-wise shifted-Popov:\n");
            _nmod_poly_mat_rand_popov_col_upper(mat, state, cpivind, cpivdeg, cshift);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, cshift, COL_UPPER))
                return 1;

            flint_printf("\nUniformly random matrix, column-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, cpivind, cpivdeg, cshift, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, cshift, COL_UPPER))
                return 1;
        }

        // general interface with square case, to test NULL arguments
        if (ncols == nrows)
        {
            flint_printf("\nUniformly random square matrix, row-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, rpivdeg, NULL, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER))
                return 1;

            flint_printf("\nUniformly random square matrix, row-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, rpivdeg, rshift, ROW_LOWER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, rshift, ROW_LOWER))
                return 1;

            flint_printf("\nUniformly random square matrix, column-wise 0-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, cpivdeg, NULL, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, NULL, COL_UPPER))
                return 1;

            flint_printf("\nUniformly random square matrix, column-wise shifted-Popov, with general interface:\n");
            nmod_poly_mat_rand_popov(mat, state, NULL, cpivdeg, cshift, COL_UPPER);
            //nmod_poly_mat_print_pretty(mat, "x");
            nmod_poly_mat_degree_matrix_print_pretty(mat);

            if (!nmod_poly_mat_is_popov(mat, cshift, COL_UPPER))
                return 1;
        }
    }
    nmod_poly_mat_clear(mat);
    flint_free(rdeg);
    flint_free(cdeg);
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 5 && argc != 6)
    {
        printf("Usage: %s [modulus] [nrows] [ncols] [len] OR %s [modulus] [nrows] [ncols] [len] [nb iterations=1]\n", argv[0], argv[0]);
        return 1;
    }

    srand(time(NULL));
    flint_rand_t state;
    flint_rand_init(state);
    flint_randseed(state, rand(), rand());

    slong field_prime = atol(argv[1]);
    slong iterations = (argc==6) ? atol(argv[5]) : 1;

    if (check(field_prime, iterations, state, atol(argv[2]), atol(argv[3]), atol(argv[4])) != 0)
        printf("BUG\n");

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
