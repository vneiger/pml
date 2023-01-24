#include <flint/fmpz_mat.h>
#include <flint/perm.h>
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
    }
    nmod_poly_mat_clear(mat);
    flint_free(rdeg);
    flint_free(cdeg);
    return 1;
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
    flint_randinit(state);
    flint_randseed(state, rand(), rand());

    slong field_prime = atol(argv[1]);
    slong iterations = (argc==6) ? atol(argv[5]) : 1;

    if (check(field_prime, iterations, state, atol(argv[2]), atol(argv[3]), atol(argv[4])) == 0)
        printf("BUG\n");

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
