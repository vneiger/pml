#include <stdlib.h>
#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_extra.h"

/** checks X is in reduced row echelon form, defining pivot as rightmost
 * nonzero entry in a row; assumes X has full row rank */
int is_rref(nmod_mat_t X)
{
    slong * pivot = malloc(X->r * sizeof(slong));
    for (slong i = 0; i < X->r; ++i)
    {
        pivot[i] = X->c - 1;
        while (pivot[i] >= 0 && nmod_mat_entry(X, i, pivot[i]) == 0)
            --pivot[i];
        if (pivot[i] < 0)
        {
            printf("found zero row in matrix\n");
            return 0;
        }
        if (i>0 && pivot[i] <= pivot[i-1])
        {
            printf("pivots not increasing\n");
            return 0;
        }
        for (slong j = i+1; j < X->r; ++j)
        {
            if (nmod_mat_entry(X,j,pivot[i]))
            {
                printf("entries not zero below pivot\n");
                return 0;
            }
        }
    }
    return 1;
}

int check(slong field_prime, slong iterations, flint_rand_t state, slong nrows, slong ncols)
{
    nmod_mat_t mat;
    nmod_mat_init(mat, nrows, ncols, field_prime);
    slong rank = FLINT_MIN(3, (FLINT_MIN(nrows,ncols)));
    slong true_rank;
    for (slong k = 0; k < iterations; ++k)
    {
        // random uniform
        flint_printf("Uniformly random matrix:\n");
        nmod_mat_rand(mat, state);
        nmod_mat_print_pretty(mat);

        // random prescribed rank

        flint_printf("\nRandom matrix of rank 2:\n");
        nmod_mat_randrank_dense(mat, state, rank);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        // random full rank

        flint_printf("\nRandom matrix with full rank:\n");
        nmod_mat_rand_fullrank(mat, state);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        // random row echelon

        flint_printf("\nRandom matrix in lower row echelon form:\n");
        nmod_mat_rand_lref(mat, state, rank, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in lower reduced row echelon form:\n");
        nmod_mat_rand_lrref(mat, state, rank);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in lower full rank row echelon form:\n");
        nmod_mat_rand_fullrank_lref(mat, state, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in lower full rank reduced row echelon form:\n");
        nmod_mat_rand_fullrank_lrref(mat, state);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper row echelon form:\n");
        nmod_mat_rand_uref(mat, state, rank, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper reduced row echelon form:\n");
        nmod_mat_rand_urref(mat, state, rank);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper full rank reduced row echelon form:\n");
        nmod_mat_rand_fullrank_uref(mat, state, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper full rank row echelon form:\n");
        nmod_mat_rand_fullrank_urref(mat, state);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        // random column echelon

        flint_printf("\nRandom matrix in lower column echelon form:\n");
        nmod_mat_rand_lcef(mat, state, rank, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in lower reduced column echelon form:\n");
        nmod_mat_rand_lrcef(mat, state, rank);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in lower full rank column echelon form:\n");
        nmod_mat_rand_fullrank_lcef(mat, state, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in lower full rank reduced column echelon form:\n");
        nmod_mat_rand_fullrank_lrcef(mat, state);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper column echelon form:\n");
        nmod_mat_rand_ucef(mat, state, rank, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper reduced column echelon form:\n");
        nmod_mat_rand_urcef(mat, state, rank);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper full rank column echelon form:\n");
        nmod_mat_rand_fullrank_ucef(mat, state, 0);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

        flint_printf("\nRandom matrix in upper full rank reduced column echelon form:\n");
        nmod_mat_rand_fullrank_urcef(mat, state);
        nmod_mat_print_pretty(mat);
        true_rank = nmod_mat_rank(mat);
        flint_printf("\t(rank = %ld)\n", true_rank);

    }
    nmod_mat_clear(mat);
    return 1;
}

int main(int argc, char *argv[])
{
    if (argc != 4 && argc != 5)
    {
        printf("Usage: %s [modulus] [nrows] [ncols] OR %s [modulus] [nrows] [ncols] [nb iterations=1]\n", argv[0], argv[0]);
        return 1;
    }

    srand(time(NULL));
    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, rand(), rand());

    slong field_prime = atol(argv[1]);
    slong iterations = (argc==5) ? atol(argv[4]) : 1;

    if (check(field_prime, iterations, state, atol(argv[2]), atol(argv[3])) == 0)
        printf("BUG\n");

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
