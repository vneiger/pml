#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"

int check(slong field_prime, slong iterations, flint_rand_t state, slong nrows, slong ncols, slong len)
{
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, nrows, ncols, field_prime);
    slong rank = FLINT_MIN(3, (FLINT_MIN(nrows,ncols)));
    slong true_rank;

    for (slong k = 0; k < iterations; ++k)
    {
        // random uniform
        flint_printf("Uniformly random matrix:\n");
        nmod_poly_mat_rand(mat, state, len);
        nmod_poly_mat_print_pretty(mat, "X");

    }
    nmod_poly_mat_clear(mat);
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
