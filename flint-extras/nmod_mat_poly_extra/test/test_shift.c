#include <time.h>
#include <stdlib.h>

#include <flint/flint.h>

#include "nmod_mat_poly.h"

int main(int argc, char *argv[])
{
    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    printf("building random matrices over GF(7)...\n");
    nmod_mat_poly_t mat1, mat2;
    nmod_mat_poly_init(mat1, 2, 3, 7);
    nmod_mat_poly_init(mat2, 4, 2, 7);
    nmod_mat_poly_rand(mat1, state, 5);
    nmod_mat_poly_rand(mat2, state, 3);
    printf("\nFirst matrix:\n");
    nmod_mat_poly_print_pretty(mat1);

    printf("\nAfter left shift by 5:\n");
    nmod_mat_poly_shift_left(mat1, mat1, 5);
    nmod_mat_poly_print_pretty(mat1);

    printf("\nSecond matrix:\n");
    nmod_mat_poly_print_pretty(mat2);

    printf("\nAfter left shift by 3:\n");
    nmod_mat_poly_shift_left(mat2, mat2, 3);
    nmod_mat_poly_print_pretty(mat2);


    printf("\nZero matrix:\n");
    nmod_mat_poly_zero(mat2);
    nmod_mat_poly_print_pretty(mat2);

    printf("\nAfter left shift by 4:\n");
    nmod_mat_poly_shift_left(mat2, mat2, 4);
    nmod_mat_poly_print_pretty(mat2);


    printf("\nOne matrix:\n");
    nmod_mat_poly_one(mat2);
    nmod_mat_poly_print_pretty(mat2);

    printf("\nAfter left shift by 2:\n");
    nmod_mat_poly_shift_left(mat2, mat2, 2);
    nmod_mat_poly_print_pretty(mat2);


    // clears memory
    nmod_mat_poly_clear(mat1);
    nmod_mat_poly_clear(mat2);
    flint_rand_clear(state);

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
