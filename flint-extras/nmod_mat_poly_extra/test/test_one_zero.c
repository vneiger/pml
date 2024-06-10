#include "nmod_mat_poly.h"
#include <flint/flint.h>

int main(int argc, char *argv[])
{
    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    printf("building random matrix over GF(7)...\n");
    nmod_mat_poly_t matp;
    nmod_mat_poly_init(matp, 4, 5, 7);
    nmod_mat_poly_rand(matp, state, 3);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if zero: %d\n", nmod_mat_poly_is_zero(matp));

    printf("\nnow making it zero\n");
    nmod_mat_poly_zero(matp);
    nmod_mat_poly_print_pretty(matp);

    printf("testing if zero: %d\n", nmod_mat_poly_is_zero(matp));
    printf("testing if one: %d\n", nmod_mat_poly_is_one(matp));

    printf("\nnow making it one\n");
    nmod_mat_poly_one(matp);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if one: %d\n", nmod_mat_poly_is_one(matp));

    printf("\nmaking it random of larger degree again\n");
    nmod_mat_poly_rand(matp, state, 3);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if one: %d\n", nmod_mat_poly_is_zero(matp));

    printf("\nnow making it one\n");
    nmod_mat_poly_one(matp);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if one: %d\n", nmod_mat_poly_is_one(matp));

    // clears memory
    nmod_mat_poly_clear(matp);
    flint_rand_clear(state);

	return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
