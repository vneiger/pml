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

    printf("building random matrix over GF(7)...\n");
    nmod_mat_poly_t matp;
    nmod_mat_poly_init(matp, 2, 3, 7);
    nmod_mat_poly_rand(matp, state, 3);

    printf("Basic print:\n");
    nmod_mat_poly_print(matp);
    printf("\nPretty print:\n");
    nmod_mat_poly_print_pretty(matp);

    printf("\nclearing, building new random matrix over GF(2)...");
    nmod_mat_poly_clear(matp);
    nmod_mat_poly_init(matp, 1, 2, 2);
    nmod_mat_poly_randtest(matp, state, 3);
	
    printf("\nPretty print:\n");
    nmod_mat_poly_print_pretty(matp);

    printf("\nRealloc (more):\n");
    printf(" before --> length: %ld, alloc: %ld\n", matp->length, matp->alloc);
    nmod_mat_poly_realloc(matp, 10);
    printf(" after  --> length: %ld, alloc: %ld\n", matp->length, matp->alloc);
    printf(" after, pretty print:\n");
    nmod_mat_poly_print_pretty(matp);

    printf("\nRealloc (less):\n");
    printf(" before --> length: %ld, alloc: %ld\n", matp->length, matp->alloc);
    nmod_mat_poly_realloc(matp, 2);
    printf(" after  --> length: %ld, alloc: %ld\n", matp->length, matp->alloc);
    printf(" after, pretty print:\n");
    nmod_mat_poly_print_pretty(matp);

    // clears memory
    nmod_mat_poly_clear(matp);
    flint_rand_clear(state);

	return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
