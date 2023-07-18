#include <flint/flint.h>
#include "nmod_mat_poly.h"

int main(int argc, char *argv[])
{
    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    const long prime = 7;
    const long rdim = 4;
    const long cdim = 2;
    const long order = 4;

    printf("building random matrix over GF(7)...\n");
    nmod_mat_poly_t mat, appbas;
    nmod_mat_poly_init(appbas, rdim, rdim, prime);
    nmod_mat_poly_init(mat, rdim, cdim, prime);
    nmod_mat_poly_rand(mat, state, order);
    printf("\nInput matrix:\n");
    nmod_mat_poly_print_pretty(mat);

    slong shift[] = {0,0,0,0};

    nmod_mat_poly_mbasis(appbas, shift, mat, 4);

    nmod_mat_t coeff;
    nmod_mat_init(coeff, mat->r, mat->c, mat->mod.n);
    printf("\nCoefficient of degree 0 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, appbas, mat, 0);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 1 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, appbas, mat, 1);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 2 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, appbas, mat, 2);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 3 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, appbas, mat, 3);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 4 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, appbas, mat, 4);
    nmod_mat_print_pretty(coeff);

    printf("\nAppbas:\n");
    nmod_mat_poly_print_pretty(appbas);

    // clears memory
    nmod_mat_poly_clear(mat);
    nmod_mat_poly_clear(appbas);
    nmod_mat_clear(coeff);
    flint_randclear(state);

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
