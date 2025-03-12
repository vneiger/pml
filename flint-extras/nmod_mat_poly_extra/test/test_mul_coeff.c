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
    nmod_mat_poly_init(mat2, 3, 2, 7);
    nmod_mat_poly_rand(mat1, state, 5);
    nmod_mat_poly_rand(mat2, state, 2);
    printf("\nFirst matrix:\n");
    nmod_mat_poly_print_pretty(mat1);
    printf("\nSecond matrix:\n");
    nmod_mat_poly_print_pretty(mat2);

    nmod_mat_t coeff;
    nmod_mat_init(coeff, mat1->r, mat2->c, mat1->mod.n);
    printf("\nCoefficient of degree 0 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 0);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 1 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 1);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 2 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 2);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 3 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 3);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 4 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 4);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 5 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 5);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 6 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 6);
    nmod_mat_print_pretty(coeff);

    printf("\nCoefficient of degree 12 of the product:\n");
    nmod_mat_poly_mul_coeff(coeff, mat1, mat2, 12);
    nmod_mat_print_pretty(coeff);

    // convert to and from to check result
    nmod_poly_mat_t pmat1, pmat2, pprod;
    nmod_poly_mat_init(pmat1, mat1->r, mat1->c, mat1->mod.n);
    nmod_poly_mat_init(pmat2, mat2->r, mat2->c, mat2->mod.n);
    nmod_poly_mat_init(pprod, mat1->r, mat2->c, mat2->mod.n);

    nmod_poly_mat_set_from_mat_poly(pmat1, mat1);
    nmod_poly_mat_set_from_mat_poly(pmat2, mat2);
    printf("\nConverted first matrix:\n");
    nmod_poly_mat_print_pretty(pmat1, "X");
    printf("\nConverted second matrix:\n");
    nmod_poly_mat_print_pretty(pmat2, "X");

    nmod_poly_mat_mul(pprod, pmat1, pmat2);
    printf("\nFull product as a polmat:\n");
    nmod_poly_mat_print_pretty(pprod, "X");

    nmod_mat_poly_t prod;
    nmod_mat_poly_init(prod, pprod->r, pprod->c, pprod->modulus);
    printf("OK1\n");
    nmod_mat_poly_set_from_poly_mat(prod, pprod);

    printf("\nFull product as a matpol:\n");
    nmod_mat_poly_print_pretty(prod);

    // clears memory
    nmod_poly_mat_clear(pmat1);
    nmod_poly_mat_clear(pmat2);
    nmod_poly_mat_clear(pprod);
    nmod_mat_poly_clear(prod);

    nmod_mat_clear(coeff);
    nmod_mat_poly_clear(mat1);
    nmod_mat_poly_clear(mat2);
    flint_rand_clear(state);

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
