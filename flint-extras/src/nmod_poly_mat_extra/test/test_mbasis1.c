/*
    Copyright (C) 2025 Vincent Neiger, Kevin Tran

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>

#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_io.h" // for print_pretty

#define PRIME_30_BITS 536870923
#define PRIME_60_BITS 576460752303423619

int test_mbasis1(void)
{
    nmod_mat_t mat;
    slong rdim = 9, cdim = 3, prime = 7;

    nmod_mat_init(mat, rdim, cdim, prime);


    flint_rand_t state;
    flint_rand_init(state);

    nmod_mat_randtest(mat, state);


    nmod_poly_mat_t res;
    nmod_poly_mat_init(res, rdim, rdim, prime);
    nmod_poly_mat_zero(res);

    slong res_shift[rdim];

    slong shift[rdim];

    for (slong i = 0; i < rdim; i++)
        shift[i] =  n_randint(state, 10) - 5;

    mbasis1(res, res_shift, mat, shift);


    printf("Matrix\n");
    nmod_mat_print_pretty(mat);
    printf("Result mbasis1\n");
    nmod_poly_mat_print_pretty(res,"x");
    nmod_poly_mat_clear(res);
    flint_rand_clear(state);

    nmod_mat_t res2;
    slong res_shift2[rdim];
    slong res_perm2[rdim];
    slong rank = mbasis1_for_mbasis(res2, res_shift2, res_perm2, mat, shift);

    nmod_mat_clear(mat);

    printf("\nMatrix A from basis = [[x,0],[A,1]]\n");
    nmod_mat_print_pretty(res2);

    printf("\nshifts = ");
    for(slong i = 0; i < rdim; i++)
        printf("%ld ", res_shift2[i]);

    printf("\npermutation = ");
    for(slong i = 0; i < rdim; i++)
        printf("%ld ", res_perm2[i]);


    printf("\nrank = %ld\n", rank);

    nmod_mat_clear(res2);
    return 0;
}

int main(void)
{
    test_mbasis1();
    return 0;
}
