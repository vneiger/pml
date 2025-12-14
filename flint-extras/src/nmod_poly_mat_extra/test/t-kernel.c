/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_extra.h"  

// test one given input
int core_test_kernel(const nmod_poly_mat_t mat)
{

    slong m=mat->r;
    slong n=mat->c;
    
    // Flint rank 
    // ----------
    
    slong rkflint;

    rkflint=nmod_poly_mat_rank(mat);
    

    // PML nullspace
    // -------------

    nmod_poly_mat_t N; 
    nmod_poly_mat_init(N, n, n, mat->modulus);

    slong nz;

    int i,j;

    // slong iz[m];

    // for (i = 0; i < m; i++) 
    // {
    //     iz[i]=0; 
    // }

    slong degN[n];

    nz=nmod_poly_mat_kernel(N, degN, mat, NULL, 2);

    int verif;
    verif =  (n-rkflint == nz);

    if (nz !=0) {

        nmod_poly_mat_t NN; 
        nmod_poly_mat_init(NN, n, nz, mat->modulus);

        for (i = 0; i < n; i++)
            for (j = 0; j < nz; j++) {
                nmod_poly_set(nmod_poly_mat_entry(NN, i, j), nmod_poly_mat_entry(N, i, j));
            }

        nmod_poly_mat_t Z;
        nmod_poly_mat_init(Z, m, nz, mat->modulus);

        nmod_poly_mat_mul(Z, mat, NN);

        verif = verif && nmod_poly_mat_is_zero(Z); 

        nmod_poly_mat_clear(Z);
        nmod_poly_mat_clear(NN);

    }
    
    return verif; 
}


TEST_FUNCTION_START(nmod_poly_mat_kernel, state)
{
    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    int i,result;

    for (i = 0; i < 16 * flint_test_multiplier(); i++)
    {

        ulong nbits = 2 + n_randint(state, 30);
        ulong rdim = 1 + n_randint(state, 60);
        ulong cdim = rdim + 1 + n_randint(state, 20);
        ulong deg = n_randint(state, 20);

        ulong prime = n_randprime(state, nbits, 1);


        nmod_poly_mat_t A;

        if (i < 4) {
            cdim=rdim;
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 1.0);
        }
        else if (i < 8) {
            cdim=rdim+1;
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.8);
        }
        else if (i < 12) { 
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.2);
        }
        else {
            nmod_poly_mat_init(A, rdim, cdim, prime);
            nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.84);
        }     

        result = core_test_kernel(A);

        nmod_poly_mat_clear(A);

        if (!result) {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

    for (i = 0; i < 2; i++)
    {
        ulong rdim = 40 + n_randint(state, 8);
        ulong cdim = rdim + n_randint(state, 8);
        ulong deg = 0;

        ulong prime = 7; 

        nmod_poly_mat_t A;

        nmod_poly_mat_init(A, rdim, cdim, prime);
        nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.8);

        result = core_test_kernel(A);

        nmod_poly_mat_clear(A);

        if (!result) {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

    for (i = 0; i < 2; i++)
    {
        ulong rdim = 20 + n_randint(state, 8);
        ulong cdim = rdim -8;
        ulong deg = 1;

        ulong prime = 2; 

        nmod_poly_mat_t A;

        nmod_poly_mat_init(A, rdim, cdim, prime);
        nmod_poly_mat_randtest_sparse(A, state, deg+1, 0.2);

        result = core_test_kernel(A);

        nmod_poly_mat_clear(A);

        if (!result) {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, degree = %wu, p = %wu\n", \
                rdim, cdim, deg, prime);
        }
    }

TEST_FUNCTION_END(state);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
