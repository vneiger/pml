/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
//#include <flint/profiler.h>
#include <time.h>
#include <flint/ulong_extras.h>
#include <flint/test_helpers.h>

#include <nmod_poly_mat_utils.h>
#include <nmod_poly_mat_io.h>

//#include "../../nmod_poly_mat_extras/test/testing_collection.h"

#include "nmod_poly_mat_dixon.h"

// Random poly_mat which is nonsingular for x=0
void nmod_poly_mat_rand_origin(nmod_poly_mat_t A, flint_rand_t state, slong order) 
{

    nmod_mat_t T1,T2;
    nmod_mat_init(T1, A->r, A->r, A->modulus);
    nmod_mat_init(T2, A->r, A->r, A->modulus);

    nmod_mat_randtril(T1,state,0);
    nmod_mat_randtriu(T2,state,0);

    nmod_mat_mul(T1,T1,T2);

    nmod_poly_mat_t B;
    nmod_poly_mat_init(B, A->r, A->r, A->modulus);

    // nmod_poly_mat_set_nmod_mat(B,T1); // GV does not work, bug ? Sam 31 mai 2025 14:29:48 CEST

    for (slong i = 0; i < A->r; i++)
    {
        for (slong j = 0; j < A->r; j++)
        {
            if (nmod_mat_entry(T1, i, j) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(B, i, j));
            else
            {
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, i, j),0,nmod_mat_entry(T1, i, j));
            }
        }
    }

    nmod_poly_mat_rand(A, state, order-1);
    nmod_poly_mat_shift_left(A,A,1);

    nmod_poly_mat_add(A,A,B);

    nmod_mat_clear(T1);
    nmod_mat_clear(T2); 
    nmod_poly_mat_clear(B); 
}


// test one given input
int core_test_dixon(const nmod_poly_mat_t A, const nmod_poly_mat_t B, slong order, slong sigma)
{

    nmod_poly_mat_t X;
    nmod_poly_mat_init(X, B->r, B->c, A->modulus);

    nmod_poly_mat_dixon(X, A, B, order, sigma);

    nmod_poly_mat_mul(X, A, X);
    nmod_poly_mat_sub(X, X, B);
    nmod_poly_mat_truncate(X,sigma);


    // printf("\n");
    // nmod_poly_mat_print_pretty(X, "x");
    // printf("\n");

    if (nmod_poly_mat_is_zero(X) !=0) 
    {
        nmod_poly_mat_clear(X);  
        return 1;
    }
    else 
    {
        nmod_poly_mat_clear(X);  
        return 0;
    }
}


//Test against several parameters 
int collection_test_dixon(flint_rand_t state)
{

    int res;

// #if FLINT64
//     slong prime[]  = 
//     {
//     2, 3, 5,    // very small ones
//     521, 524309, 536870923,  // 10, 20, 30 bits
//     549755813911, 562949953421381, 576460752303423619  // 40, 50, 60 bits
//     };
// #else 
// slong prime[]  = 
//     {
//     2, 3, 5,    // very small ones
//     521, 524309, 536870923  // 10, 20, 30 bits
//     };
// #endif 

//     slong rdim[] = {2, 5, 10, 20, 30};

//     slong Bcdim[] = {1, 2, 10};

//     slong degA[] = {1, 2, 10};

//     slong degB[] = {1, 2, 100};

#if FLINT64
    slong prime[]  = 
    {
    3,    // very small ones
    521,  // 10, 20, 30 bits
    576460752303423619  // 40, 50, 60 bits
    };
#else 
slong prime[]  = 
    {
    2, 5,    // very small ones
    521, 524309 // 10, 20, 30 bits
    };
#endif 

    slong rdim[] = {2, 10, 20};

    slong Bcdim[] = {1, 10};

    slong degA[] = {1, 12};

    slong degB[] = {1, 8};

    slong order;
    slong sigma;

    //printf("Launching testing collection\n");

    for (slong ip=0; ip<sizeof(prime)/sizeof(prime[0]); ip++)
        for (slong in=0; in<sizeof(rdim)/sizeof(rdim[0]); in++)
            for (slong ida=0; ida<sizeof(degA)/sizeof(degA[0]); ida++)
                for (slong bdim=0; bdim<sizeof(Bcdim)/sizeof(Bcdim[0]); bdim++)
                    for (slong db=0; db<sizeof(degB)/sizeof(degB[0]); db++)
                    {

                        sigma = 2*rdim[in]*degA[ida];
                        order = degA[ida];

                        /* //printf("prime %ld, rdim %ld, degA %ld, order %ld, sigma %ld, Bcdim %ld, degB %ld. \n",\ */
                        /*     //prime[ip], rdim[in], degA[ida], order,\ */
                        /*     //sigma, Bcdim[bdim], degB[db]); */

                        nmod_poly_mat_t A;
                        nmod_poly_mat_t B;

                        nmod_poly_mat_init(A, rdim[in], rdim[in], prime[ip]);
                        nmod_poly_mat_init(B, rdim[in], Bcdim[bdim], prime[ip]);

                        nmod_poly_mat_rand_origin(A, state, degA[ida]+1);
                        nmod_poly_mat_rand(B, state, degB[db]+1);

                        res = core_test_dixon(A, B, order, sigma);

                        if (res == 0 ) 
                        {
                            printf("failed ... exiting\n"); 
                            return 0;
                        }

                        nmod_poly_mat_clear(A);
                        nmod_poly_mat_clear(B);
                    }

    for (slong ip=0; ip<sizeof(prime)/sizeof(prime[0]); ip++)
        for (slong in=0; in<sizeof(rdim)/sizeof(rdim[0]); in++)
            for (slong ida=0; ida<sizeof(degA)/sizeof(degA[0]); ida++)
                for (slong bdim=0; bdim<sizeof(Bcdim)/sizeof(Bcdim[0]); bdim++)
                    for (slong db=0; db<sizeof(degB)/sizeof(degB[0]); db++)
                    {

                        sigma = 2*rdim[in]*degA[ida];
                        order = 2;

                        /* //printf("prime %ld, rdim %ld, degA %ld, order %ld, sigma %ld, Bcdim %ld, degB %ld. \n",\ */
                        /*     //prime[ip], rdim[in], degA[ida], order,\ */
                        /*     //sigma, Bcdim[bdim], degB[db]); */

                        nmod_poly_mat_t A;
                        nmod_poly_mat_t B;

                        nmod_poly_mat_init(A, rdim[in], rdim[in], prime[ip]);
                        nmod_poly_mat_init(B, rdim[in], Bcdim[bdim], prime[ip]);

                        nmod_poly_mat_rand_origin(A, state, degA[ida]+1);
                        nmod_poly_mat_rand(B, state, degB[db]+1);

                        res = core_test_dixon(A, B, order, sigma);

                        if (res == 0 ) 
                        {
                            printf("failed ... exiting\n"); 
                            return 0;
                        }

                        nmod_poly_mat_clear(A);
                        nmod_poly_mat_clear(B);
                    }
return 1;
    
}


TEST_FUNCTION_START(nmod_poly_mat_dixon, state)
{

    int res=0;  

    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    res=collection_test_dixon(state);

    if (res == 0)
    {
       TEST_FUNCTION_FAIL("");
    }
    else
    {
    TEST_FUNCTION_END(state);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
