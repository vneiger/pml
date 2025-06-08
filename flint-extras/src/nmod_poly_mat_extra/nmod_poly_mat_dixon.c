#include <stdlib.h>
#include <math.h> 

#include <flint/nmod_poly.h> 

#include "nmod_poly_mat_dixon.h"


/**
 *  TODO: probably only degree > 0 currently 
 * 
 *  Truncated inverse of A mod x^order
 *  A is assumed to be invertible for x=0 (will produce an error otherwise) 
 * 
 *  CHECK TODO: uses a method based on a shifted approximant computation 
 *              check whether the specification of nmod_poly_mat_pmbasis 
 *              is as expected 
 * 
 */

void nmod_poly_mat_inv_trunc(nmod_poly_mat_t S, 
                            const nmod_poly_mat_t A, 
                            ulong order)

{

    slong i,j;

    slong n=A->r;


    // H = [A -Id]

    nmod_poly_mat_t H;
    nmod_poly_mat_init(H, 2*n, n, A->modulus);

    nmod_poly_mat_set(H, A);  // GV, correct even the dimension of H is larger ? 

    nmod_poly_t mone;
    nmod_poly_init(mone, A->modulus);
    nmod_poly_set_coeff_ui(mone, 0, A->modulus -1); // 1


    for (i = n; i < 2*n; i++)
        nmod_poly_set(nmod_poly_mat_entry(H, i, i-n), mone);
   
    // printf("\n");
    // nmod_poly_mat_print_pretty(H, "x");
    // printf("\n");

    // Will be the approximant basis 

    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, 2*n, 2*n, A->modulus);

    // Appropriate shift
    slong shift[2*n];

    for (i = 0; i < n; i++) 
    {
        shift[i]=0; 
        shift[i+n]=order; 
    }

    // Shifted approximant computation 

    nmod_poly_mat_pmbasis(tsf, shift, H, order); 

    // Check whether the lower right part of the approximant is the odentity matrix 

    for (i = 0; i < n; i++)
        if (nmod_poly_is_one(nmod_poly_mat_entry(tsf, n+i, n+i)) != 1) 
        {
            printf("error in nmod_poly_mat_inv_trunc: check of identity, singular at origin?\n");
            exit(-1);
        }

    for (i = 1; i < n; i++)
        for (j = 0; j < i; j++)
            if (nmod_poly_is_zero(nmod_poly_mat_entry(tsf, n+i, n+j)) != 1) 
            {
                printf("error in nmod_poly_mat_inv_trunc: check of identity, singular at origin?\n");
                exit(-1);
            }

    for (i = 1; i < n; i++)
        for (j = 0; j < i; j++)
            if (nmod_poly_is_zero(nmod_poly_mat_entry(tsf, n+i, n+j)) != 1) 
            {
                printf("error in nmod_poly_mat_inv_trunc: check of identity, singular at origin?\n");
                exit(-1);
            }
    
    // Extract the truncated inverse from the approximant basis 

    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(S, i, j),nmod_poly_mat_entry(tsf, n+i, j));


    nmod_poly_mat_clear(H);
    nmod_poly_mat_clear(tsf);

}



/**
 * 
 * TODO: probably only degree > 0 currently 
 *
 * x-adic iterations à la Dixon : AX=B mod x^sigma
 * A is nxn, B is nxm 
 * A is assumed to be invertible for x=0 (will produce an error otherwise) 
 * 
 * Iterations mod x^order to obtain a total approximation mod x^sigma
 * 
 */


void nmod_poly_mat_dixon(nmod_poly_mat_t X, 
                            const nmod_poly_mat_t A, 
                            const nmod_poly_mat_t B, 
                            ulong order,
                            ulong sigma)

{

    slong l;

    slong n=A->r;

    // C: inverse of A mod x^order 
    nmod_poly_mat_t C;
    nmod_poly_mat_init(C, n, n, A->modulus);

    nmod_poly_mat_inv_trunc(C, A, order);

    // Number of iterations 
    slong nbiter=ceil((double) sigma/order);


    // Temporary matrices 
    nmod_poly_mat_t BB;
    nmod_poly_mat_init(BB, n, B->c, A->modulus);
    nmod_poly_mat_set(BB,B);

    nmod_poly_mat_t S;
    nmod_poly_mat_init(S, n, B->c, A->modulus);

    nmod_poly_mat_t T;
    nmod_poly_mat_init(T, n, B->c, A->modulus);

    // The solution 
    nmod_poly_mat_zero(X);


    // printf("\n");
    // printf("sigma: %ld  order: %ld  nbiter: %ld \n",sigma,order,nbiter);
    // //nmod_poly_mat_print_pretty(A, "x");
    // printf("\n");

    // Main loop for Dixons's iterations 

    for (l=0; l<nbiter; l++) 
    {

        nmod_poly_mat_mul(S, C, BB);
        nmod_poly_mat_truncate(S,order);

        // X is constructed matrix digit by matrix digit 
        nmod_poly_mat_shift_left(T,S,l*order);
        nmod_poly_mat_add(X,X,T);

        // New residue for next iteration 
        nmod_poly_mat_mul(T, A, S);
        nmod_poly_mat_sub(T,BB,T);
        nmod_poly_mat_shift_right(BB,T,order);

    }

    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(BB);
    nmod_poly_mat_clear(S);
    nmod_poly_mat_clear(T);


}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
