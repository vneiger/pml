/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/
#include <math.h>

#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h" // for nmod_find_root

#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_forms.h"



/** Spreads the columns of A into chunks of length 'chunk' 
 *   for each column j of A there are colAL[j] new created columns,
 *   lin_A has been initialized outside with the correct dimensions
 */


void static spread_columns(nmod_poly_mat_t lin_A, const slong* colAL,\
                const nmod_poly_mat_t A, const slong chunk)
{
    int i,j,l;

    slong m = A->r;
    slong n = A->c; 

    nmod_poly_t tpol;
    nmod_poly_init(tpol,A->modulus);

    slong dec=0;

    /** on the columns of A */
    for (j=0; j<n; j++) 
    {
       /** to see: loop in another way for efficiency */
       if (colAL[j] > 0)
       {
            for (i=0; i<m; i++) 
            {
                nmod_poly_set(tpol, nmod_poly_mat_entry(A, i, j));
                nmod_poly_set_trunc(nmod_poly_mat_entry(lin_A, i, dec),tpol,chunk);

                for (l=1; l<colAL[j]; l++) 
                {
                    if (chunk <= nmod_poly_degree(tpol)+1)
                    {
                        nmod_poly_shift_right(tpol,tpol,chunk);
                        nmod_poly_set_trunc(nmod_poly_mat_entry(lin_A, i, dec+l),tpol,chunk);
                    }
                    else 
                    {
                        nmod_poly_zero(nmod_poly_mat_entry(lin_A, i, dec+l));
                    }
                }
            }
            dec+=colAL[j];
        }
    }

    nmod_poly_clear(tpol);
}

           
/** The columns of lin_A are compressed to give A  
 *    colAL[j] chunks of length 'chunk' give the column j in A 
 */

void static compress_columns(nmod_poly_mat_t A, const slong* colAL,\
                const nmod_poly_mat_t lin_A, const slong chunk)
{
    int i,j,l;

    slong m = A->r;
    slong n = A->c; 

    nmod_poly_t tpol,taccu;
    nmod_poly_init(tpol,A->modulus);
    nmod_poly_init(taccu,A->modulus);

    slong dec=0;

    nmod_poly_mat_zero(A);

    for (j=0; j<n; j++) 
    {
       /** to see: loop in another way for efficiency */
       if (colAL[j] > 0)
       {
            for (i=0; i<m; i++) 
            {
                nmod_poly_set(taccu, nmod_poly_mat_entry(lin_A, i, dec));

                for (l=1; l<colAL[j]; l++) 
                {
                    nmod_poly_shift_left(tpol,nmod_poly_mat_entry(lin_A, i, dec+l), l*chunk);
                    nmod_poly_add(taccu, taccu, tpol);
                }

                nmod_poly_set(nmod_poly_mat_entry(A, i, j), taccu);
            }

            dec+=colAL[j];
        }
    }

    nmod_poly_clear(tpol);
    nmod_poly_clear(taccu);
}


/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input 
 *  uses mul_geometric 
 * 
 *  Linearizes B by columns into chunks of length deg(A) + 1
 * 
 *  Todo ASSUMPTION (not checked): existence of element of "large enough" order
 *           and fail flag when element not found
 *  Todo  no linearization in certain cases ?  
 *  
 */

void nmod_poly_mat_mul_linearized(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    slong m = A->r;
    slong l = A->c;
    slong n = B->c;

    if (m < 1 || (A->c) < 1 || n < 1 || nmod_poly_mat_is_zero(B))
    {
        nmod_poly_mat_zero(C);
        return;
    }

    /** Really needed for aliasing?
     */
    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, A->modulus);
        nmod_poly_mat_mul_linearized(T, A, B);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    /** the chunk could be adapted 
     */
    slong chunk; 
    chunk = (nmod_poly_mat_degree(A)+1);

    if ((nmod_poly_mat_degree(B)+1) <= chunk)
    {  
        nmod_poly_mat_mul(C,A,B);
        return;
    }

    slong dB[n];
    nmod_poly_mat_column_degree(dB, B, NULL);

    slong colBL[m];  // column j of B gives colBL[j] columns in output  
    slong ln=0;   // new column dimension 

    for (int j=0; j<n; j++)
    {
        colBL[j] = ceil((double) (dB[j]+1)/chunk);
        ln += colBL[j];
    }

    nmod_poly_mat_t lin_B;
    nmod_poly_mat_init(lin_B, l, ln, A->modulus);

    spread_columns(lin_B, colBL, B, chunk); 

    nmod_poly_mat_t lin_C;
    nmod_poly_mat_init(lin_C, m, ln, A->modulus);

    nmod_poly_mat_mul_geometric(lin_C,A,lin_B);

    compress_columns(C, colBL, lin_C, chunk); 

    nmod_poly_mat_clear(lin_B);
    nmod_poly_mat_clear(lin_C);
}



