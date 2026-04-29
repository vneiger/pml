/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h" // for nmod_find_root

#include "nmod_poly_mat_multiply.h"



/** Spreads the columns of A into chunks of length 'chunk' 
 *   for each column j of A there are colAL[j] new created columns,
 *   lin_A has been initialized outside with the correct dimensions
 */


void spread_columns(nmod_poly_mat_t lin_A, const slong* colAL,\
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
}

           
/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  
 *  TODO outputt can alias input 
 *  TODO ASSUMPTION (not checked): existence of element of "large enough" order
 *  TODO -> fail flag when element not found
 *  
 */





