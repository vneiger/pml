/* 
   Copyright (C) 2024 Gilles Villard

   This file is part of mapml. mapml is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   mapml is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with mapml. If not, see <http://www.gnu.org/licenses/>. */

// *******************************************************


#ifndef MAPML_DUMMY_EXPORT_C
#define MAPML_DUMMY_EXPORT_C


#include "mapml.h"
#include "conversion.h"



/**********************************************************
 * 
 * maple polynomial round trip (for cheks) 
 * 
 * * Converts a maple vector 
 *   to a maple list polynomial
 *
 *  ALGEB args[1]: polynomial string 
 *        args[2]: modulus 
 * 
 ***********************************************************/


ALGEB polynomial_rt(MKernelVector kv, ALGEB *args){

    ALGEB vect=args[1];

    ulong modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_t p;  
    // The polynomial is initialized by the call below

    get_nmod_poly(p, modulus, kv, vect);
    
    return nmod_poly_to_algeb(kv,p);

}

/**********************************************************
 * 
 * Maple matrix polynomial round trip (for cheks) 
 * 
 * Converts a maple vector polynomial matrix 
 *   to a maple list polynomial matrix 
 *
 * Internal: mat_poly
 *  
 *  ALGEB args[1]: polynomial matrix: matrix of vectors 
 *        args[2]: modulus 
 * 
 * 
 ***********************************************************/


ALGEB matpoly_rt(MKernelVector kv, ALGEB *args){

    ALGEB vect_mat=args[1];

    ulong modulus = MapleToInteger64(kv,args[2]);

    nmod_mat_poly_t A;

    get_nmod_mat_poly(A, modulus, kv, vect_mat);
    
    return nmod_mat_poly_to_algeb(kv,A);

}


/**********************************************************
 * 
 * Maple matrix polynomial round trip (for cheks) 
 * 
 * Converts a maple vector polynomial matrix 
 *   to a maple list polynomial matrix 
 *
 * Internal: poly_mat
 *  
 *  ALGEB args[1]: polynomial matrix: matrix of vectors 
 *        args[2]: modulus 
 * 
 * 
 ***********************************************************/

ALGEB polymat_rt(MKernelVector kv, ALGEB *args){

    ALGEB vect_mat=args[1];

    ulong modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_mat_t A;

    get_nmod_poly_mat(A, modulus, kv, vect_mat);
    
    return nmod_poly_mat_to_algeb(kv,A);

}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
