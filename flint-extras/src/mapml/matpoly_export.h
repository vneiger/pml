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


#ifndef MAPML_MATPOLY_EXPORT_H
#define MAPML_MATPOLY_EXPORT_H


#include "mapml.h"
#include "conversion.h"


ALGEB pm_coeffs(MKernelVector kv, ALGEB *args);


/**********************************************************
 * 
 * modulo matrix polynomial determinant  
 * 
 *  ALGEB args[1]: matrix polynomial, vector entries 
 *        args[2]: modulus 
 * 
 * Returns the determinant as a list of coefficients 
 * 
 ***********************************************************/


ALGEB pm_determinant(MKernelVector kv, ALGEB *args);



/*******************************************************************
 * 
 * Approximant row basis for a vector (m x 1 matrix) of derivatives 
 * 
 *  ALGEB args[1]: shift
 *        args[2]: polynomial (vector of coefficients)  
 *        args[3]: order of approximation 
 *        args[4]: modulus 
 *        args[5]: method, "MBasis" or "PMBasis" 
 * 
 *  Returns ALGEB  res:=[M,dct] 
 *    M: a polynomial matrix, list entries, m x m 
 *    dct: list of out defects  
 * 
 *  Eventually, be careful with the sign between either 
 *      the defect (e.g. in gfun) or shifts = -dct in pml
 * 
 *******************************************************************/


ALGEB pm_diff_mbasis(MKernelVector kv, ALGEB *args);


/**********************************************************
 * 
 * modulo matrix polynomial mbasis  
 * 
 *  ALGEB args[1]: shift
 *        args[2]: matrix polynomial, vector entries 
 *        args[3]: order
 *        args[4]: modulus 
 * 
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

ALGEB pm_matrix_mbasis(MKernelVector kv, ALGEB *args);



/**********************************************************
 * 
 * modulo matrix polynomial pmbasis  
 * 
 *  ALGEB args[1]: shift
 *        args[2]: matrix polynomial, vector entries 
 *        args[3]: order
 *        args[4]: modulus 
 * 
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

ALGEB pm_matrix_pmbasis(MKernelVector kv, ALGEB *args);


/**********************************************************
 * 
 * modulo polynomial matrix weak Popov form   
 * 
 * ++++++++ TODO COMMENT
 * 
 *  
 * 
 ***********************************************************/


// SEE WHICH ENTRIES

 ALGEB pm_weakpopov(MKernelVector kv, ALGEB *args);

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
