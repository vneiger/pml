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


#ifndef MAPML_CONVERSION_H
#define MAPML_CONVERSION_H

#include "mapml.h"



/**********************************************************
 * 
 * Converts an nmod_poly_t to a maple list of coefficients 
 *   !! no modulus is transmitted
 * 
 *********************************************************/


ALGEB nmod_poly_to_algeb(MKernelVector kv, const nmod_poly_t p); 


/**********************************************************
 * 
 * Converts an nmod_mat_poly to a maple matrix of lists  
 *   !! no modulus is transmitted
 * 
 * Todo: convert matrix coefficients directly, 
 *    not by global matrix conversion ? 
 * 
 * For the moment: simply goes via a poly_mat 
 * 
 ***********************************************************/

ALGEB nmod_mat_poly_to_algeb(MKernelVector kv, const nmod_mat_poly_t Ain);


/**********************************************************
 * 
 * Converts an nmod_poly_mat to a maple matrix of lists  
 *   !! no modulus is transmitted
 * 
 ***********************************************************/

ALGEB nmod_poly_mat_to_algeb(MKernelVector kv, const nmod_poly_mat_t A);


/************************************************************
 * 
 * Converts a maple coefficient vector 
 *  to an nmod_poly_t   
 * 
 *  ALGEB vect_A: a maple vector
 * 
 *  Initialization included
 * 
 *************************************************************/


void get_nmod_poly(nmod_poly_t p,   const ulong modulus, MKernelVector kv, ALGEB vect);


/************************************************************
 * 
 * Converts a maple matrix of coefficient vectors  
 *  to an nmod_mat_poly_t   
 * 
 *  ALGEB vect_A: a maple matrix of vectors
 *  
 *  Initialization included
 * 
 *************************************************************/


void get_nmod_mat_poly(nmod_mat_poly_t Aout,   const ulong modulus, MKernelVector kv, ALGEB vect_A);


/************************************************************
 * 
 * Converts a maple matrix of coefficient vectors  
 *  to an nmod_poly_mat_t   
 * 
 *  ALGEB vect_A: a maple matrix of vectors
 *  
 *  Initialization included
 * 
 *************************************************************/

void get_nmod_poly_mat(nmod_poly_mat_t A,   const ulong modulus, MKernelVector kv, ALGEB vect_A);


void get_fmpq_poly(fmpq_poly_t p, MKernelVector kv, ALGEB vect);

void get_fmpq_poly_array(fmpq_poly_t *vp, MKernelVector kv, ALGEB maple);

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
