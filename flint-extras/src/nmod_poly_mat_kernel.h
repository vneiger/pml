/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_KERNEL_H
#define NMOD_POLY_MAT_KERNEL_H

#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been 
 *     sorted by shifted degree 
 * 
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 * 
 *  Calls nmod_poly_mat_zls_sorted after an initial sorting 
 * 
 *  TODO/TO SEE: 
 *    
 * Input: 
 *    iA in m x n 
 *     ishift[n], NULL (the degrees are computed) or initialized outside, 
 *      the shift for the kernel    
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the 
 *       corresponding column degrees of A." 
 *    kappa, a double >= 2, for the order of the order bases 
 *              kappa * s instead of 3 *s in ZLS  
 *
 *  Output:
 *    returns the dimension w of the kernel, which may be zero 
 *    N, is initialized  n x n outside  
 *       its first w columns give a minimal basis of the kernel   
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shifted degrees of the kernel basis 
 * 
 */

int nmod_poly_mat_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                         const slong *ishift, const double kappa); 



/**
 * Experimental, should not be really considered  
 *
 */

int nmod_poly_mat_approximant_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                                 const slong *ishift);


#ifdef __cplusplus
}
#endif

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


