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

/** \file nmod_poly_mat_kernel.h
 * Definition (shifted minimal kernel basis).
 * ------------------------------------------
 * Consider:
 *   - an m x n matrix of univariate polynomials F,
 *   - a degree shift s (a list of m integers).
 *
 * A (left) kernel basis for F is a matrix over the univariate polynomials
 * whose rows form a basis for the module
 *   { p in K[X]^{1 x m}  |  p * F == 0 },
 * whose rank is nz = m - rank(F), called the nullity of F. Such a basis
 * matrix has dimensions nz x m and has full row rank.
 *
 * A kernel basis for (F,d) is said to be <em>a shift-minimal</em> (resp.
 * <em>a shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>)
 * kernel basis if it is in shift-reduced form (resp. in shift-ordered weak
 * Popov form, resp. in shift-Popov form). See nmod_poly_mat_forms.h for
 * definitions of these forms.
 */

/** \file nmod_poly_mat_kernel.h
 * Conventions.
 * ------------
 * Apart from the general interfaces (TODO) which offer several choices of
 * orientation, all other functions compute left kernel bases and use the
 * following parameters:
 *
 * \param[out] ker the output kernel basis (cannot alias `pmat`)
 * \param[in] pmat the input polynomial matrix (no restriction)
 * \param[in,out] shift in: the input shift; and out: the output shifted row
 * degree of `ker` (list of integers, length must be the number of rows of
 * `pmat`)
 *
 * The computed `ker` is in shifted ordered weak Popov form, or in the canonical
 * shifted Popov form when the name of the function indicates so.
 */









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

/* int nmod_poly_mat_kernel(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \ */
/*                          const slong *ishift, const double kappa); */

int nmod_poly_mat_kernel_zls(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
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
