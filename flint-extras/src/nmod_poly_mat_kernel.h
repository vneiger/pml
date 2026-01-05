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

#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_forms.h"

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
 * \anchor kernel_basis
 *
 * Apart from the general interfaces which offer several choices of
 * orientation, all other functions compute left kernel bases and use the
 * following parameters:
 *
 * \return slong, the nullity of `pmat` (called `nz` below)
 *
 * \param[out] ker the output kernel basis (cannot alias `pmat`, must be
 * initialized with at least `nz` rows)
 * \param[out] pivind the pivot index of `ker` (list of integers, length must
 * be the number of rows of `pmat`)
 * \param[in,out] shift in: the input shift; and out: the output shifted row
 * degree of `ker` (list of integers, length must be the number of rows of
 * `pmat`)
 * \param[in] pmat the input polynomial matrix (no restriction)
 *
 * The computed `ker` is in shifted ordered weak Popov form, or in the
 * canonical shifted Popov form when the name of the function indicates so.
 *
 * Some functions may overwrite data in pmat (check the presence of `const`).
 *
 * `ker` is not resized/reallocated: it keeps its original shape, but only the
 * first `nz` rows contain useful entries. Similarly, for `pivind` and `shift`
 * only the first `nz` entries matter in output.
 */

/** general interface */
/* TODO `form` not really implemented yet, only shifted weak Popov */
slong nmod_poly_mat_kernel(nmod_poly_mat_t ker,
                           slong * pivind,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           poly_mat_form_t form,
                           orientation_t orient);

/** Computes a `shift`-ordered weak Popov left kernel basis `ker` for `pmat`,
 * through a minimal approximant basis at sufficiently large order.
 */
slong nmod_poly_mat_kernel_via_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      const nmod_poly_mat_t pmat);

/** Computes a `shift`-ordered weak Popov left kernel basis `ker` for `pmat`,
 * using the Zhou-Labahn-Storjohann algorithm:
 *   Wei Zhou, George Labahn, Arne Storjohann -- "Computing Minimal Nullspace Bases" 
 *   Proceedings ISSAC 2012 -- https://doi.org/10.1145/2442829.2442881
 */
/* TODO middle_product currently naive */
slong nmod_poly_mat_kernel_zls_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      nmod_poly_mat_t pmat);


/* FIXME not implemented yet: using interpolant basis at sufficiently many points */
/* slong nmod_poly_mat_kernel_via_interp(nmod_poly_mat_t ker, */
/*                                       slong * shift, */
/*                                       const nmod_poly_mat_t pmat); */

/* FIXME not implemented yet: ZLS algo using interpolant basis */
/* slong nmod_poly_mat_kernel_zls_interp(nmod_poly_mat_t ker, */
/*                                       slong * shift, */
/*                                       const nmod_poly_mat_t pmat); */

/* FIXME could add a "prepare pmat" function:
 * A/ compute rdeg and cdeg
 * ->suppress zero columns
 * ->trivialize zero rows
 * B/ handle case of constant matrices
 * C/ go for actual work
 * */



/** Verifying if a matrix is a minimal kernel basis.
 *
 * This checks whether the matrix `ker` is a `shift`-minimal kernel basis for
 * `pmat` for the required form `form`, with orientation described by `orient`.
 *
 * \param[in] ker kernel basis
 * \param[in] shift shift
 * \param[in] pmat polynomial matrix
 * \param[in] form required form for `kerbas` (see #poly_mat_form_t)
 * \param[in] orient indicates the orientation (left/right kernel) and the definition of pivots
 * \param[in] randomized [not implemented yet] if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return int, result of the verification
 *
 */
/* TODO WARNING! for the moment, does not really check generation! */
int nmod_poly_mat_is_kernel(const nmod_poly_mat_t ker,
                            const slong * shift,
                            const nmod_poly_mat_t pmat,
                            poly_mat_form_t form,
                            orientation_t orient);


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
