/*
    Copyright (C) 2025 Gilles Villard
    Copyright (C) 2026 Vincent Neiger

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
 * whose rank is nullity = m - rank(F), called the nullity of F. Such a basis
 * matrix has dimensions nullity x m and has full row rank.
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
 * \return slong, the nullity of `pmat`
 *
 * \param[out] ker the output kernel basis (cannot alias `pmat`, must be
 * initialized with at least `nullity` rows)
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
 * first `nullity` rows contain useful entries. Similarly, for `pivind` and `shift`
 * only the first `nullity` entries matter in output.
 */

/** general interface
 * -> all four orientations are supported (ROW_LOWER|ROW_UPPER yield a left
 *  kernel basis while COL_LOWER|COL_UPPER yield a right kernel basis)
 * -> pivind may be NULL if that information is not wanted
 * -> shift may be NULL in input, in which case it is left NULL in output, and
 *  the basis is computed in weak Popov form with respect to row degrees of
 *  `pmat` (if left kernel) or column degrees of `pmat` (if right kernel)
 * -> size of the integers (see @ref IntegerSizes): this relies on
 *  ::nmod_poly_mat_kernel_zls_approx, whose condition is assumed, with `m`
 *  being the number of rows of `pmat` for a left kernel, and its number of
 *  columns for a right kernel. When `shift` is NULL, the shift used is the
 *  row or column degree of `pmat`, and this condition only involves degrees.
 */
/* TODO `form` not implemented yet, only shifted weak Popov */
slong nmod_poly_mat_kernel(nmod_poly_mat_t ker,
                           slong * pivind,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           poly_mat_form_t form,
                           orientation_t orient);

/** Computes a `shift`-ordered weak Popov left kernel basis `ker` for `pmat`,
 * through a minimal approximant basis at sufficiently large order.
 *
 * The approximation order is `min(m,n) * deg(pmat) + max(shift) - min(shift)
 * + deg(pmat) + 1`, where `pmat` is `m x n`: the cost grows with the
 * amplitude `max(shift) - min(shift)` of the shift.
 *
 * Size of the integers (see @ref IntegerSizes): the shift entries then
 * increase by at most this order in the approximant basis computation. It is
 * assumed that `|shift[i]| <= WORD_MAX / 8` for all `i` and that
 * `(min(m,n) + 1) * deg(pmat) + 1 <= WORD_MAX / 8`; this ensures that no
 * overflow occurs.
 */
slong nmod_poly_mat_kernel_via_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      const nmod_poly_mat_t pmat);

/** Computes a `shift`-ordered weak Popov left kernel basis `ker` for `pmat`,
 * using the Zhou-Labahn-Storjohann algorithm:
 *   Wei Zhou, George Labahn, Arne Storjohann -- "Computing Minimal Nullspace Bases" 
 *   Proceedings ISSAC 2012 -- https://doi.org/10.1145/2442829.2442881
 *
 * \todo unbalanced row degrees (e.g. a few rows of large degree) are partly
 * handled (in particular with linear system solving instances in mind, one row
 * of very large degree), but we expect some improvements to still be possible
 *
 * \todo unbalanced column degrees (e.g. a few columns of large degree) are not
 * handled specifically
 */
/** Size of the integers (see @ref IntegerSizes): at each step, the
 * approximation order is computed from the sum of the shift entries (and `m`
 * times their minimum, up to degrees), where `m` is the number of rows of
 * `pmat`; then the shift entries increase by the degrees of the computed
 * approximant bases, and are used as shifts in the subsequent recursive calls.
 * It is assumed that `|shift[i]| <= WORD_MAX / (8 m)` for all `i`, and that
 * `m` times the degrees of the polynomial matrices computed along the way
 * (`pmat`, intermediate approximant bases, the output kernel basis) is at most
 * `WORD_MAX / 8`; this ensures that no overflow occurs.
 */
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


/** @name Kernel of a Krylov matrix
 *
 * \todo this is not at the right place: where should it go?
 * a general relation basis header?
 */
//@{

/** Sets `P` to the basis in `shifts`-Popov form (row-wise, with ROW_LOWER
 * pivots, see @ref MatrixForms), which is unique, of the module of vectors
 * `p` of `m` polynomials such that `E_0 p_0(M) + ... + E_{m-1} p_{m-1}(M) =
 * 0`, where `E` is `m x n` and `M` is `n x n` (see @ref KrylovMatrices).
 * The output `P` must be initialized with dimensions `m x m`. Returns 1, or
 * 0 if the degree bounds `degrees` are found not to be valid (`P` is then
 * undefined). The other arguments and the conditions on them are as in
 * ::nmod_mat_krylov_kernel, which gives the same basis in compact form.
 */
int nmod_poly_mat_krylov_kernel_basis(nmod_poly_mat_t P,
                                      const nmod_mat_t E,
                                      const nmod_mat_t M,
                                      const slong * shifts,
                                      const slong * degrees);

/** Same as ::nmod_poly_mat_krylov_kernel_basis, for a Jordan matrix `M`
 * given by `eig`, `sizes` and `nblocks` (as in
 * ::nmod_mat_krylov_kernel_jordan). */
int nmod_poly_mat_krylov_kernel_basis_jordan(nmod_poly_mat_t P,
                                             const nmod_mat_t E,
                                             nn_srcptr eig,
                                             const slong * sizes,
                                             slong nblocks,
                                             const slong * shifts,
                                             const slong * degrees);

//@} // doxygen group: Kernel of a Krylov matrix

#ifdef __cplusplus
}
#endif

#endif
