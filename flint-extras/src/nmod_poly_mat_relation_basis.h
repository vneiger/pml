/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_RELATION_BASIS_H
#define NMOD_POLY_MAT_RELATION_BASIS_H

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/** \file nmod_poly_mat_relation_basis.h
 * Relation bases.
 * ---------------
 * Given an `m x n` matrix `E` with rows `E_0, ..., E_{m-1}` and an `n x n`
 * matrix `M`, the vectors `p` of `m` polynomials such that
 * `E_0 p_0(M) + ... + E_{m-1} p_{m-1}(M) = 0` form a module of rank `m`,
 * the module of relations of `E` with respect to `M`. Its bases are the
 * kernel bases of the Krylov matrix of `E` and `M` (see @ref KrylovMatrices).
 *
 * Approximant, M-Pade and interpolant bases are relation bases with respect
 * to Jordan matrices (nilpotent blocks for approximants, a diagonal matrix
 * for interpolants): see ::nmod_poly_mat_approximant_basis_krylov and
 * ::nmod_poly_mat_mpade_basis_krylov in nmod_poly_mat_approximant.h, and
 * ::nmod_poly_mat_interpolant_basis_krylov in nmod_poly_mat_interpolant.h.
 */

/** @name Relation bases via Krylov matrices
 */
//@{

/** Sets `P` to the basis in `shift`-Popov form (row-wise, with ROW_LOWER
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
                                      const slong * shift,
                                      const slong * degrees);

/** Same as ::nmod_poly_mat_krylov_kernel_basis, for a Jordan matrix `M`
 * given by `eig`, `sizes` and `nblocks` (as in
 * ::nmod_mat_krylov_kernel_jordan). */
int nmod_poly_mat_krylov_kernel_basis_jordan(nmod_poly_mat_t P,
                                             const nmod_mat_t E,
                                             nn_srcptr eig,
                                             const slong * sizes,
                                             slong nblocks,
                                             const slong * shift,
                                             const slong * degrees);

//@} // doxygen group: Relation bases via Krylov matrices

#ifdef __cplusplus
}
#endif

#endif
