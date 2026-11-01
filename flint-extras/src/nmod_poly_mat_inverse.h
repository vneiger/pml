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

#ifndef NMOD_POLY_MAT_INVERSE_H
#define NMOD_POLY_MAT_INVERSE_H

/** \file nmod_poly_mat_inverse.h
 *
 * Inverse of a univariate polynomial matrix, and related operations such as
 * the truncated inverse (expansion of the inverse as a power series).
 *
 * \todo implementation of inversion algorithms (see the corresponding todo in
 * ntl-extras, mat_lzz_pX_inverse.h)
 * \todo high-order lifting (see ntl-extras, mat_lzz_pX_inverse.h)
 */

#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                TRUNCATED INVERSE EXPANSION                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Truncated inverse expansion
 * \anchor TruncatedInverse
 *
 * For a square polynomial matrix `pmat` and a nonnegative integer `order`,
 * the truncated inverse of `pmat` at order `order` is the unique polynomial
 * matrix `imat` of degree less than `order` such that `imat * pmat` and
 * `pmat * imat` are both the identity matrix modulo `x^order`; in other
 * words, `imat` is the inverse of `pmat` in the ring of matrices over the
 * power series, truncated at order `order`. It exists if and only if the
 * constant coefficient `pmat(0)` is invertible (or `order` is zero, in which
 * case it is zero).
 *
 * In the functions below:
 *  - `pmat` must be square, and `imat` must be initialized with the same
 *  dimensions and modulus as `pmat`;
 *  - the functions which take as input the inverse `inv0` of the constant
 *  coefficient of `pmat` (constant matrix of the same dimensions) have no
 *  failure case, and `imat` cannot alias `pmat` in these functions; for the
 *  other functions, `imat` may alias `pmat`, and the return value is `1` if
 *  `pmat(0)` is invertible (or if `order` is zero, or if `pmat` is the empty
 *  matrix), and `0` otherwise, in which case `imat` is not modified;
 *  - only the coefficients of `pmat` of degree less than `order` are used;
 *  - nonpositive values of `order` are considered as `0`.
 *
 * Size of the integers (see @ref IntegerSizes): `order` is assumed to be at
 * most `WORD_MAX / 2`.
 */
//@{

/** Threshold on truncation order (switching from basecase to Newton).
 * \todo provisional values, to be tuned (also depending on the modulus) */
slong _nmod_poly_mat_inv_series_cutoff(slong m);

/** Threshold on the length of the input: if `pmat` has length at most this,
 * the general interface uses the basecase algorithm, whatever the order.
 * \todo provisional values, to be tuned (also depending on the modulus) */
slong _nmod_poly_mat_inv_series_short_length(slong m);

/** Computes the truncated inverse `imat` of `pmat` at order `order` (see
 * @ref TruncatedInverse). This is the main interface, which chooses between
 * basecase and Newton, depending on `order` and on the degree of `pmat` (see
 * ::_nmod_poly_mat_inv_series_cutoff and
 * ::_nmod_poly_mat_inv_series_short_length). */
int nmod_poly_mat_inv_series(nmod_poly_mat_t imat,
                             const nmod_poly_mat_t pmat,
                             slong order);

/** Same as ::nmod_poly_mat_inv_series, with `inv0` the inverse of `pmat(0)`
 * given as input (see @ref TruncatedInverse). */
void _nmod_poly_mat_inv_series(nmod_poly_mat_t imat,
                               const nmod_poly_mat_t pmat,
                               const nmod_mat_t inv0,
                               slong order);

/** Computes the truncated inverse `imat` of `pmat` at order `order` (see
 * @ref TruncatedInverse) by the basecase algorithm. The cost is quadratic in
 * `order` if `deg(pmat) >= order`, and linear in `order` for a fixed degree of
 * `pmat`. */
int nmod_poly_mat_inv_series_basecase(nmod_poly_mat_t imat,
                                      const nmod_poly_mat_t pmat,
                                      slong order);

/** Same as ::nmod_poly_mat_inv_series_basecase, with `inv0` the inverse of
 * `pmat(0)` given as input (see @ref TruncatedInverse). */
void _nmod_poly_mat_inv_series_basecase(nmod_poly_mat_t imat,
                                        const nmod_poly_mat_t pmat,
                                        const nmod_mat_t inv0,
                                        slong order);

/** Computes the truncated inverse `imat` of `pmat` at order `order` (see
 * @ref TruncatedInverse) by Newton iteration, whose cost is quasi-linear in
 * `order`. If `cutoff` is nonpositive, the default
 * ::_nmod_poly_mat_inv_series_cutoff is used; values `1` and `2` are
 * equivalent (Newton iteration starting from `pmat(0)^{-1}`). */
int nmod_poly_mat_inv_series_newton(nmod_poly_mat_t imat,
                                    const nmod_poly_mat_t pmat,
                                    slong order,
                                    slong cutoff);

/** Same as ::nmod_poly_mat_inv_series_newton, with `inv0` the inverse of
 * `pmat(0)` given as input (see @ref TruncatedInverse), and `cutoff` used as
 * given (values less than `2` are considered as `2`). */
void _nmod_poly_mat_inv_series_newton(nmod_poly_mat_t imat,
                                      const nmod_poly_mat_t pmat,
                                      const nmod_mat_t inv0,
                                      slong order,
                                      slong cutoff);

/** Truncated inverse of `A` modulo `x^order`, via a shifted approximant
 * basis of `[A; -I]`.
 *
 * Returns `1` if the bottom-right block of the computed approximant basis is
 * the identity matrix, in which case `S` is the truncated inverse (the
 * bottom-left block of that basis); otherwise returns `0`, and `S` is not
 * modified. In particular, `0` is returned if `A(0)` is not invertible. This
 * function is kept for comparison with ::nmod_poly_mat_inv_series, which is
 * faster.
 *
 * CHECK TODO: check whether the specification of nmod_poly_mat_pmbasis
 * ensures that the bottom-right block is the identity when `A(0)` is
 * invertible
 *
 * Size of the integers (see @ref IntegerSizes): the shift used in the
 * approximant basis computation has entries `0` and `order`, which then
 * increase by at most `order`; it is assumed that `order <= WORD_MAX / 2`.
 */
slong nmod_poly_mat_inv_trunc(nmod_poly_mat_t S,
                              const nmod_poly_mat_t A,
                              ulong order);

//@} // doxygen group: Truncated inverse expansion

#ifdef __cplusplus
}
#endif

#endif /* NMOD_POLY_MAT_INVERSE_H */
