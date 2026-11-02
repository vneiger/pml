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

#ifndef NMOD_POLY_MAT_LINSOLVE_H
#define NMOD_POLY_MAT_LINSOLVE_H

/** \file nmod_poly_mat_linsolve.h
 *
 * Solving linear systems with univariate polynomial matrices: power series
 * solution (truncated at some order), and solution over the polynomials (with
 * a denominator).
 *
 * \todo series solving at low precision, high-order lifting, solving via
 * series and rational reconstruction (see ntl-extras, mat_lzz_pX_linsolve.h)
 */

#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                  POWER SERIES SOLUTION                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Power series solution
 * \anchor SeriesSolve
 *
 * For an `m x m` polynomial matrix `A` whose constant coefficient `A(0)` is
 * invertible, an `n x m` polynomial matrix `B`, and a nonnegative integer
 * `prec`, the power series solution of `X A = B` at precision `prec` is the
 * unique `n x m` polynomial matrix `X` of degree less than `prec` such that
 * `X A = B mod x^prec`; in other words, `X = B A^{-1} mod x^prec`, where
 * `A^{-1}` is the inverse of `A` in the ring of matrices over the power
 * series.
 */
//@{

/** Computes the power series solution `X = B A^{-1} mod x^prec` (see @ref
 * SeriesSolve) by x-adic lifting à la Dixon: the truncated inverse of `A` at
 * order `order` is computed by ::nmod_poly_mat_inv_series, and then each
 * iteration finds `order` new coefficients of `X`, by a low product with this
 * truncated inverse and a middle product with `A` (to update the residue).
 * The residue is kept of length less than `length(A)`, and each iteration
 * reads only the `order` coefficients of `B` it needs, so that the data
 * movements are linear in `prec` even if `B` has large length.
 *
 *  - `X` must be initialized with the same dimensions as `B`; it may alias
 *  `A` or `B`;
 *  - the return value is `1` if `A(0)` is invertible (or if `prec` is zero,
 *  or if `m` is zero), and `0` otherwise, in which case `X` is not modified;
 *  - if `order` is nonpositive, the default `order = length(A) = deg(A) + 1`
 *  is used (or `1` if `A` is zero); in any case, `order` is reduced to `prec`
 *  if it is larger;
 *  - only the coefficients of `A` and `B` of degree less than `prec` are used;
 *  - nonpositive values of `prec` are considered as `0`.
 *
 * Size of the integers (see @ref IntegerSizes): `min(order, prec)` is
 * assumed to be at most `WORD_MAX / 2` (see ::nmod_poly_mat_inv_series).
 *
 * \todo choice of `order` and of the algorithm: tuning (the iterations are
 * cheaper with a longer truncated inverse if `A` has small degree, while for
 * `n` close to `m` and `prec` not much larger than `deg(A)`, it is faster to
 * compute the truncated inverse of `A` at order `prec` and multiply by `B`)
 * \todo precomputed transforms of `A` and of its truncated inverse, shared by
 * all iterations
 */
int nmod_poly_mat_solve_series_dixon(nmod_poly_mat_t X,
                                     const nmod_poly_mat_t A,
                                     const nmod_poly_mat_t B,
                                     slong order,
                                     slong prec);

//@} // doxygen group: Power series solution

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                SOLUTION OVER THE POLYNOMIALS               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Solution over the polynomials
 * \anchor PolySolve
 *
 * For an `m x n` polynomial matrix `A` and a `1 x n` polynomial vector `b`,
 * the functions below compute a `1 x m` polynomial vector `u` and a nonzero
 * polynomial `den` such that `u A = den b`, that is, a solution `u / den` of
 * `y A = b` over the rational fractions, if such a solution exists (that is,
 * if `b` is in the row space of `A` over the rational fractions).
 */
//@{

/** Computes `u` and `den` such that `u A = den b` (see @ref PolySolve), via a
 * kernel basis. Returns `1` if a solution exists, and `0` otherwise, in which
 * case `u` and `den` are not modified.
 *
 * The polynomial vector `u` must be initialized with `1` row and `m` columns;
 * it may alias `A` or `b`. The vector `b` must have `1` row (not checked).
 *
 * The output `den` is monic, and `den` and the entries of `u` have no
 * nontrivial common divisor. If `A` has full row rank (for example, if `A` is
 * square and nonsingular), this specifies a unique solution: `u / den` is the
 * solution of `y A = b` written in lowest terms, with a common denominator.
 *
 * Algorithm: a vector `[-u, den]` with `den` nonzero is sought in the left
 * kernel of the `(m+1) x n` matrix `[A ; b]` (`A` with `b` appended as last
 * row). A kernel basis is computed by ::nmod_poly_mat_kernel_zls_approx, with
 * the shift `s` being the row degree of this matrix, that is, the row degrees
 * of `A` and the degree of `b` (zero rows counting as degree `0`). A solution
 * exists if and only if some row of this basis has a nonzero last entry, and
 * the returned solution is the normalization of such a row with minimal
 * `s`-row degree. Since this basis is `s`-reduced, the predictable degree
 * property shows that this solution has minimal `s`-row degree among all
 * solutions: it minimizes `max(max_i (deg(u_i) + s_i), deg(den) + s_m)`.
 * This is relevant only when `A` does not have full row rank, and in general,
 * several solutions reach this minimal degree.
 *
 * Note: the solution is not necessarily the row of this basis whose `s`-pivot
 * index is the last one; such a row may not exist. For example, for `A = [[x,
 * 1], [x^2 + 1, x]]` and `b = [1, 0]`, the kernel is generated by `[x, -1,
 * 1]`, whose `s`-row degree `2` (for `s = [1, 2, 0]`) is reached at indices
 * `0` and `1`, so that its `s`-pivot index is `1`.
 *
 * Size of the integers (see @ref IntegerSizes): the condition of
 * ::nmod_poly_mat_kernel_zls_approx is assumed, for the matrix `[A ; b]` and
 * the shift `s` described above.
 */
int nmod_poly_mat_linsolve_via_kernel(nmod_poly_mat_t u,
                                      nmod_poly_t den,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t b);

//@} // doxygen group: Solution over the polynomials

#ifdef __cplusplus
}
#endif

#endif  /* NMOD_POLY_MAT_LINSOLVE_H */
