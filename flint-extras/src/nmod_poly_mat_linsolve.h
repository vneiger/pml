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
 * Orientation: the systems below are row-wise, that is, `X A = B` or `u A =
 * den b` with the unknown on the left of `A`, as for the left kernel and
 * approximant bases elsewhere in PML. This differs from FLINT's solving
 * functions (`nmod_mat_solve`, `nmod_poly_mat_solve`, ...), which solve `A X =
 * B`. 
 *
 * \todo Warning: these solving functions are not finalized yet, and this
 * convention (as well as their interfaces) may change in the future.
 *
 * \todo series solving at low precision, high-order lifting (see ntl-extras,
 * mat_lzz_pX_linsolve.h); general interface for the solution over the
 * polynomials, choosing among the algorithms below
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
 * The truncated inverse and `A` are the right operands of all the products:
 * except for small degrees and dimensions (see the code for the rule), their
 * fft_small transforms are computed once and stored (see
 * ::nmod_poly_mat_mulmid_precomp_right_init), so that an iteration only transforms
 * the left operands and the results. The default order is then short, about
 * `deg(A)` for one row in `B` (for one row, the stored transforms make each
 * iteration cost about `m^2` pointwise products per coefficient of `X`, and a
 * longer order only makes the truncated inverse more expensive), longer for
 * more rows. Without stored transforms, the default order is `8 length(A)`.
 *
 *  - `X` must be initialized with the same dimensions as `B`; it may alias
 *  `A` or `B`;
 *  - the return value is `1` if `A(0)` is invertible (or if `prec` is zero,
 *  or if `m` is zero), and `0` otherwise, in which case `X` is not modified;
 *  - if `order` is nonpositive, the default order is used: `N - w + 1` with
 *  stored transforms, where `N` is the smallest power of two at least
 *  `max(16, 2 f deg(A))` (`f = 1` for one row in `B`, 2 up to 4 rows, 4 up to
 *  16 rows, 8 beyond) and `w = max(deg(A), min(length(B), (N+1)/2))`, so that
 *  the transforms of the low product have length `N`; and `8 length(A)` (or
 *  `1` if `A` is zero) otherwise; in any case, `order` is reduced to `prec`
 *  if it is larger;
 *  - only the coefficients of `A` and `B` of degree less than `prec` are used;
 *  - nonpositive values of `prec` are considered as `0`.
 *
 * Size of the integers (see @ref IntegerSizes): `min(order, prec)` is
 * assumed to be at most `WORD_MAX / 2` (see ::nmod_poly_mat_inv_series).
 *
 * Memory: the stored transforms take about `np` times the size of the
 * truncated inverse and of `A` at the transform length (np = 1 to 3 primes),
 * that is, a few times the size of `A` with the default order.
 *
 * \todo tuning, on several machines and moduli: the rule for stored
 * transforms and the default orders were measured on one machine (Intel Xeon
 * with AVX-512); on it, against the version without stored transforms (and
 * its default order `8 length(A)`), for one row and `prec = 2 m deg(A) + 1`:
 * 1.5 to 6 times faster from `deg(A)^2 m` of a few thousands (e.g. `m = 64`,
 * `deg(A) = 32`, 60-bit prime: 0.18 s against 1.24 s), 1.2 to 2.2 times
 * faster for 2 to 64 rows, 2 to 6 times faster for `B` of length `prec`.
 * Besides, for `n` close to `m` and `prec` not much larger than `deg(A)`, it
 * may be faster to compute the truncated inverse of `A` at order `prec` and
 * multiply by `B`.
 */
int nmod_poly_mat_solve_series_dixon(nmod_poly_mat_t X,
                                     const nmod_poly_mat_t A,
                                     const nmod_poly_mat_t B,
                                     slong order,
                                     slong prec);

/** Same as ::nmod_poly_mat_solve_series_dixon, with the use of stored
 * transforms given by `precomp`: 1 to use them (when available, see
 * ::nmod_poly_mat_mulmid_precomp_right_init), 0 not to, -1 for the automatic
 * choice; the default order (nonpositive `order`) follows this choice. For
 * tests and profiles (with `precomp = 0`, this is the version without stored
 * transforms). */
int _nmod_poly_mat_solve_series_dixon(nmod_poly_mat_t X,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t B,
                                      slong order,
                                      slong prec,
                                      int precomp);

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

/** Computation of the denominator in ::nmod_poly_mat_linsolve_via_series */
typedef enum
{
    LINSOLVE_RECONSTRUCTION = 0, /**< by rational reconstruction */
    LINSOLVE_DETERMINANT = 1,    /**< as the determinant, then reduced */
} linsolve_den_t;

/** Computes `u` and `den` such that `u A = den b` (see @ref PolySolve), for a
 * square matrix `A`, from a power series solution. Same specification as
 * ::nmod_poly_mat_linsolve_via_kernel (return value, normalization of the
 * output, aliasing), except that `A` must be square: an exception is thrown
 * otherwise. If `A` is nonsingular, `u / den` is the solution of `y A = b`
 * in lowest terms, with `den` monic.
 *
 * Algorithm, for `A` of size `n x n` with `n >= 1`, and `b` nonzero:
 *  - degree bounds: `deg(det(A)) <= D = min(|rdeg(A)|, |cdeg(A)|)`, and the
 *  entries of `b adj(A)` (the numerators of Cramer's rule) have degree at
 *  most `N = min(deg(b) + |rdeg(A)| - min(rdeg(A)), max_k (deg(b_k) -
 *  cdeg_k(A)) + |cdeg(A)|)`: the `j`-th one is the determinant of `A` with
 *  row `j` replaced by `b`, that is, `sum_k b_k C_{j,k}` with `C_{j,k}` the
 *  cofactor of `A` at `(j,k)`, of degree at most `|rdeg(A)| - rdeg_j(A)` and
 *  at most `|cdeg(A)| - cdeg_k(A)`. The solution in lowest terms satisfies
 *  the same bounds;
 *  - expansion point: if `A(0)` is singular, the system is solved for `A(x +
 *  a)` and `b(x + a)` instead, for some `a` such that `A(a)` is invertible,
 *  and the solution is shifted back. The point `a` is chosen at random among
 *  a few values; if none is suitable (this is the case for all of them if `A`
 *  is singular, and may also happen for small moduli), the solution is
 *  computed by ::nmod_poly_mat_linsolve_via_kernel;
 *  - with `den_method` being ::LINSOLVE_RECONSTRUCTION: the power series
 *  solution `X = b A^{-1} mod x^prec` (see @ref SeriesSolve) is computed by
 *  ::nmod_poly_mat_solve_series_dixon at precision `prec = N' + D + 1`, where
 *  `N' = max(N, D-1)`. The denominator `q` of a random linear combination
 *  `s` of the entries of `X` is found by rational reconstruction with degree
 *  bounds `N'` and `D`, by a half-gcd on the `2 D` coefficients of `s` of
 *  degree `N' - D + 1` to `N' + D` (see the code); `q` divides the
 *  denominator of the solution in lowest terms. If `w = q X mod x^prec` has
 *  degree at most `N'`, then `w / q` is the solution and `q` is a multiple
 *  of its denominator: this certifies that `w / q` is the solution in lowest
 *  terms. Generically, this holds for the first linear combination;
 *  otherwise, this is repeated with another one, and after a few failures,
 *  the denominator is computed as below. The output is always correct; only
 *  the running time depends on the random choices;
 *  - with `den_method` being ::LINSOLVE_DETERMINANT: `X` is computed at
 *  precision `N + 1`, the determinant of `A` by
 *  ::nmod_poly_mat_determinant, and `u = det(A) X mod x^(N+1)` is `b
 *  adj(A)`; then `u` and `det(A)` are divided by their gcd. This avoids the
 *  rational reconstruction, which needs a precision larger by about `D`.
 *
 * The random choices use a local state with a fixed seed, so that the
 * results and the running times are reproducible.
 *
 * Size of the integers (see @ref IntegerSizes): the conditions of
 * ::nmod_poly_mat_solve_series_dixon at precision `N' + D + 1`, of
 * ::nmod_poly_mat_determinant, and of ::nmod_poly_mat_linsolve_via_kernel
 * (in case it is called) are assumed.
 */
int nmod_poly_mat_linsolve_via_series(nmod_poly_mat_t u,
                                      nmod_poly_t den,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t b,
                                      linsolve_den_t den_method);

//@} // doxygen group: Solution over the polynomials

#ifdef __cplusplus
}
#endif

#endif  /* NMOD_POLY_MAT_LINSOLVE_H */
