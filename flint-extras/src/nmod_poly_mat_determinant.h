/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_DETERMINANT_H
#define NMOD_POLY_MAT_DETERMINANT_H

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file nmod_poly_mat_determinant.h
 * Determinant algorithms for square polynomial matrices over `nmod`, for a
 * prime modulus.
 *
 * In all functions, `det` is the output determinant and `pmat` is the input
 * square matrix; `det` cannot alias an entry of `pmat`.
 */

/** Mulders-Storjohann triangularization of the square matrix `mat`, in place,
 * by unimodular row operations (iterative triangularization via weak Popov
 * forms of leading submatrices). If `mat` is nonsingular, this returns the
 * determinant (+1 or -1) of the transformation and swaps into `diag[i]` the
 * `i`-th diagonal entry of the obtained upper triangular matrix (`diag` must
 * have `mat->r` initialized entries, which are not entries of `mat`; the
 * entries of `mat` are meaningless afterwards); otherwise, it returns 0 and
 * `diag` is not modified. A `0 x 0` matrix is nonsingular (return value 1). */
slong _nmod_poly_mat_upper_triangularize_iter(nmod_poly_struct * diag, nmod_poly_mat_t mat);

/** Mulders-Storjohann determinant algorithm: the determinant is the product
 * of the diagonal entries given by ::_nmod_poly_mat_upper_triangularize_iter,
 * times the determinant of the transformation. The input matrix is modified
 * (see ::_nmod_poly_mat_upper_triangularize_iter) and must be square (this is
 * not tested). */
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat);

/** Determinant by expansion by minors. For dimension at most 4, this uses
 * explicit formulas which share 2 x 2 minors (as in the NTL-based part of
 * PML). For larger dimensions, this expands along successive rows and
 * stores all the `k x k` minors of the first `k` rows, which costs
 * `n * 2^(n-1)` polynomial multiplications and stores up to `binomial(n,
 * n/2)` polynomials: this is meant for small dimensions only (an exception
 * is raised beyond dimension 16). */
void nmod_poly_mat_det_minors(nmod_poly_t det, const nmod_poly_mat_t pmat);

/** Determinant by evaluation at the points `1, r^2, r^4, ...` of a geometric
 * progression, computation of the determinants of the evaluated constant
 * matrices, and interpolation. The number of points is `D + 1`, where `D =
 * min(|rdeg(pmat)|, |cdeg(pmat)|)` is an upper bound on `deg(det(pmat))`.
 * This requires a prime modulus `p` with `p >= 10 * (D+1)` (see
 * #NMOD_POLY_CAN_USE_GEOMETRIC): if this is not satisfied, the function
 * returns 0 and `det` is left unchanged; otherwise, the function returns 1.
 * Exception: if `pmat` has a zero row or a zero column, `det` is set to zero
 * and 1 is returned, whatever the modulus. The memory used is about `n^2
 * (D+1)` words (all evaluations are stored).
 *
 * \todo ::nmod_poly_mat_degree_det_bound gives a sharper bound `D`, but then
 * the entries may have degree larger than `D` and the evaluation must
 * handle inputs longer than the number of points
 */
int nmod_poly_mat_det_geometric(nmod_poly_t det, const nmod_poly_mat_t pmat);

/** Determinant via the triangularization of Labahn, Neiger, Zhou (Journal
 * of Complexity 2017, Algorithm 2), in a row-wise version (see
 * ::nmod_poly_mat_hnf_diagonal_lnz). Write pmat = [F_L | F_R] with F_L having
 * ceil(n/2) columns, s = rdeg(pmat), K an s-minimal left kernel basis of
 * F_R, and F_R = C * R with R a row basis of F_R; then det(pmat) = kappa *
 * det(K * F_L) * det(R), where kappa is a nonzero constant obtained from
 * K(0) and C(0) only [LNZ17, Lemma 4.4, transposed], and the two
 * determinants are computed recursively. As for the Hermite diagonal, the row
 * basis R is not computed when it is certified to be unimodular (then R = I
 * and C = F_R), which is the generic case; otherwise, R and C are computed by
 * the Zhou-Labahn approach reusing K (see ::_nmod_poly_mat_row_basis_zl). This
 * is deterministic and works for any prime modulus.
 *
 * Size of the integers (see @ref IntegerSizes): as for
 * ::nmod_poly_mat_hnf_diagonal_lnz, the conditions only involve degrees and
 * dimensions, and hold in practice.
 *
 * \todo as for the Hermite diagonal, inputs with unbalanced degrees (e.g. one
 * row of much higher degree) are slow, since the approximation order of the
 * kernel basis follows the maximum row degree
 *
 * \todo for degree-1 inputs, the kernel step could be done by linear
 * algebra: for F_R = [A2 ; x I - A1] with A2 invertible, K = [x I - S | -A2]
 * with S = A2 A1 A2^{-1} is a minimal kernel basis and K F_L is monic
 * (n = 800, 60-bit: first step 1000 ms -> 163 ms); in general, kernel rows
 * correspond to linear relations among the rows of A2, A2 A1, A2 A1^2, ...
 * (observability indices of (A1, A2)), and the same holds at deeper levels
 * for [N ; D] with D monic, via the block companion of D (an all-linear-
 * algebra LNZ for characteristic matrices); see also the Hermite diagonal
 * of characteristic matrices
 */
void nmod_poly_mat_det_lnz(nmod_poly_t det, const nmod_poly_mat_t pmat);

/** Same as ::nmod_poly_mat_det_lnz if `hybrid` is zero. Otherwise, the
 * determinants of the two matrices obtained at the first step of the
 * triangularization are computed by the general interface
 * ::nmod_poly_mat_determinant (which may itself use this function). */
void _nmod_poly_mat_det_lnz(nmod_poly_t det, const nmod_poly_mat_t pmat, int hybrid);

/** General interface for the determinant: chooses among the above
 * algorithms depending on the dimension `n`, the degree `d`, and the
 * modulus:
 * - expansion by minors for `n <= 5`;
 * - determinant of a constant matrix if `d == 0`;
 * - zero if `pmat` has a zero row or column;
 * - for `d == 1`, writing `pmat = x A + B`: characteristic polynomial of
 *   `-A^{-1} B` if `A` is invertible, or reversed characteristic polynomial
 *   of `-B^{-1} A` if `B` is invertible (for a characteristic matrix `x I -
 *   M`, calling a characteristic polynomial function directly avoids the
 *   factorization of `A = I`);
 * - otherwise, the transpose is used if it has smaller sum of row degrees,
 *   and then, with `D = min(|rdeg(pmat)|, |cdeg(pmat)|)`: evaluation /
 *   interpolation at a geometric progression for `n <= 24` and large `d`
 *   (`d >= 100`, or `d >= 30` for `n <= 12`, or `d >= 20` for `n <= 8`) if
 *   `n*d <= 3*D`, the modulus allows it and the memory needed is at most 256
 *   MiB; Mulders-Storjohann's iterative algorithm if `n*d > 3*D` (strongly
 *   unbalanced degrees, e.g. a few rows of much larger degree), or `n <= 16`
 *   and `n*d <= 320`, or `d <= 1` and `n <= 32`; and otherwise, one step of
 *   the Labahn-Neiger-Zhou triangularization with the two resulting
 *   determinants computed by this same interface.
 *
 * \todo thresholds to be retuned on a few machines and moduli: they come
 * from timings on an Intel Xeon with AVX-512 with a 60-bit prime, before the
 * speed-up of the weak Popov iteration; the ratio 3 between `n*d` and `D` is
 * the crossover observed between Mulders-Storjohann and LNZ for `n <= 256`
 *
 * \todo a well-made degree-1 path: the characteristic polynomial of `-A^{-1}
 * B` is computed by FLINT's `nmod_mat_charpoly` (Danilevsky), whereas
 * ::_nmod_mat_charpoly_lnz_block_krylov (nmod_mat_extra.h) is 1.1 to 11
 * times faster on generic matrices of dimension 20 to 1024; this needs a
 * charpoly interface with a fallback for the non-generic matrices on which
 * it fails, and `A^{-1} B` could be applied to the Krylov blocks instead of
 * being formed
 */
void nmod_poly_mat_determinant(nmod_poly_t det, const nmod_poly_mat_t pmat);


#ifdef __cplusplus
}
#endif

#endif /* NMOD_POLY_MAT_DETERMINANT_H */
