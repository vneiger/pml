/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/** \file nmod_mat_extra.h */

#ifndef __NMOD_MAT_EXTRA__H
#define __NMOD_MAT_EXTRA__H

#include <flint/flint.h>
#include <flint/perm.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* RANDOM MATRICES                                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Random
 *
 *  These functions are for creating random matrices. They take as input
 *  `mat`, already allocated with some row and column dimensions, and fill them
 *  with random entries with some constraints: prescribed rank, row echelon
 *  form, ...
 */
//@{

/*------------------------------------------------------------*/
/* Uniform Random                                             */
/*------------------------------------------------------------*/

/** Combines Flint's `randrank` and `randops` functions to obtain a random
 * dense matrix `mat` having specified rank */
void nmod_mat_randrank_dense(nmod_mat_t mat,
                             flint_rand_t state,
                             slong rank);

/** Combines Flint's `randrank` and `randops` functions to obtain a random
 * dense matrix having full rank */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank(nmod_mat_t mat, flint_rand_t state)
{
    nmod_mat_randrank_dense(mat, state, FLINT_MIN(mat->r,mat->c));
}

/*------------------------------------------------------------*/
/* Random Lower Row Echelon                                   */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "lower" row echelon form of specified rank. Lower
 * means it has a lower triangular shape, that is, pivots are defined as
 * the indices of rightmost nonzero entries. If `unit` is `1`, pivot entries
 * are ones, otherwise they are random nonzero entries. */
void nmod_mat_rand_lref(nmod_mat_t mat,
                        flint_rand_t state,
                        slong rank,
                        int unit);

/** Sets `mat` to a random full rank "lower" row echelon form. This is
 * the same as calling ::nmod_mat_rand_lref with rank the minimum of row and
 * column dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_lref(nmod_mat_t mat,
                            flint_rand_t state,
                            int unit)
{
    nmod_mat_rand_lref(mat, state, FLINT_MIN(mat->r,mat->c), unit);
}

/*------------------------------------------------------------*/
/* Random Reduced Lower Row Echelon                           */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "lower" reduced row echelon form of specified rank.
 * Lower means it has a lower triangular shape, that is, pivots are defined as
 * the indices of rightmost nonzero entries. */
void nmod_mat_rand_lrref(nmod_mat_t mat,
                         flint_rand_t state,
                         slong rank);

/** Sets `mat` to a random full rank "lower" reduced row echelon form. This is
 * the same as calling ::nmod_mat_rand_lrref with rank the minimum of row and
 * column dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_lrref(nmod_mat_t mat,
                             flint_rand_t state)
{
    nmod_mat_rand_lrref(mat, state, FLINT_MIN(mat->r,mat->c));
}

/*------------------------------------------------------------*/
/* Random Upper Row Echelon                                   */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "upper" row echelon form of specified rank. Upper
 * means it has an upper triangular shape, that is, pivots are defined as the
 * indices of leftmost nonzero entries. If `unit` is `1`, pivots are ones,
 * otherwise they are random nonzero entries. */
void nmod_mat_rand_uref(nmod_mat_t mat,
                        flint_rand_t state,
                        slong rank,
                        int unit);

/** Sets `mat` to a random full rank "upper" row echelon form. This is
 * the same as calling ::nmod_mat_rand_uref with rank the minimum of row and
 * column dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_uref(nmod_mat_t mat,
                            flint_rand_t state,
                            int unit)
{
    nmod_mat_rand_uref(mat, state, FLINT_MIN(mat->r,mat->c), unit);
}

/*------------------------------------------------------------*/
/* Random Reduced Upper Row Echelon                           */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "upper" reduced row echelon form of specified rank.
 * Upper means it has an upper triangular shape, that is, pivots are defined as
 * the indices of leftmost nonzero entries. */
void nmod_mat_rand_urref(nmod_mat_t mat,
                         flint_rand_t state,
                         slong rank);

/** Sets `mat` to a random full rank "upper" reduced row echelon form. This is
 * the same as calling ::nmod_mat_rand_urref with rank the minimum of row and
 * column dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_urref(nmod_mat_t mat,
                             flint_rand_t state)
{
    nmod_mat_rand_urref(mat, state, FLINT_MIN(mat->r,mat->c));
}

/*------------------------------------------------------------*/
/* Random Lower Column Echelon                                */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "lower" column echelon form of specified rank. Lower
 * means it has a lower triangular shape, that is, pivots are defined as the
 * indices of uppermost nonzero entries.  If `unit` is `1`, pivots are ones,
 * otherwise they are random nonzero entries. */
NMOD_MAT_INLINE void
nmod_mat_rand_lcef(nmod_mat_t mat,
                   flint_rand_t state,
                   slong rank,
                   int unit)
{
    nmod_mat_t tmp;
    nmod_mat_init(tmp, mat->c, mat->r, mat->mod.n);
    nmod_mat_rand_uref(tmp, state, rank, unit);
    nmod_mat_transpose(mat, tmp);
    nmod_mat_clear(tmp);
}

/** Sets `mat` to a random full rank "lower" column echelon form. This is the
 * same as calling ::nmod_mat_rand_lcef with rank the minimum of row and column
 * dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_lcef(nmod_mat_t mat,
                            flint_rand_t state,
                            int unit)
{
    nmod_mat_rand_lcef(mat, state, FLINT_MIN(mat->r,mat->c), unit);
}

/*------------------------------------------------------------*/
/* Random Reduced Lower Column Echelon                        */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "lower" reduced column echelon form of specified
 * rank. Lower means it has a lower triangular shape, that is, pivots are
 * defined as the indices of uppermost nonzero entries. */
NMOD_MAT_INLINE void
nmod_mat_rand_lrcef(nmod_mat_t mat,
                    flint_rand_t state,
                    slong rank)
{
    nmod_mat_t tmp;
    nmod_mat_init(tmp, mat->c, mat->r, mat->mod.n);
    nmod_mat_rand_urref(tmp, state, rank);
    nmod_mat_transpose(mat, tmp);
    nmod_mat_clear(tmp);
}

/** Sets `mat` to a random full rank "lower" reduced column echelon form. This
 * is the same as calling ::nmod_mat_rand_lrcef with rank the minimum of row
 * and column dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_lrcef(nmod_mat_t mat,
                             flint_rand_t state)
{
    nmod_mat_rand_lrcef(mat, state, FLINT_MIN(mat->r,mat->c));
}

/*------------------------------------------------------------*/
/* Random Upper Column Echelon                                */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "upper" column echelon form of specified rank. Upper
 * means it has an upper triangular shape, that is, pivots are defined as the
 * indices of bottommost nonzero entries.  If `unit` is `1`, pivots are ones,
 * otherwise they are random nonzero entries. */
NMOD_MAT_INLINE void
nmod_mat_rand_ucef(nmod_mat_t mat,
                   flint_rand_t state,
                   slong rank,
                   int unit)
{
    nmod_mat_t tmp;
    nmod_mat_init(tmp, mat->c, mat->r, mat->mod.n);
    nmod_mat_rand_lref(tmp, state, rank, unit);
    nmod_mat_transpose(mat, tmp);
    nmod_mat_clear(tmp);
}

/** Sets `mat` to a random full rank "upper" column echelon form. This is the
 * same as calling ::nmod_mat_rand_ucef with rank the minimum of row and column
 * dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_ucef(nmod_mat_t mat,
                            flint_rand_t state,
                            int unit)
{
    nmod_mat_rand_ucef(mat, state, FLINT_MIN(mat->r,mat->c), unit);
}

/*------------------------------------------------------------*/
/* Random Reduced Upper Column Echelon                        */
/*------------------------------------------------------------*/

/** Sets `mat` to a random "upper" reduced column echelon form of specified
 * rank. Upper means it has an upper triangular shape, that is, pivots are
 * defined as the indices of bottommost nonzero entries. */
NMOD_MAT_INLINE void
nmod_mat_rand_urcef(nmod_mat_t mat,
                    flint_rand_t state,
                    slong rank)
{
    nmod_mat_t tmp;
    nmod_mat_init(tmp, mat->c, mat->r, mat->mod.n);
    nmod_mat_rand_lrref(tmp, state, rank);
    nmod_mat_transpose(mat, tmp);
    nmod_mat_clear(tmp);
}

/** Sets `mat` to a random full rank "upper" reduced column echelon form. This
 * is the same as calling ::nmod_mat_rand_urcef with rank the minimum of row
 * and column dimensions. */
NMOD_MAT_INLINE void
nmod_mat_rand_fullrank_urcef(nmod_mat_t mat,
                             flint_rand_t state)
{
    nmod_mat_rand_urcef(mat, state, FLINT_MIN(mat->r,mat->c));
}

//@} // doxygen group: Random

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ROW/COLUMN ROTATIONS                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Rotations
 *
 *  These functions are for row rotations and column rotations.
 *  \todo could be in generics? (same functions for nmod_poly_mat, see
 *  the row variants therein: almost zero modification...)
 */
//@{

// rotate rows of mat from i to j (requirement: 0 <= i <= j < mat->r)
// and apply the corresponding transformation to vec (requirement: j < len(vec))
// If i == j, then nothing happens.
// vec can be NULL, in case it is omitted
// More precisely this performs simultaneously:
//      mat[i,:]     <--    mat[j,:]
//      mat[i+1,:]   <--    mat[i,:]
//      mat[i+2,:]   <--    mat[i+1,:]
//         ...       <--       ...
//      mat[j-1,:]   <--    mat[j-2,:]
//      mat[j,:]     <--    mat[j-1,:]
// as well as
//      vec[i]     <--    vec[j]
//      vec[i+1]   <--    vec[i]
//      vec[i+2]   <--    vec[i+1]
//        ...      <--      ...
//      vec[j-1]   <--    vec[j-2]
//      vec[j]     <--    vec[j-1]
void _nmod_mat_rotate_rows_downward(nmod_mat_t mat, slong * vec, slong i, slong j);

// rotate rows of mat from i to j (requirement: 0 <= i <= j < mat->r)
// and apply the corresponding transformation to vec (requirement: j < len(vec))
// If i == j, then nothing happens.
// vec can be NULL, in case it is omitted
// More precisely this performs simultaneously:
//      mat[i,:]     <--    mat[i+1,:]
//      mat[i+1,:]   <--    mat[i+2,:]
//      mat[i+2,:]   <--       ...
//         ...       <--    mat[j-1,:]
//      mat[j-1,:]   <--    mat[j,:]
//      mat[j,:]     <--    mat[i,:]
// as well as
//      vec[i]     <--    vec[i+1]
//      vec[i+1]   <--    vec[i+2]
//      vec[i+2]   <--      ...
//        ...      <--    vec[j-1]
//      vec[j-1]   <--    vec[j]
//      vec[j]     <--    vec[i]
void _nmod_mat_rotate_rows_upward(nmod_mat_t mat, slong * vec, slong i, slong j);


// rotate columns of mat from i to j (requirement: 0 <= i <= j < mat->c)
// and apply the corresponding transformation to vec (requirement: j < len(vec))
// If i == j, then nothing happens.
// vec can be NULL, in case it is omitted
// More precisely this performs simultaneously:
//      mat[:,i]     <--    mat[:,j]
//      mat[:,i+1]   <--    mat[:,i]
//      mat[:,i+2]   <--    mat[:,i+1]
//         ...       <--       ...
//      mat[:,j-1]   <--    mat[:,j-2]
//      mat[:,j]     <--    mat[:,j-1]
// as well as
//      vec[i]     <--    vec[j]
//      vec[i+1]   <--    vec[i]
//      vec[i+2]   <--    vec[i+1]
//        ...      <--      ...
//      vec[j-1]   <--    vec[j-2]
//      vec[j]     <--    vec[j-1]
void _nmod_mat_rotate_columns_rightward(nmod_mat_t mat, slong * vec, slong i, slong j);

// rotate columns of mat from i to j (requirement: 0 <= i <= j < mat->c)
// and apply the corresponding transformation to vec (requirement: j < len(vec))
// If i == j, then nothing happens.
// vec can be NULL, in case it is omitted
// More precisely this performs simultaneously:
//      mat[:,i]     <--    mat[:,i+1]
//      mat[:,i+1]   <--    mat[:,i+2]
//      mat[:,i+2]   <--       ...
//         ...       <--    mat[:,j-1]
//      mat[:,j-1]   <--    mat[:,j]
//      mat[:,j]     <--    mat[:,i]
// as well as
//      vec[i]     <--    vec[i+1]
//      vec[i+1]   <--    vec[i+2]
//      vec[i+2]   <--      ...
//        ...      <--    vec[j-1]
//      vec[j-1]   <--    vec[j]
//      vec[j]     <--    vec[i]
void _nmod_mat_rotate_columns_leftward(nmod_mat_t mat, slong * vec, slong i, slong j);

/** Permute columns of a matrix `mat` according to `perm_act`, and propagate
 * the action on `perm_store`.
 * That is, performs for each appropriate index `j`, the operations
 * `perm_store[j] <- perm_store[perm_act[j]]`
 * `mat[i][j] <- mat[i][perm_act[j]] for all row indices i`
 **/
void nmod_mat_permute_columns(nmod_mat_t mat, const slong * perm_act, slong * perm_store);

//@} // doxygen group: Rotations

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CHARACTERISTIC POLYNOMIAL                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Characteristic polynomial
 *
 * Characteristic polynomial `cp = det(x I - A)` of a square matrix `A`, for
 * a prime modulus, through the determinant of a smaller polynomial matrix
 * (see nmod_poly_mat_determinant.h). The `_lnz` family uses the
 * Labahn-Neiger-Zhou determinant, after a reduction by linear algebra; it
 * may fail on non-generic matrices, and choosing another method then is left
 * to the caller.
 */
//@{

/** Characteristic polynomial `cp = det(x I - A)` of the `n x n` matrix `A`,
 * for a prime modulus, through a block Krylov reduction followed by a
 * polynomial matrix determinant. Let `k = ceil(n / nlevels)` (with `nlevels`
 * first clamped to `[1, n]`), and write `n = q*k + r` with `0 <= r < k`.
 *
 * If `k == n` (in particular if `nlevels <= 1`), there is no reduction:
 * `cp` is the determinant of `x I - A` computed by
 * ::_nmod_poly_mat_det_lnz (first step of the Labahn-Neiger-Zhou
 * triangularization, then the general interface), and the function returns
 * 1.
 *
 * Otherwise, the Krylov matrix `K = [E | A E | ... | A^(q-1) E | A^q E_r]`
 * of the first `k` unit vectors is formed, where `E` is made of the first
 * `k` columns of the identity and `E_r` of its first `r` columns, so that
 * the column `A^j e_i` is at position `j*k + i`. If `K` is invertible, then
 * `K^{-1} A K` is a block companion matrix, and `cp` is the determinant of
 * the `k x k` polynomial matrix `P` with column degrees `q+1` (`r` times) and
 * `q`, whose column `i` is `x^(delta_i) e_i - sum_j x^j z_(j,i)`, where
 * `delta_i` is the column degree and `z_(j,i)` is the block of `k`
 * coefficients of `A^(delta_i) e_i` on the vectors `A^j e_0, ..., A^j
 * e_(k-1)` in the basis `K`. This determinant is computed by
 * ::nmod_poly_mat_determinant, and the function returns 1. For `nlevels =
 * 2` (resp. 4), `P` is about `n/2 x n/2` of degree 2 (resp. `n/4 x n/4` of
 * degree 4).
 *
 * If `K` is singular, the function returns 0 and `cp` is unchanged. Writing
 * `K = [[I_k, X], [0, Y]]`, this happens if and only if `Y` (dimension
 * `n-k`) is singular, which is not the case for a generic `A`, but is the
 * case for example for triangular, diagonal, block diagonal (with blocks
 * larger than `k`), or companion matrices, and for any matrix with more than
 * `k` invariant factors. When `K` has at least two nontrivial blocks of
 * columns (`n - k > k`), the rank of the first one is checked before
 * computing the next ones, which detects most structured cases at little
 * cost.
 *
 * Cost: `n^2 (n - k)` operations for the Krylov matrix (as `q` products by
 * `A` of `n x k` matrices), `(n-k)^3/3 + (n-k)^2 k + (n-k) k^2` for the
 * coefficients `z_(j,i)` by linear system solving, and then the determinant
 * of `P`. On random matrices, the latter dominates for small `nlevels`
 * (e.g. 2 or 4) while it is negligible for `nlevels` around 64, which is
 * the fastest choice in practice: for `n >= 20` and `nlevels = 64`, this is
 * 1.1 to 4 times faster than FLINT's `nmod_mat_charpoly` for a 60-bit prime,
 * and 3 to 11 times faster for a 20-bit prime (Intel Xeon with AVX-512,
 * dimension up to 1024).
 *
 * \todo nlevels threshold to be finely tuned
 */
int _nmod_mat_charpoly_lnz_block_krylov(nmod_poly_t cp, const nmod_mat_t A, slong nlevels);

//@} // doxygen group: Characteristic polynomial

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* KRYLOV MATRICES, BASES AND KERNELS                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Krylov matrices, bases and kernels
 * \anchor KrylovMatrices
 *
 * <b>Krylov matrix.</b> Let `E` be an `m x n` matrix with rows `E_0, ...,
 * E_{m-1}`, `M` an `n x n` matrix, `s = shift` a tuple of `m` integers (the
 * shift) and `d = degrees` a tuple of `m` nonnegative degree bounds. The row
 * `E_i M^j` is designated by the pair `(i, j)`, and the pairs are sorted by
 * increasing priority `(s[i] + j, i)`, lexicographically: a smaller `s[i] +
 * j` comes first, and ties are broken by the row index. The Krylov matrix
 * `K_{s,d}(E, M)` is formed by the rows `(i, j)` for `0 <= j <= d[i]`, in
 * this order. For the zero shift and `d[i] = d` for all `i`, it is `[E ; E M ;
 * ... ; E M^d]`; for `s = (0, k, 2k, ...)` with `k > max(d)`, the iterates of
 * `E_0` come first, then those of `E_1`, and so on. Only the differences
 * between the entries of `s` matter.
 *
 * <b>Krylov basis.</b> If `E_i M^j` is a linear combination of the rows which
 * precede it, then so is `E_i M^{j+1}` (staircase property, [JNSV17, Lemma
 * 7.6]). Hence the row rank profile of the Krylov matrix (its first `r`
 * linearly independent rows, `r` being its rank) consists of the rows `(i,
 * j)` with `j < delta[i]`, for some integers `delta[i]` called the minimal
 * degrees, with `delta[0] + ... + delta[m-1] = r <= n`, provided that the
 * degree bounds are valid: `delta[i] <= d[i]` for all `i`, that is, each `E_i
 * M^{d[i]}` is a linear combination of the rows which precede it. The bounds
 * `d[i] = n` are always valid (Cayley-Hamilton), and so is the degree of the
 * minimal polynomial of `M`. These rows form the Krylov basis.
 *
 * <b>Kernel.</b> The vectors `p = (p_0, ..., p_{m-1})` of polynomials such
 * that `E_0 p_0(M) + ... + E_{m-1} p_{m-1}(M) = 0` form a free module of rank
 * `m` (the coefficients of such a `p` form a vector in the left kernel of a
 * Krylov matrix). Its basis in `s`-Popov form (see @ref MatrixForms;
 * row-wise, with ROW_LOWER pivots) is unique: its pivots are on the diagonal,
 * of degrees `delta[0], ..., delta[m-1]`, and its row `i` is
 * `x^{delta[i]} e_i - sum_k R[i,k] x^{exps[k]} e_{rows[k]}`, where `(rows[k],
 * exps[k])` runs over the Krylov basis and `E_i M^{delta[i]} = sum_k R[i,k]
 * E_{rows[k]} M^{exps[k]}` [JNSV17, Lemma 7.11]. The functions below give
 * this basis in compact form (`R`, `rows`, `exps`, `deltas`);
 * ::nmod_poly_mat_krylov_kernel_basis gives it as a polynomial matrix.
 *
 * <b>Jordan matrices.</b> The functions with suffix `_jordan` take for `M` a
 * Jordan matrix, block diagonal with `nblocks` blocks `x_b I_{t_b} +
 * Z_{t_b}`, where `Z_t` is the `t x t` matrix with ones on its superdiagonal
 * (so that the row vector `v Z_t` is `v` shifted by one position to the
 * right), given by `eig[b] = x_b` (`eig == NULL`: all zero) and `sizes[b] =
 * t_b > 0` (`sizes == NULL`: all one). A row is then multiplied by `M` in
 * `O(n)` operations. Identifying the block `b` of columns of `E_i` with the
 * coefficients of a polynomial `f_{i,b}` of degree less than `t_b`, the row
 * `E_i p(M)` corresponds to the polynomials `p(x + x_b) f_{i,b} mod x^{t_b}`
 * [JNSV17, Section 1.2]: all eigenvalues zero gives approximants
 * (Hermite-Pade approximation), all sizes one gives interpolants, and the
 * general case M-Pade approximants; see
 * ::nmod_poly_mat_approximant_basis_krylov,
 * ::nmod_poly_mat_interpolant_basis_krylov and
 * ::nmod_poly_mat_mpade_basis_krylov.
 *
 * <b>Conventions.</b>
 * - The modulus must be prime.
 * - `shift == NULL` means the zero shift; `degrees == NULL` means `d[i] = n`
 *   for all `i`. Except for the Krylov matrix functions, degree bounds larger
 *   than `n` are replaced by `n`, and the bounds must be valid; the kernel
 *   functions detect when they are not.
 * - Size of the integers (see @ref IntegerSizes): there is no condition on
 *   the shift, whose entries may be arbitrary `slong` values (they are only
 *   compared, after an internal normalization); degree bounds must be
 *   nonnegative (an exception is raised otherwise).
 *
 * Algorithm: the rows are processed in priority order, generated only when
 * needed and stopped at the first dependent row of each index, and they are
 * eliminated by chunks of at most `c` consecutive rows (`c <= 128`) using
 * `nmod_mat_ple` and matrix products. At most `r + m` rows are needed: the
 * `r` rows of the basis and, for the kernels, the first dependent row of each
 * index. When `r = n`, the elimination stops, since all later rows are
 * dependent: the basis alone may then cost about `n` rows (for instance,
 * generically, the first `n` rows of `E` when `m >= n`), and the kernels
 * compute the remaining dependent rows `E_i M^{delta[i]}` without
 * eliminating them. A few more rows may be computed: the rows of an index
 * computed in the same chunk after its first dependent row are discarded, at
 * most `min(c, delta[i]) - 1` for index `i` (fewer than `r` in total). The
 * cost is `O((r + m) n^2)` plus `O(r^3 + m r^2)` for the kernel; see
 * `nmod_mat_extra/krylov_basis.c` for details. The interface follows the
 * methods `krylov_matrix`, `krylov_basis` and `krylov_kernel_basis` of
 * SageMath matrices.
 *
 * Reference: [JNSV17] C.-P. Jeannerod, V. Neiger, E. Schost, G. Villard.
 * Computing minimal interpolation bases. Journal of Symbolic Computation 83
 * (2017), Section 7.
 */
//@{

/** Sets `K` to the Krylov matrix `K_{s,d}(E, M)` (see @ref KrylovMatrices),
 * where `s = shift` and `d = degrees`. The matrix `K` must have `N` rows and
 * `n` columns, where `N` is the sum of the `d[i] + 1` for `0 <= i < m`. The degree bounds need not be
 * valid, and they are not replaced by `n` when larger.
 */
void nmod_mat_krylov_matrix(nmod_mat_t K,
                            const nmod_mat_t E,
                            const nmod_mat_t M,
                            const slong * shift,
                            const slong * degrees);

/** Same as ::nmod_mat_krylov_matrix, for a Jordan matrix `M` given by `eig`,
 * `sizes` and `nblocks` (see @ref KrylovMatrices). */
void nmod_mat_krylov_matrix_jordan(nmod_mat_t K,
                                   const nmod_mat_t E,
                                   nn_srcptr eig,
                                   const slong * sizes,
                                   slong nblocks,
                                   const slong * shift,
                                   const slong * degrees);

/** Computes the Krylov basis of `K_{s,d}(E, M)` (see @ref KrylovMatrices),
 * where `s = shift` and `d = degrees`, and returns its rank `r`. The first
 * `r` rows of `B` are set to the Krylov basis, in priority order, and
 * `(rows[k], exps[k])` to the pair `(i, j)` of its row `k`, for `0 <= k < r`.
 * The matrix `B` must have `n` columns and at least `min(n, N)` rows, where
 * `N` is the sum of the `d[i] + 1` (`n` rows always suffice); `rows` and `exps` must
 * have room for `r` entries (`n` always suffices).
 *
 * The validity of the degree bounds is not checked. If they are not valid,
 * the output is the row rank profile of a matrix which differs from
 * `K_{s,d}(E, M)`, and some pair `(i, d[i])` appears in it (this can be
 * checked in `O(r)`).
 */
slong nmod_mat_krylov_basis(nmod_mat_t B,
                            slong * rows,
                            slong * exps,
                            const nmod_mat_t E,
                            const nmod_mat_t M,
                            const slong * shift,
                            const slong * degrees);

/** Same as ::nmod_mat_krylov_basis, for a Jordan matrix `M` given by `eig`,
 * `sizes` and `nblocks` (see @ref KrylovMatrices). */
slong nmod_mat_krylov_basis_jordan(nmod_mat_t B,
                                   slong * rows,
                                   slong * exps,
                                   const nmod_mat_t E,
                                   nn_srcptr eig,
                                   const slong * sizes,
                                   slong nblocks,
                                   const slong * shift,
                                   const slong * degrees);

/** Computes the kernel of the Krylov matrix `K_{s,d}(E, M)` in compact form
 * (see @ref KrylovMatrices), where `s = shift` and `d = degrees`: on
 * success, returns the rank `r` of the Krylov matrix, sets `(rows[k],
 * exps[k])` to the Krylov basis as in ::nmod_mat_krylov_basis, `deltas[i]` to
 * the minimal degrees for `0 <= i < m`, and the first `r` columns of `R` to
 * the relations
 * `E_i M^{deltas[i]} = sum_k R[i,k] E_{rows[k]} M^{exps[k]}`. The matrix `R`
 * must have `m` rows and at least `r` columns (`n` columns always suffice);
 * `rows` and `exps` must have room for `r` entries and `deltas` for `m`
 * entries.
 *
 * Returns `-1` if the degree bounds are found not to be valid (some row
 * `(i, d[i])` is in the row rank profile); the outputs are then undefined.
 * This detection is exact.
 */
slong nmod_mat_krylov_kernel(nmod_mat_t R,
                             slong * rows,
                             slong * exps,
                             slong * deltas,
                             const nmod_mat_t E,
                             const nmod_mat_t M,
                             const slong * shift,
                             const slong * degrees);

/** Same as ::nmod_mat_krylov_kernel, for a Jordan matrix `M` given by `eig`,
 * `sizes` and `nblocks` (see @ref KrylovMatrices). */
slong nmod_mat_krylov_kernel_jordan(nmod_mat_t R,
                                    slong * rows,
                                    slong * exps,
                                    slong * deltas,
                                    const nmod_mat_t E,
                                    nn_srcptr eig,
                                    const slong * sizes,
                                    slong nblocks,
                                    const slong * shift,
                                    const slong * degrees);

//@} // doxygen group: Krylov matrices, bases and kernels

#ifdef __cplusplus
}
#endif

#endif  // __NMOD_MAT_EXTRA__H
