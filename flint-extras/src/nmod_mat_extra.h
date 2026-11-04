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

#ifdef __cplusplus
}
#endif

#endif  // __NMOD_MAT_EXTRA__H
