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
#include <flint/flint-config.h>
#include <flint/perm.h>
#include <flint/nmod_mat.h>

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

#ifdef __cplusplus
}
#endif

#endif  // __NMOD_MAT_EXTRA__H
