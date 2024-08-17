/** \file nmod_mat_extra.h
 *
 * \todo
 * - helper: pivots | nonpivots of reduced row echelon form (rightmost/leftmost) ???
 * - PLUQ factorization (try naive + Crout)
 * - left nullspace via PLUQ (or at least more direct than transposition)
 * - row/column rank profile (extract more information from Flint's nullspace?)
 *
 */

#ifndef __NMOD_MAT_EXTRA__H
#define __NMOD_MAT_EXTRA__H

#include <flint/flint.h>
#include <flint/perm.h>
#include <flint/nmod_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MATRIX MULTIPLICATION                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** matrix multiplication using AVX2 instructions for moduli less than 2^30 (TODO refine bound) */
void nmod_mat_mul_small_modulus(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
/** matrix-vector multiplication using AVX2 instructions for moduli less than 2^30 (TODO refine bound) */
// v[:A->r] = A[:,:len]*u[:len], v may not alias u
void nmod_mat_mul_nmod_vec_small_modulus(nn_ptr v, const nmod_mat_t A, nn_srcptr u, ulong len);

/** matrix multiplication not using AVX2 instructions */
// C = A*B
void nmod_mat_mul_newdot(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
// v[:A->r] = A[:,:len]*u[:len], v may not alias u
void nmod_mat_mul_nmod_vec_newdot(nn_ptr v, const nmod_mat_t A, nn_srcptr u, ulong len);

NMOD_MAT_INLINE
void nmod_mat_mul_pml(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    if (A->mod.n < (1L < 29))
        nmod_mat_mul_small_modulus(C, A, B);
    else
        nmod_mat_mul(C, A, B);
}

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

/** Fills matrix with uniformly random entries */
void nmod_mat_rand(nmod_mat_t mat, flint_rand_t state);

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
/* PLUQ DECOMPOSITION                                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name PLUQ
 *
 *  These functions are for PLUQ decomposition and related tasks.
 *  \todo doc
 */
//@{


/** PLUQ decomposition.
 *
 * For an m x n matrix A, computes a PLUQ decomposition
 * --> returns P, LU, Q, rank such that:
 *  - LU is L and U compactly stored ## TODO be more precise
 *  - P is a row permutation (list [0...m-1], permuted)
 *  - Q is a column permutation (list [0...n-1], permuted)
 *  - rank is the rank of A
 *  - A = P*L*U*Q
 *
 * input: P and Q have been allocated as identity permutations of the right
 * length, typically with
 *  slong * P = _perm_init(A->r);
 *  slong * Q = _perm_init(A->c);
 *
 * Rotations are performed to preserve the row rank profile, so that in the end
 * we have P = rrp+nrrp, where rrp is the row rank profile of A and nrrp is the
 * list of row indices not in rrp (both rrp and nrrp are increasing).
 *
 * Using lexicographic order, using row rotations + column rotations
 * ==> preserves the rank profile matrix and row (resp.column) relations precedence
 *    
 * \todo with the current computation of P,Q,
 * the functions nmod_permute_{rows,columns} are such that one should
 * first invert P and Q before using them:
 * _perm_inv(P, P, LU->r);
 * nmod_permute_rows(LU, P, NULL);
 */
slong nmod_mat_pluq(nmod_mat_t A, slong * P, slong * Q);
slong nmod_mat_pluq_crout(nmod_mat_t A, slong * P, slong * Q);

/** Expand LU stored in compact format into two triangular matrices L and U.
 * Conventions, for LU of dimensions m x n of rank rk:
 * - L is m x m unit lower triangular (with the 1's not stored in LU)
 * - only first rk columns of L are stored in LU, the rest are implicit
 *   identity columns
 * - U is m x n upper triangular, with only the first `rk` rows nonzero
 * Requirements:
 * - the rank rk has to be provided
 * - U may alias LU, L may not alias LU
 * - for "zeroin" variant only: on input, L is zero and (either U aliases LU or U is zero)
 **/
void _nmod_mat_l_u_from_compactlu_zeroin(nmod_mat_t L,
                                         nmod_mat_t U,
                                         const nmod_mat_t LU,
                                         slong rank);
void nmod_mat_l_u_from_compactlu(nmod_mat_t L,
                                 nmod_mat_t U,
                                 const nmod_mat_t LU,
                                 slong rank);

//@} // doxygen group: PLUQ

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


/** Left nullspace of A.
 *
 *  Computes a basis X for the left nullspace of A, in reduced row echelon form
 *  with pivots being the rightmost nonzero entries. X should be given
 *  uninitialized, it will be initialized during the call with the right number
 *  of rows.
 *
 *  This first calls nmod_mat_left_nullspace_compact(), and expands the compact
 *  nullspace representation given by this call into the complete dense
 *  nullspace representation.
 *
 * \param [out] X matrix where the nullspace will be stored (uninitialized)
 * \param [in] A input matrix
 * \return nullity of A (i.e. rank of X)
 *
 * @see nmod_mat_left_nullspace_compact
 */
FLINT_DLL slong nmod_mat_left_nullspace(nmod_mat_t X, const nmod_mat_t A);

/** Left nullspace of A in compact form.
 *
 *  Computes a basis X for the left nullspace of A, in reduced row echelon form
 *  with pivots being the rightmost nonzero entries. Only the nonpivot columns
 *  of X are stored, in the order they appear in the nullspace basis. The list
 *  permutation contains the concatenation of two lists, each in increasing
 *  order: the positions of the columns without pivots in the rref nullspace
 *  basis, and the positions of the columns with pivots in the rref nullspace
 *  basis (the former being also the row rank profile of A).
 *
 * \param [out] X matrix where the nullspace will be stored (uninitialized)
 * \param [out] permutation list, allocated with A->r elements
 * \param [in] A input matrix
 * \return nullity of A (i.e. rank of X)
 *
 * \todo efficiency is probably not best for small matrices: this uses Flint's
 * right nullspace and matrix transposition
 */
FLINT_DLL slong nmod_mat_left_nullspace_compact(
                                                nmod_mat_t X,
                                                slong * permutation,
                                                const nmod_mat_t A
                                                );

#ifdef __cplusplus
}
#endif

#endif  // __NMOD_MAT_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
