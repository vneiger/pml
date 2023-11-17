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

/** matrix multiplication using AVX2 instructions for moduli less than 2^30 */
void nmod_mat_mul_small_modulus(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);



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
