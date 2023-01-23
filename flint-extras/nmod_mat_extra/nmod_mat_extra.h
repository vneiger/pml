/** \file nmod_mat_extra.h
 *
 * \todo
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

/** RANDOM */

/** Fills matrix with uniformly random entries */
void nmod_mat_rand(nmod_mat_t mat, flint_rand_t state);














/** Permute rows of a matrix `mat` according to `perm_act`, and propagate the
 * action on `perm_store`. Namely, performs for each appropriate index `i`, the
 * operations `perm_store[i] <- perm_store[perm_act[i]]` and
 * `rows[i] <- rows[perm_act[i]]`.
 * \todo Should be present in a future release of Flint (?) --> then remove */
NMOD_MAT_INLINE void
nmod_mat_permute_rows(nmod_mat_t mat,
                      slong * perm_store,
                      const slong * perm_act)
{
    // perm_store[i] <- perm_store[perm_act[i]] 
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    // rows[i] <- rows[perm_act[i]] 
    mp_limb_t ** mat_tmp = flint_malloc(mat->r * sizeof(mp_limb_t *));
    for (slong i = 0; i < mat->r; ++i)
        mat_tmp[i] = mat->rows[perm_act[i]];
    for (slong i = 0; i < mat->r; ++i)
        mat->rows[i] = mat_tmp[i];

    flint_free(mat_tmp);
}


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
 * \param [in] A input matrix
 * \param [out] X matrix where the nullspace will be stored (uninitialized)
 * \return nullity of A (i.e. rank of X)
 *
 * @see nmod_mat_left_nullspace_compact
 */
FLINT_DLL slong nmod_mat_left_nullspace(nmod_mat_t X, nmod_mat_t A);

/** Left nullspace of A in compact form.
 *
 *  Computes a basis X for the left nullspace of A, in reduced row echelon form
 *  with pivots being the rightmost nonzero entries. Only the nonpivot columns
 *  of X are stored, in the order they appear in the nullspace basis. The list
 *  permutation contains the concatenation of two lists, each in increasing
 *  order: the positions of the columns with pivots in the nullspace, and the
 *  positions of the columns without pivots in the nullspace (the latter being
 *  also the row rank profile of A).
 *
 * \param [in] A input matrix
 * \param [out] X matrix where the nullspace will be stored (uninitialized)
 * \param [out] permutation list, allocated with A->r elements
 * \return nullity of A (i.e. rank of X)
 *
 * \todo better implementation based on PLUQ decomposition; currently
 * this relies on Flint's nullspace and matrix transposition
 */
FLINT_DLL slong nmod_mat_left_nullspace_compact(
                                                nmod_mat_t X,
                                                slong * permutation,
                                                nmod_mat_t A
                                                );


/** Left lower triangular solving: X = B * L^{-1} */
// TODO not yet implemented
FLINT_DLL void nmod_mat_solve_left_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_left_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_left_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);

/** Left upper triangular solving: X = B * U^{-1} */
// TODO not yet implemented
FLINT_DLL void nmod_mat_solve_left_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_left_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);
FLINT_DLL void nmod_mat_solve_left_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit);

// TODO:
// - helper: pivots | nonpivots of reduced row echelon form (rightmost/leftmost) ???
// - PLUQ factorization (try naive + Crout)
// - left nullspace via PLUQ
// - row/column rank profile

#ifdef __cplusplus
}
#endif

#endif  // __NMOD_MAT_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
