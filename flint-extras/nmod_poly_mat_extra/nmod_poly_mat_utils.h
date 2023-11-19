#ifndef NMOD_POLY_MAT_UTILS_H
#define NMOD_POLY_MAT_UTILS_H

/** \brief Basic routines for univariate polynomial matrices over `nmod`
 *
 * \file nmod_poly_mat_utils.h
 * \version 0.0
 * \date 2023-01-25
 *
 */

#include <flint/nmod_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_mat_poly.h" // TODO remove
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_forms.h"

#ifdef __cplusplus
extern "C" {
#endif

/** The input `window` is some window (r1,c1,r2,c2) of some polynomial matrix
 * `mat`. This function changes it into the window (r1,c1+cc1,r2,c2+cc2) of the
 * same `mat`. It is assumed that the latter window is indeed a valid window of
 * the original matrix `mat`; in particular one must provide cc1,cc2 such that
 * c1+cc1 < mat->c and c2+cc2 < mat->c.  **/
// TODO currently unused, remove?
NMOD_POLY_MAT_INLINE void
_nmod_poly_mat_window_update_columns(nmod_poly_mat_t window, slong cc1, slong cc2)
{
    // if original mat->c was 0, each window->rows[i] is and should remain NULL
    // otherwise, window->rows[i] are all non-NULL and should be updated
    if (window->rows)  // NULL if window->r <= 0
        for (slong i = 0; i < window->r; i++)
            if (window->rows[i]) 
                window->rows[i] += cc1;
    window->c = cc2-cc1;
}

/** The input `window` is some window (r1,c1,r2,c2) of some polynomial matrix
 * `mat`. This function changes it into the window (r1,c1,r2,c2+c) of the same
 * `mat`. The provided c can be positive or negative and must be such that the
 * latter window is a valid window of `mat`: 0 <= c2+c <= mat->c.  **/
NMOD_POLY_MAT_INLINE void
_nmod_poly_mat_window_resize_columns(nmod_poly_mat_t window, slong c)
{
    window->c += c;
}


/** Tests whether `pmat` is a unimodular matrix, that is, a square matrix with
 * determinant being a nonzero constant. Rectangular matrices are accepted as
 * input (and are not unimodular). */
// deterministic test, slow; and randomized test via evaluation (not reliable for small fields)
int nmod_poly_mat_is_unimodular(const nmod_poly_mat_t pmat);
int nmod_poly_mat_is_unimodular_randomized(const nmod_poly_mat_t pmat, flint_rand_t state);

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
void _nmod_poly_mat_rotate_rows_downward(nmod_poly_mat_t mat, slong * vec, slong i, slong j);

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
void _nmod_poly_mat_rotate_rows_upward(nmod_poly_mat_t mat, slong * vec, slong i, slong j);

// . mat is the input matrix
// . r is an integer between 0 and mat->r -1
// . vec has >= r entries (signed integers)
// . perm must be allocated with at least mat->r entries (no need to fill it
// with values)
// This finds the unique permutation which stable-sorts vec[0:r] into
// non-decreasing order, records this permutation in perm[0:r], and applies it
// both to vec and to the first r rows of mat (the out-vec[i] is the
// in-vec[perm[i]], the out-mat[i,:] is the in-mat[perm[i],:]). The entries
// perm[r:mat->r] are filled with iota (perm[i] = i), to ensure that perm
// actually represents the used row permutation on the whole set of rows.
//
// Example uses:
// - weak Popov -> ordered weak Popov: vec is the s-pivot index of mat in
// s-weak Popov form with r nonzero rows; then this puts mat in s-ordered weak
// Popov form and ensures that vec is updated accordingly; also, the output
// perm can be used to apply the corresponding row permutation to another
// matrix with as many rows as mat (e.g. some unimodular transformation)
// - step in weak Popov -> Popov: vec is the s-pivot degree of mat in s-weak
// Popov form with r nonzero rows; then this arranges the s-pivot degree in
// nondecreasing order, and the relative position of rows with identical
// s-pivot degree is preserved.
void _nmod_poly_mat_permute_rows_by_sorting_vec(nmod_poly_mat_t mat,
                                                slong r,
                                                slong * vec,
                                                slong * perm);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SWAP, PERMUTE, TRANSPOSE                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Swap, permute, transpose
 *
 *  These functions allow one to swap or permute the rows or the columns a
 *  given polynomial matrix, and to get its transpose.
 */
//@{

/** Swap two rows of a polynomial matrix. It is assumed that `r` and `s` are
 * valid row indices for `mat` (this is not checked). The case of equality
 * `r==s` is allowed. This swaps pointers to rows. If `perm` is not `NULL`,
 * it should be an array for which `r` and `s` are valid indices; then the
 * corresponding elements of this array will be swapped. */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_swap_rows(nmod_poly_mat_t mat,
                        slong * perm,
                        slong r, slong s)
{
    if (r != s)
    {
        if (perm)
        {
            slong t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        nmod_poly_struct * tmp = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = tmp;
    }
}

/** Invert rows of a polynomial matrix. If `mat` has `m` rows, then its row `0`
 * is swapped with its row `m-1`, its row `1` is swapped with its row `m-2`,
 * etc. If `perm` is not `NULL`, it should be an array of length `m`; its
 * elements of this array will be inverted accordingly. */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_invert_rows(nmod_poly_mat_t mat, slong * perm)
{
    for (slong i = 0; i < mat->r/2; i++)
        nmod_poly_mat_swap_rows(mat, perm, i, mat->r - i - 1);
}

/** Swap two columns of a polynomial matrix. It is assumed that `r` and `s` are
 * valid column indices for `mat` (this is not checked). The case of equality
 * `r==s` is allowed. This swaps pointers to polynomial entries. If `perm` is
 * not `NULL`, it should be an array for which `r` and `s` are valid indices;
 * then the corresponding elements of this array will be swapped. */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_swap_columns(nmod_poly_mat_t mat,
                           slong * perm,
                           slong r, slong s)
{
    if (r != s)
    {
        if (perm)
        {
            slong t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        for (slong t = 0; t < mat->r; t++)
            nmod_poly_swap(mat->rows[t] + r, mat->rows[t] + s);
    }
}

/** Invert columns of a polynomial matrix. If `mat` has `n` columns, then its
 * column `0` is swapped with its column `n-1`, its column `1` is swapped with
 * its column `n-2`, etc. If `perm` is not `NULL`, it should be an array of
 * length `n`; its elements of this array will be inverted accordingly. */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_invert_columns(nmod_poly_mat_t mat, slong * perm)
{
    if (perm)
    {
        for (slong j = 0; j < mat->c/2; j++)
        {
            slong t = perm[j];
            perm[j] = perm[mat->c - j - 1];
            perm[mat->c - j - 1] = t;
        }
    }

    for (slong t = 0; t < mat->r; t++)
        for (slong j = 0; j < mat->c/2; j++)
            nmod_poly_swap(mat->rows[t]+j, mat->rows[t]+ mat->c - j - 1);
}


/** Permute rows of a polynomial matrix `mat` according to `perm_act`, and
 * propagate the action on `perm_store`. Namely, performs for each appropriate
 * index `i`, the operations `perm_store[i] <- perm_store[perm_act[i]]` and
 * `rows[i] <- rows[perm_act[i]]`.  */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_permute_rows(nmod_poly_mat_t mat,
                           const slong *perm_act,
                           slong * perm_store)
{
    // perm_store[i] <- perm_store[perm_act[i]] 
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    // rows[i] <- rows[perm_act[i]] 
    nmod_poly_struct ** mat_tmp = flint_malloc(mat->r * sizeof(nmod_poly_struct *));
    for (slong i = 0; i < mat->r; ++i)
        mat_tmp[i] = mat->rows[perm_act[i]];
    for (slong i = 0; i < mat->r; ++i)
        mat->rows[i] = mat_tmp[i];

    flint_free(mat_tmp);
}

/** Permute columns of a polynomial matrix `mat` according to `perm_act`, and
 * propagate the action on `perm_store`. Namely, performs for each appropriate
 * index `i` and `j`, the operations `perm_store[j] <- perm_store[perm_act[j]]`
 * and `rows[i][j] <- rows[i][perm_act[j]]`.  */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_permute_columns(nmod_poly_mat_t mat,
                              const slong * perm_act,
                              slong * perm_store)
{
    // perm_store[i] <- perm_store[perm_act[i]] 
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->c);


    // simulate columns[j] <- columns[perm_act[j]]
    nmod_poly_struct * row_tmp = flint_malloc(mat->c * sizeof(nmod_poly_struct));
    for (slong i = 0; i < mat->r; i++)
    {
        for (slong j = 0; j < mat->c; j++)
            row_tmp[j] = mat->rows[i][perm_act[j]];
        for (slong j = 0; j < mat->c; j++)
            mat->rows[i][j] = row_tmp[j];
    }
    flint_free(row_tmp);
}

/** Computes the transpose `tmat` or the polynomial matrix `pmat`; `tmat`
 * may alias `pmat` */
// TODO
//void transpose(nmod_poly_mat tmat, nmod_poly_mat pmat);

/** Computes the mirror `mvec` of the polynomial vector `pvec`, meaning that
 * `mvec[i]` is `pvec[n-1-i]` for `0<=i<n`, where `n` is the length of `pvec`
 */
// TODO
//void mirror( mvec,  pvec);

//@} // doxygen group: Swap, permute, transpose



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TRANSPOSE, TRUNC, SHIFT, REVERSE                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/** @name Truncate
 *
 *  Truncate a polynomial vector, a polynomial matrix, or a specific row or
 *  column of a polynomial matrix.
 *
 *  In all the following functions which involve an OUT parameter (`tvec` or
 *  `tmat`), this parameter may alias the IN parameter (`pvec` or `pmat`).
 *
 * \todo different truncation orders on columns/rows
 */
//@{

/** Computes the matrix `tmat` which is the polynomial matrix `pmat` with its
 * row `i` truncated at order `len` */
// TODO
//void nmod_poly_mat_truncate_row(nmod_poly_mat_t tmat, const nmod_poly_mat_t pmat, long len, long i);

/** Computes the matrix `tmat` which is the polynomial matrix `pmat` with its
 * column `j` truncated at order `len` */
// TODO
//void nmod_poly_mat_truncate_column(nmod_poly_mat_t tmat, const nmod_poly_mat_t pmat, long len, long j);

//@} // doxygen group: Truncate


/** @name Shifts (multiplication by powers of the variable)
 *
 * Left and right shift operations (X is the variable):
 *  - LeftShift by n means multiplication by X^n
 *  - RightShift by n means division by X^n
 *  - shift `n` has to be nonnegative in both cases.
 *
 *  In all the following functions which involve an OUT parameter (`svec` or
 *  `smat`), this parameter may alias the corresponding IN parameter (`pvec` or
 *  `pmat`).
 *
 * \todo 
 *   - versions with different shifting orders on rows/columns; 
 *   - shiftAdd, shiftSub
 */
//@{

/** Computes the left `n`-shift `svec` of the polynomial vector `pvec` */
// TODO
//void LeftShift( svec,  pvec, long n);

/** Computes the right `n`-shift `svec` of the polynomial vector `pvec` */
// TODO
//void RightShift( svec,  pvec, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its left `n`-shift
 */
// TODO
//void LeftShiftRow(nmod_poly_mat smat, nmod_poly_mat pmat, const long i, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its right `n`-shift
 */
// TODO
//void RightShiftRow(nmod_poly_mat smat, nmod_poly_mat pmat, const long i, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its left `n`-shift
 */
// TODO
//void LeftShiftCol(nmod_poly_mat_t smat, nmod_poly_mat_t pmat, const long j, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its right `n`-shift
 */
// TODO
//void RightShiftCol(nmod_poly_mat_t smat, nmod_poly_mat_t pmat, const long j, long n);

//@} // doxygen group: Shifts (multiplication by powers of the variable)


/** @name Reverse
 * \anchor Reverse
 *
 *  The reverse, sometimes called mirror, of a polynomial is the operation of
 *  reversing the order of its coefficients. This is usually either done with
 *  respect to the degree of that polynomial, or more generally with respect to
 *  some specified bound. The reverse of a polynomial matrix is defined
 *  analogously, with either the degree of the matrix or a specified bound. We
 *  also consider row-wise and column-wise variants, in which such a bound is
 *  specified for each row, or for each column.
 *
 *  Explicitly, for a polynomial `f` of degree `d` with coefficients `f[i]`,
 *  and for a given nonnegative integer bound `hi`, the reverse of `f` with
 *  respect to `hi` is the polynomial of degree at most `hi` whose coefficients
 *  are `f[hi]`, `f[hi-1]`, ..., `f[0]`. This definition can be extended for
 *  negative `hi` by considering that this reverse polynomial is zero.
 *
 *  In the functions below, when there is an OUT parameter `rmat`, it may alias
 *  the IN parameter `pmat`.
 *
 */
//@{

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a bound `hi` (see @ref Reverse for more details).
 */
//TODO
//void reverse(nmod_poly_mat_t rmat, nmod_poly_mat_t pmat, long hi);

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to the degree of `pmat` (see @ref Reverse for more details).
 */
//TODO
//inline void reverse(nmod_poly_mat_t rmat, nmod_poly_mat_t pmat)
//{ reverse(rmat, pmat, deg(pmat)); }

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each row (see @ref Reverse for more
 * details). The length of `hi` must be the number of rows of `pmat`.
 */
//TODO
//void row_reverse(nmod_poly_mat_t rmat, nmod_poly_mat_t  pmat, long * hi);

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each column (see @ref Reverse for more
 * details). The length of `hi` must be the number of columns of `pmat`.
 */
//TODO
//void col_reverse(nmod_poly_mat_t rmat, nmod_poly_mat_t pmat, long * hi);

//@} // doxygen group: Reverse

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Generation of random matrices with degree bounds
 *
 *  Functions for generating a random polynomial matrix of given
 *  dimensions and degree; the degree bound can either be a global
 *  degree bound, or a list of bounds for the degree of each row
 *  (respectively, of each column), or a matrix of bounds.
 *
 *  \todo complete this doc with ref to function for random matrix with given form
 */
//@{

/** Computes a random polynomial vector `pvec` of length `n` and degree less
 * than `d`
 */
// TODO vector version
//void nmod_poly_vec_rand(nmod_poly_vec_t pvec, flint_rand_t state, slong len);

/** Fills the polynomial matrix `pmat` with dense polynomials of length `len`
 * with coefficients taken uniformly at random.  */
void nmod_poly_mat_rand(nmod_poly_mat_t mat,
                        flint_rand_t state,
                        slong len);

/** Fills polynomial matrix `mat` with random dense polynomial entries such
 * that `mat[i,j] has length up to `rdeg[i]+1` for all `i`. Assumes `rdeg` has
 * the right length, i.e. the number of rows of mat. */
void nmod_poly_mat_rand_row_degree(nmod_poly_mat_t mat,
                                   flint_rand_t state,
                                   const slong * rdeg);

/** Fills polynomial matrix `mat` with random dense polynomial entries such
 * that `mat[i,j] has length up to `cdeg[j]+1` for all `j`. Assumes `cdeg` has
 * the right length, i.e. the number of columns of mat. */
void nmod_poly_mat_rand_column_degree(nmod_poly_mat_t mat,
                                      flint_rand_t state,
                                      const slong * cdeg);

/** Fills polynomial matrix `mat` with random dense polynomial entries such
 * that `mat[i,j] has length up to `dmat[i,j]+1` for all `i,j`. Assumes `dmat`
 * has the right number of rows and columns (i.e. the same as those of mat). */
void nmod_poly_mat_rand_degree_matrix(nmod_poly_mat_t mat,
                                      flint_rand_t state,
                                      const fmpz_mat_t dmat);

//@} // doxygen group: Generation of random matrices with degree bounds

/** @name Generation of random matrices with specific forms
 *
 *  Functions for generating a random polynomial matrix of a given
 *  form: reduced, ordered weak Popov, Popov, Hermite; their shifted
 *  variants and with choice of row/column orientation.
 *
 *  For definitions of pivot profiles, pivot indices and pivot degrees,
 *  see @ref Pivots and @ref MatrixForms.
 */
//@{

/** Fills polynomial matrix `mat` with random dense polynomial entries such
 * that `mat` is in `shift`-Popov form row-wise, with `shift`-pivot profile
 * specified by `pivind` and `pivdeg`. If `m` and `n` are the number of rows
 * and columns of `mat`, the input must satisfy `m <= n`, the length of
 * `pivind` and `pivdeg` must be `m`, `pivind` must consist of increasing
 * integers between `0` and `n-1`, `pivdeg` must be nonnegative, `shift` must
 * have length `n`. */
void _nmod_poly_mat_rand_popov_row_lower(nmod_poly_mat_t mat,
                                         flint_rand_t state,
                                         const slong * pivind,
                                         const slong * pivdeg,
                                         const slong * shift);

/** Fills polynomial matrix `mat` with random dense polynomial entries such
 * that `mat` is in `shift`-Popov form column-wise, with `shift`-pivot profile
 * specified by `pivind` and `pivdeg`. If `m` and `n` are the number of rows
 * and columns of `mat`, the input must satisfy `n <= m`, the length of
 * `pivind` and `pivdeg` must be `n`, `pivind` must consist of increasing
 * integers between `0` and `m-1`, `pivdeg` must be nonnegative, `shift` must
 * have length `m`.  */
void _nmod_poly_mat_rand_popov_col_upper(nmod_poly_mat_t mat,
                                         flint_rand_t state,
                                         const slong * pivind,
                                         const slong * pivdeg,
                                         const slong * shift);

/** Fills polynomial matrix `mat` with random dense polynomial entries such
 * that `mat` is in `shift`-Popov form, with `shift`-pivot profile specified by
 * `pivind` and `pivdeg`; the orientation is specified by the input `orient`.
 * The input requirements depend on the orientation, see the documentation of
 * ::nmod_poly_mat_rand_popov_row_lower or ::nmod_poly_mat_rand_popov_col_upper
 * . Here `shift` can be `NULL` (equivalent to `[0,...,0]`), and if `mat` is
 * square, `pivind` can be `NULL` (note that if it is not, the input
 * requirement means that it has to be the list of successive integers
 * [0,1,2,..], of the right length).
 * \todo ROW_UPPER / COL_LOWER not implemented
 **/
void nmod_poly_mat_rand_popov(nmod_poly_mat_t mat,
                              flint_rand_t state,
                              const slong * pivind,
                              const slong * pivdeg,
                              const slong * shift,
                              orientation_t orient);

//@} // doxygen group: Generation of random matrices with specific forms


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM CONSTANT OR MATP                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Set from constant matrix and matrix polynomial.
 *
 * Provides `set` and `init_set`, from either a constant matrix,
 * or from a `nmod_mat_poly_t`.
 *
 * Two main representations are used for polynomial matrices:
 *    - a polynomial matrix stored as a matrix with polynomial entries, that
 *    is, of type `nmod_poly_mat_t`
 *    - a polynomial matrix stored as a polynomial with matrix coefficients,
 *    that is, of type `nmod_mat_poly_t`
 *
 * In the latter representation, the array of matrices has length
 * `degree(mat)+1`; in particular, the zero matrix may be represented by an
 * array of length 0. Note also the following requirement: in the second
 * representation, all the matrix coefficients have the same row and column
 * dimensions; this is currently assumed to hold, and is not checked by the
 * algorithms.
 *
 * The following functions perform conversions between these types, with two
 * variants: either the degree is deduced from the input, or a user-provided
 * truncation order is used.
 *
 * \todo init-set variants
 */
//@{

/** Initialize a polynomial matrix `pmat` of the same dimensions and modulus as
 * a constant matrix `mat`, and set its constant coefficient to `mat`.
 * \todo to be written
 * careful: probably better to use nmod_poly_init with preinv...
 * but this should have been done in Flint's native poly_mat_init too?
 **/
//void nmod_poly_mat_init_set_from_nmod_mat(nmod_poly_mat_t pmat, const nmod_mat_t cmat);

// TODO remove
void nmod_poly_mat_set_from_mat_poly0(nmod_poly_mat_t pmat,
                                     const nmod_mat_poly0_t matp);

/** Set from polynomial with matrix coefficients `matp`, truncated at the
 * specified `order` (a nonnegative integer). */
// TODO benchmark and try variants if needed
FLINT_DLL void
nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                      const nmod_mat_poly_t matp,
                                      slong order);

/** Set from polynomial with matrix coefficients `matp`. */
// TODO benchmark and try variants if needed
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_set_from_mat_poly(nmod_poly_mat_t pmat,
                                const nmod_mat_poly_t matp)
{
    nmod_poly_mat_set_trunc_from_mat_poly(pmat, matp, matp->length);
}

//@} // doxygen group: Conversion from constant matrix and matrix polynomial

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_UTILS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
