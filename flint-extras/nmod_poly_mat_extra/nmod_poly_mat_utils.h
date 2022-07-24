#ifndef NMOD_POLY_MAT_UTILS_H
#define NMOD_POLY_MAT_UTILS_H

/** \brief Basic routines for univariate polynomial matrices over `nmod`
 *
 * \file nmod_poly_mat_utils.h
 * \version 0.0
 * \date 2022-06-25
 *
 */

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_mat_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MATRIX DEGREE                                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Degree
 *
 *  The degree of a polynomial matrix is the maximum of the degrees of all its
 *  entries; by convention, the zero matrix has degree -1 (thus, except for the
 *  zero matrix, the degree of a polynomial matrix is always nonnegative).
 */
//@{

/** Compute and return the degree of a polynomial vector `pvec` */
// TODO

/** Compute and return the degree of a polynomial matrix `pmat` */
NMOD_POLY_MAT_INLINE slong
nmod_poly_mat_degree(const nmod_poly_mat_t pmat)
{
    return nmod_poly_mat_max_length(pmat)-1;
}

/** Tests whether `pmat` is a constant matrix, that is, of degree 0 */
NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_constant(const nmod_poly_mat_t pmat)
{
    return (nmod_poly_mat_max_length(pmat) == 1);
}

/** Compute and return the degree of a vector polynomial `vecp` */
// TODO

/** Compute and return the degree of a matrix polynomial `matp` */
NMOD_POLY_MAT_INLINE slong
nmod_mat_poly_degree(const nmod_mat_poly_t matp)
{
    // TODO not guaranteed to be the actual degree?
	return matp->degree;
}


/** Tests whether `matp` is a constant matrix, that is, of degree 0 */
NMOD_POLY_MAT_INLINE int
nmod_mat_poly_is_constant(const nmod_mat_poly_t matp)
{
    // TODO not guaranteed to be the actual degree?
    return matp->degree == 0;
}

//@} // doxygen group: Degree

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SETTING AND GETTING COEFFICIENTS                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Setting and getting coefficients
 *
 *  Seeing a polynomial matrix in `Mat<zz_pX>` as a univariate polynomial with
 *  matrix coefficients in `Mat<zz_p>`, these functions allow one to get or set
 *  one of these coefficients, for a given polynomial matrix.
 */
//@{

/** Sets `coeff` to be the coefficient of `pmat` of degree `degree` */
void coefficient_matrix(nmod_mat_t coeff,
                        const nmod_poly_mat_t pmat,
                        slong degree);

/** Sets the coefficient of `pmat` of degree `degree` to be `coeff` */
// TODO
//void SetCoeff(Mat<zz_pX> & pmat, long i, const Mat<zz_p> & coeff);

//@} // doxygen group: Setting and getting coefficients


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TRANSPOSE, TRUNC, SHIFT, REVERSE                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Transpose and mirror
 *
 *  Transpose of a polynomial matrix, and mirror of a vector.
 */
//@{
/** Computes the transpose `tmat` or the polynomial matrix `pmat`; `tmat`
 * may alias `pmat` */
// TODO
//void transpose(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat);

/** Computes the mirror `mvec` of the polynomial vector `pvec`, meaning that
 * `mvec[i]` is `pvec[n-1-i]` for `0<=i<n`, where `n` is the length of `pvec`
 */
// TODO
//void mirror(Vec<zz_pX> & mvec, const Vec<zz_pX> & pvec);

//@} // doxygen group: Transpose and mirror


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

/** Computes the truncation `tvec` of a polynomial vector `pvec` at order `n` */
// TODO
//void trunc(Vec<zz_pX> & tvec, const Vec<zz_pX> & pvec, long n);

/** Computes the truncation `tmat` of a polynomial matrix `pmat` at order `n` */
// TODO
//void trunc(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat, long n);

/** Computes the matrix `tmat` which is the polynomial matrix `pmat` with its
 * row `i` truncated at order `n` */
// TODO
//void truncRow(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat, long i, long n);

/** Computes the matrix `tmat` which is the polynomial matrix `pmat` with its
 * column `j` truncated at order `n` */
// TODO
//void truncCol(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat, long j, long n);

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
//void LeftShift(Vec<zz_pX> & svec, const Vec<zz_pX> & pvec, long n);

/** Computes the right `n`-shift `svec` of the polynomial vector `pvec` */
// TODO
//void RightShift(Vec<zz_pX> & svec, const Vec<zz_pX> & pvec, long n);

/** Computes the left `n`-shift `smat` of the polynomial matrix `pmat` */
// FIXME investigate what happens with Flint's shift when polynomial is zero
void nmod_poly_mat_shift_left(nmod_poly_mat_t res,
                              const nmod_poly_mat_t pmat, slong k);

/** Computes the right `n`-shift `smat` of the polynomial matrix `pmat` */
// FIXME investigate what happens with Flint's shift when polynomial is zero
void nmod_poly_mat_shift_right(nmod_poly_mat_t res,
                               const nmod_poly_mat_t pmat, slong k);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its left `n`-shift
 */
// TODO
//void LeftShiftRow(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long i, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its right `n`-shift
 */
// TODO
//void RightShiftRow(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long i, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its left `n`-shift
 */
// TODO
//void LeftShiftCol(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long j, long n);

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its right `n`-shift
 */
// TODO
//void RightShiftCol(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long j, long n);

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
//void reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat, long hi);

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to the degree of `pmat` (see @ref Reverse for more details).
 */
//TODO
//inline void reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat)
//{ reverse(rmat, pmat, deg(pmat)); }

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each row (see @ref Reverse for more
 * details). The length of `hi` must be the number of rows of `pmat`.
 */
//TODO
//void row_reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat, const VecLong & hi);

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each column (see @ref Reverse for more
 * details). The length of `hi` must be the number of columns of `pmat`.
 */
//TODO
//void col_reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat, const VecLong & hi);

//@} // doxygen group: Reverse


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* EVALUATE (SINGLE POINT)                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Evaluation at one point
 *
 *  Functions for computing the evaluation of a polynomial matrix at a single
 *  point (element of the base field, that is, of type `nmod`).
 *
 */
//@{

// /** Computes the evaluation `evmat` of the polynomial matrix `pmat` at the
//  * point `pt` (`pmat` stored as a matrix of polynomials) */
// Flint native, as nmod_poly_mat_evaluate_nmod

/** Computes the evaluation `evmat` of the polynomial matrix `matp` at the
 * point `pt` (polynomial matrix stored as a vector of constant matrices). If
 * `matp` has length zero, `eval` is the `0x0` matrix. */
// TODO
//void eval(Mat<zz_p> & evmat, const Vec<Mat<zz_p>> & matp, const zz_p & pt);

//@} // doxygen group: Evaluation at one points


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Generation of random matrices
 *
 *  Functions for generating a random polynomial matrix of given
 *  dimensions and degree; the degree bound can either be a global
 *  degree bound, or a list of bounds for the degree of each row
 *  (respectively, of each column).
 *
 *  \todo complete this doc with ref to function for random matrix with given form
 */
//@{

/** Computes a random polynomial vector `pvec` of length `n` and degree less
 * than `d`
 */
// TODO
//void random(Vec<zz_pX> & pvec, long n, long d);

/** Computes a random polynomial matrix `pmat` with `m` rows, `n` columns, and
 * degree less than `d`
 */
// TODO
//void random(Mat<zz_pX> & pmat, long m, long n, long d);

/** Computes a random polynomial matrix `pmat` with `m` rows, `n` columns, and
 * degree of `i`th row less than `rdeg[i]` for all `i`
 */
// TODO
//void random_mat_zz_pX_rdeg(Mat<zz_pX> & pmat, long m, long n, VecLong rdeg);

/** Computes a random polynomial matrix `pmat` with `m` rows, `n` columns, and
 * degree of `j`th column less than `cdeg[j]` for all `j`
 */
// TODO
//void random_mat_zz_pX_cdeg(Mat<zz_pX> & pmat, long m, long n, VecLong cdeg);

/** Computes a random polynomial matrix `pmat` with the same dimensions
 * as `dmat`, and degree of `(i,j)` entry less than `dmat[i][j]` for
 * all `i` and `j`
 */
// TODO
//void random(Mat<zz_pX> & pmat, Mat<long> dmat);

/** Computes a random polynomial matrix `matp` with `m` rows, `n` columns,
 * and degree less than `d`, represented as a vector of constant matrices.
 */
// TODO
//void random(Vec<Mat<zz_p>> & matp, long m, long n, long d);

//@} // doxygen group: Generation of random matrices



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CONVERT                                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Conversion from constant matrix.
 */
//@{

/** Initialize a polynomial matrix `pmat` of the same dimensions and modulus as
 * a constant matrix `mat`, and set its constant coefficient to `mat`.
 * \todo careful: probably better to use nmod_poly_init with preinv...
 * but this should have been done in Flint's native poly_mat_init too?
 **/
//void nmod_poly_mat_init_set_from_nmod_mat(nmod_poly_mat_t pmat, const nmod_mat_t cmat);

/** Set the polynomial matrix `pmat` to be a constant polynomial matrix whose
 * constant coefficient is a copy of `cmat`. This assume `pmat` is already
 * initialized with the same modulus and dimensions of `cmat`.
 **/
void nmod_poly_mat_set_from_nmod_mat(nmod_poly_mat_t pmat,
                                     const nmod_mat_t cmat);

/** Conversion from a constant matrix `mat` into a matrix polynomial `matp`
 * (which is constant, equal to `mat`).
 * \todo
 */


//@} // doxygen group: Conversion from constant matrix


/** @name Conversion nmod_poly_mat <-> nmod_mat_poly
 *
 *  Two main representations are used for polynomial matrices:
 *     - a polynomial matrix stored as a matrix with polynomial entries, that
 *     is, of type `nmod_poly_mat_t`
 *     - a polynomial matrix stored as a polynomial with matrix coefficients,
 *     that is, of type `nmod_mat_poly_t`
 *
 *  In the latter representation, the array of matrices has length at least
 *  `degree(mat)+1`; in particular, the zero matrix may be represented by an
 *  array of length 0. Note also the following requirement: in the second
 *  representation, all the matrix coefficients have the same row and column
 *  dimensions; this is currently assumed to hold, and is not checked by the
 *  algorithms.
 *
 *  The following functions perform conversions between these types, with two
 *  variants: either the degree is deduced from the input, or a user-provided
 *  truncation order is used.
 */
//@{

/** Converts from matrix with polynomial entries `pmat` to polynomial with
 * matrix coefficients `matp`.
 */
//void nmod_poly_mat_to_mat_poly(nmod_mat_poly_t matp,
			    //const nmod_poly_mat_t pmat);

/** Converts from polynomial with matrix coefficients `matp` to matrix with
 * polynomial entries `pmat`.
 * \todo improvements in implementation
 **/
void nmod_poly_mat_set_from_mat_poly(nmod_poly_mat_t pmat,
                                     const nmod_mat_poly_t matp);

/** Converts from matrix with polynomial entries `pmat`, truncated at the
 * specified `order`, to polynomial with matrix coefficients `matp`. The
 * integer `order` must be nonnegative (this is currently not verified by the
 * code), and it is guaranteed that the vector `matp` has length `order`.
 */
// TODO
//void conv(Vec<Mat<zz_p>> & matp, const Mat<zz_pX> & pmat, const long order);

/** Converts from polynomial with matrix coefficients `matp`, truncated at the
 * specified `order` (a nonnegative integer), to matrix with polynomial entries
 * `pmat`; if `matp` has length 0 then `pmat` is cleared (set to zero without
 * changing its dimensions).
 */
// TODO
//void conv(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp, const long order);

//@} // doxygen group: Conversion nmod_poly_mat <-> nmod_mat_poly

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_UTILS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
