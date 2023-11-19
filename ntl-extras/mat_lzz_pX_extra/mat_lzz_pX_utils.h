#ifndef MAT_LZZ_PX_UTILS__H
#define MAT_LZZ_PX_UTILS__H

/** \brief Basic routines for univariate polynomial matrices over `zz_p`
 *
 * \file mat_lzz_pX_utils.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost.
 * \version 0.1
 * \date 2018-12-09
 *
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <vector>

/** Vectors of long's, for representing lists of degrees, list of indices,
 * degree shifts, etc.
 */
typedef std::vector<long> VecLong;

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ZERO MATRIX AND IDENTITY MATRIX                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Zero matrix and identity matrix.
 *
 *  Functions to deal with the zero matrix (clear, test if zero)
 *  and the identity matrix (set to identity, test if identity). 
 */
//@{

/** Clears the matrix, that is, sets `pmat` to be zero without modifying its
 * row and column dimensions.
 */
void clear(Mat<zz_pX> & pmat);

/** Clears the submatrix of the matrix `pmat` starting at `(r_offset,c_offset)`
 * and with dimensions `nrows x ncols`. If this involves indices that are out
 * of the bounds defined by the dimensions of `pmat`, then we discard them and
 * restrict to the submatrix indeed contained in `pmat`. The four integer
 * parameters should be nonnegative (this is not checked by the function). */
void clear(Mat<zz_pX> & pmat, long r_offset, long c_offset, long nrows, long ncols);

/** Sets the matrix `pmat` to be the (square) identity matrix with `dim`
 * rows and `dim` columns.
 */
void ident(Mat<zz_pX> & pmat, long dim);

/** Build the identity matrix of the specified dimension `dim`, and return it.
 */
Mat<zz_pX> ident_mat_zz_pX(long dim);

/** Test whether `pvec` is the zero vector (whatever its dimension). */
bool IsZero(const Vec<zz_pX> & pvec);

/** Test whether `pmat` is the zero matrix (whatever its dimensions). */
bool IsZero(const Mat<zz_pX> & pmat);

/** Tests whether `pmat` is a square matrix of some dimension, and is the
 * identity matrix of that dimension.
 */
bool IsIdent(const Mat<zz_pX> & pmat);

/** Tests whether `pmat` is the identity matrix of the specified dimension. */
bool IsIdent(const Mat<zz_pX> & pmat, long dim);

//@} // doxygen group: Zero matrix and identity matrix

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
long deg(const Vec<zz_pX> & pvec);

/** Compute and return the degree of a polynomial matrix `pmat` */
long deg(const Mat<zz_pX> & pmat);

/** Tests whether `pmat` is a constant matrix, that is, of degree 0. */
bool IsConstant(const Mat<zz_pX> & pmat);

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

/** Sets `coeff` to be the coefficient of `pmat` of degree `i` */
void GetCoeff(Mat<zz_p> & coeff, const Mat<zz_pX> & pmat, long i);

/** Builds and returns the coefficient of `pmat` of degree `i` */
inline Mat<zz_p> coeff(const Mat<zz_pX> & pmat, long i)
{ Mat<zz_p> coeff; GetCoeff(coeff, pmat, i); return coeff; }

/** Sets the coefficient of `pmat` of degree `i` to be `coeff` */
void SetCoeff(Mat<zz_pX> & pmat, long i, const Mat<zz_p> & coeff);

//@} // doxygen group: Setting and getting coefficients


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TRANSPOSE, TRUNC, SHIFT, REVERSE                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Transpose                                                  */
/*------------------------------------------------------------*/
/** @name Transpose and mirror
 *
 *  Transpose of a polynomial matrix, and mirror of a vector.
 */
//@{
/** Computes the transpose `tmat` or the polynomial matrix `pmat`; `tmat`
 * may alias `pmat` */
void transpose(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat);

/** Computes and returns the transpose or the polynomial matrix `pmat` */
inline Mat<zz_pX> transpose(const Mat<zz_pX> & pmat)
{ Mat<zz_pX> tmat; transpose(tmat, pmat); return tmat; }

/** Computes the mirror `mvec` of the polynomial vector `pvec`, meaning that
 * `mvec[i]` is `pvec[n-1-i]` for `0<=i<n`, where `n` is the length of `pvec`
 */
void mirror(Vec<zz_pX> & mvec, const Vec<zz_pX> & pvec);

/** Computes and returns the mirror of the polynomial vector `pvec`, as
 * defined above (see #mirror(Vec<zz_pX> & mvec, const Vec<zz_pX> & pvec))
 */
inline Vec<zz_pX> mirror(const Vec<zz_pX> & pvec)
{ Vec<zz_pX> mvec; mirror(mvec, pvec); return mvec; }

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
void trunc(Vec<zz_pX> & tvec, const Vec<zz_pX> & pvec, long n);
/** Computes and returns the truncation of a polynomial vector `pvec` at order `n` */
inline Vec<zz_pX> trunc(const Vec<zz_pX> & pvec, long n)
{ Vec<zz_pX> tvec; trunc(tvec, pvec, n); return tvec; }

/** Computes the truncation `tmat` of a polynomial matrix `pmat` at order `n` */
void trunc(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat, long n);
/** Computes and returns the truncation of a polynomial matrix `pmat` at order `n` */
inline Mat<zz_pX> trunc(const Mat<zz_pX> & pmat, long n)
{ Mat<zz_pX> tmat; trunc(tmat, pmat, n); return tmat; }

/** Computes the matrix `tmat` which is the polynomial matrix `pmat` with its
 * row `i` truncated at order `n` */
void truncRow(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat, long i, long n);
/** Computes and returns (a copy of) the polynomial matrix `pmat` with its row
 * `i` truncated at order `n` */
inline Mat<zz_pX> truncRow(const Mat<zz_pX> & pmat, long i, long n)
{ Mat<zz_pX> tmat; truncRow(tmat, pmat, i, n); return tmat; }

/** Computes the matrix `tmat` which is the polynomial matrix `pmat` with its
 * column `j` truncated at order `n` */
void truncCol(Mat<zz_pX> & tmat, const Mat<zz_pX> & pmat, long j, long n);
/** Computes and returns (a copy of) the polynomial matrix `pmat` with its
 * column `j` truncated at order `n` */
inline Mat<zz_pX> truncCol(const Mat<zz_pX> & pmat, long j, long n)
{ Mat<zz_pX> tmat; truncCol(tmat, pmat, j, n); return tmat; }

//@} // doxygen group: Truncate


/** @name Shifts (multiplication by powers of the variable)
 *
 * Left and right shift operations (X is the variable):
 *  - LeftShift by n means multiplication by X^n
 *  - RightShift by n means division by X^n
 *  - a negative shift reverses the direction of the shift.
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
void LeftShift(Vec<zz_pX> & svec, const Vec<zz_pX> & pvec, long n);
/** Computes and returns the left `n`-shift of the polynomial vector `pvec` */
inline Vec<zz_pX> LeftShift(const Vec<zz_pX> & pvec, long n)
{ Vec<zz_pX> svec; LeftShift(svec, pvec, n); return svec; }
/** Operator equivalent of #LeftShift(const Vec<zz_pX> & pvec, long n) */
inline Vec<zz_pX> operator<<(const Vec<zz_pX> & pvec, long n)
{ Vec<zz_pX> svec; LeftShift(svec, pvec, n); return svec; }
/** Operator: `pvec <<= n;` performs a left `n`-shift on the polynomial vector `pvec` */
inline Vec<zz_pX> & operator<<=(Vec<zz_pX> & pvec, long n)
{ LeftShift(pvec, pvec, n); return pvec; }

/** Computes the right `n`-shift `svec` of the polynomial vector `pvec` */
void RightShift(Vec<zz_pX> & svec, const Vec<zz_pX> & pvec, long n);
/** Computes and returns the right `n`-shift of the polynomial vector `pvec` */
inline Vec<zz_pX> RightShift(const Vec<zz_pX> & pvec, long n)
{ Vec<zz_pX> svec; RightShift(svec, pvec, n); return svec; }
/** Operator equivalent of #RightShift(const Vec<zz_pX> & pvec, long n) */
inline Vec<zz_pX> operator>>(const Vec<zz_pX> & pvec, long n)
{ Vec<zz_pX> svec; RightShift(svec, pvec, n); return svec; }
/** Operator: `pvec >>= n;` performs a right `n`-shift on the polynomial vector `pvec` */
inline Vec<zz_pX> & operator>>=(Vec<zz_pX> & pvec, long n)
{ RightShift(pvec, pvec, n); return pvec; }

/** Computes the left `n`-shift `smat` of the polynomial matrix `pmat` */
void LeftShift(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, long n);
/** Computes and returns the left `n`-shift of the polynomial matrix `pmat` */
inline Mat<zz_pX> LeftShift(const Mat<zz_pX> & pmat, long n)
{ Mat<zz_pX> smat; LeftShift(smat, pmat, n); return smat; }
/** Operator equivalent of #LeftShift(const Mat<zz_pX> & pmat, long n) */
inline Mat<zz_pX> operator<<(const Mat<zz_pX> & pmat, long n)
{ Mat<zz_pX> smat; LeftShift(smat, pmat, n); return smat; }
/** Operator: `pmat <<= n;` performs a left `n`-shift on the polynomial matrix `pmat` */
inline Mat<zz_pX> & operator<<=(Mat<zz_pX> & pmat, long n)
{ LeftShift(pmat, pmat, n); return pmat; }

/** Computes the right `n`-shift `smat` of the polynomial matrix `pmat` */
void RightShift(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, long n);
/** Computes and returns the right `n`-shift of the polynomial matrix `pmat` */
inline Mat<zz_pX> RightShift(const Mat<zz_pX> & pmat, long n)
{ Mat<zz_pX> smat; RightShift(smat, pmat, n); return smat; }
/** Operator equivalent of #RightShift(const Mat<zz_pX> & pmat, long n) */
inline Mat<zz_pX> operator>>(const Mat<zz_pX> & pmat, long n)
{ Mat<zz_pX> smat; RightShift(smat, pmat, n); return smat; }
/** Operator: `pmat >>= n;` performs a right `n`-shift on the polynomial matrix `pmat` */
inline Mat<zz_pX> & operator>>=(Mat<zz_pX> & pmat, long n)
{ RightShift(pmat, pmat, n); return pmat; }


/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its left `n`-shift
 */
void LeftShiftRow(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long i, long n);

/** Computes and returns the matrix which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its right `n`-shift
 */
inline Mat<zz_pX> LeftShiftRow(const Mat<zz_pX> & pmat, const long i, long n)
{ Mat<zz_pX> smat; LeftShiftRow(smat, pmat, i, n); return smat; }

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its right `n`-shift
 */
void RightShiftRow(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long i, long n);

/** Computes and returns the matrix which is the same as the polynomial matrix
 * `pmat` but with its `i`-th row replaced by its right `n`-shift
 */
inline Mat<zz_pX> RightShiftRow(const Mat<zz_pX> & pmat, const long i, long n)
{ Mat<zz_pX> smat; RightShiftRow(smat, pmat, i, n); return smat; }

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its left `n`-shift
 */
void LeftShiftCol(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long j, long n);

/** Computes and returns the matrix which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its left `n`-shift
 */
inline Mat<zz_pX> LeftShiftCol(const Mat<zz_pX> & pmat, const long j, long n)
{ Mat<zz_pX> smat; LeftShiftCol(smat, pmat, j, n); return smat; }

/** Computes the matrix `smat` which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its right `n`-shift
 */
void RightShiftCol(Mat<zz_pX> & smat, const Mat<zz_pX> & pmat, const long j, long n);

/** Computes and returns the matrix which is the same as the polynomial matrix
 * `pmat` but with its `j`-th column replaced by its right `n`-shift
 */
inline Mat<zz_pX> RightShiftCol(const Mat<zz_pX> & pmat, const long j, long n)
{ Mat<zz_pX> smat; RightShiftCol(smat, pmat, j, n); return smat; }

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
void reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat, long hi);

/** Computes and returns the matrix reverse of a polynomial matrix `pmat`, with
 * respect to a bound `hi` (see @ref Reverse for more details).
 */
inline Mat<zz_pX> reverse(const Mat<zz_pX> & pmat, long hi)
{ Mat<zz_pX> rmat; reverse(rmat, pmat, hi); return rmat; }

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to the degree of `pmat` (see @ref Reverse for more details).
 */
inline void reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat)
{ reverse(rmat, pmat, deg(pmat)); }

/** Computes and returns the matrix reverse of a polynomial matrix `pmat`, with
 * respect to the degree of `pmat` (see @ref Reverse for more details).
 */
inline Mat<zz_pX> reverse(const Mat<zz_pX> & pmat)
{ Mat<zz_pX> rmat; reverse(rmat, pmat, deg(pmat)); return rmat; }

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each row (see @ref Reverse for more
 * details). The length of `hi` must be the number of rows of `pmat`.
 */
void row_reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat, const VecLong & hi);

/** Computes and returns the matrix reverse of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each row (see @ref Reverse for more
 * details). The length of `hi` must be the number of rows of `pmat`.
 
 */
inline Mat<zz_pX> row_reverse(const Mat<zz_pX> & pmat, const VecLong & hi)
{ Mat<zz_pX> rmat; row_reverse(rmat, pmat, hi); return rmat; }

/** Computes the matrix reverse `rmat` of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each column (see @ref Reverse for more
 * details). The length of `hi` must be the number of columns of `pmat`.
 */
void col_reverse(Mat<zz_pX> & rmat, const Mat<zz_pX> & pmat, const VecLong & hi);

/** Computes and returns the matrix reverse of a polynomial matrix `pmat`, with
 * respect to a list of bounds `hi` for each column (see @ref Reverse for more
 * details). The length of `hi` must be the number of columns of `pmat`.
 */
inline Mat<zz_pX> col_reverse(const Mat<zz_pX> & pmat, const VecLong & hi)
{ Mat<zz_pX> rmat; col_reverse(rmat, pmat, hi); return rmat; }



//@} // doxygen group: Reverse


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* EVALUATE (SINGLE/MULTIPOINT)                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Evaluation at one or multiple points
 *
 *  Functions for computing the evaluation of a polynomial matrix at a single
 *  point (element of the base field, that is, of type `zz_p`) or at multiple
 *  points.
 *
 * \todo add multipoint evaluation (general points, geometric points, ..)
 */
//@{

/** Computes the evaluation `evmat` of the polynomial matrix `pmat` at the
 * point `pt` (`pmat` stored as a matrix of polynomials) */
void eval(Mat<zz_p> & evmat, const Mat<zz_pX> & pmat, const zz_p & pt);

/** Computes and returns the evaluation of the polynomial matrix `pmat` at the
 * point `pt` (`pmat` stored as a matrix of polynomials) */
inline Mat<zz_p> eval(const Mat<zz_pX> & pmat, const zz_p & pt)
{ Mat<zz_p> evmat; eval(evmat, pmat, pt); return evmat; }

/** Computes the evaluation `evmat` of the polynomial matrix `matp` at the
 * point `pt` (polynomial matrix stored as a vector of constant matrices). If
 * `matp` has length zero, `eval` is the `0x0` matrix. */
void eval(Mat<zz_p> & evmat, const Vec<Mat<zz_p>> & matp, const zz_p & pt);

/** Computes and returns the evaluation of the polynomial matrix `matp` at the
 * point `pt` (polyomial matrix stored as a vector of constant matrices). If
 * `matp` has length zero, this returns the `0x0` matrix.  */
inline Mat<zz_p> eval(const Vec<Mat<zz_p>> & matp, const zz_p & pt)
{ Mat<zz_p> evmat; eval(evmat, matp, pt); return evmat; }

//@} // doxygen group: Evaluation at one or multiple points


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
void random(Vec<zz_pX> & pvec, long n, long d);

/** Computes and returns a random polynomial vector of length `n` and degree
 * less than `d`
 */
inline Vec<zz_pX> random_vec_zz_pX(long n, long d)
{ Vec<zz_pX> pvec; random(pvec, n, d); return pvec; }

/** Computes a random polynomial matrix `pmat` with `m` rows, `n` columns, and
 * degree less than `d`
 */
void random(Mat<zz_pX> & pmat, long m, long n, long d);

/** Computes and returns a random polynomial matrix with `m` rows, `n` columns,
 * and degree less than `d`
 */
inline Mat<zz_pX> random_mat_zz_pX(long m, long n, long d)
{ Mat<zz_pX> pmat; random(pmat, m, n, d); return pmat; }

/** Computes a random polynomial matrix `pmat` with `m` rows, `n` columns, and
 * degree of `i`th row less than `rdeg[i]` for all `i`
 */
void random_mat_zz_pX_rdeg(Mat<zz_pX> & pmat, long m, long n, VecLong rdeg);

/** Computes a random polynomial matrix `pmat` with `m` rows, `n` columns, and
 * degree of `j`th column less than `cdeg[j]` for all `j`
 */
void random_mat_zz_pX_cdeg(Mat<zz_pX> & pmat, long m, long n, VecLong cdeg);

/** Computes a random polynomial matrix `pmat` with the same dimensions
 * as `dmat`, and degree of `(i,j)` entry less than `dmat[i][j]` for
 * all `i` and `j`
 */
void random(Mat<zz_pX> & pmat, Mat<long> dmat);

/** Computes and returns a random polynomial matrix `pmat` with the same
 * dimensions as `dmat`, and degree of `(i,j)` entry less than `dmat[i][j]` for
 * all `i` and `j`
 */
inline Mat<zz_pX> random_mat_zz_pX(Mat<long> dmat)
{ Mat<zz_pX> pmat; random(pmat, dmat); return pmat; }


/** Computes a random polynomial matrix `matp` with `m` rows, `n` columns,
 * and degree less than `d`, represented as a vector of constant matrices.
 */
void random(Vec<Mat<zz_p>> & matp, long m, long n, long d);

/** Computes and returns a random polynomial matrix with `m` rows, `n` columns,
 * and degree less than `d`, represented as a vector of constant matrices.
 */
inline Vec<Mat<zz_p>> random_vec_mat_zz_p(long m, long n, long d)
{ Vec<Mat<zz_p>> matp; random(matp, m, n, d); return matp; }

//@} // doxygen group: Generation of random matrices



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CONVERT                                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Conversion from constant matrix.
 */
//@{

/** Conversion from a constant matrix `mat` into a polynomial matrix `pmat`
 * (which is constant, equal to `mat`).
 */
void conv(Mat<zz_pX> & pmat, const Mat<zz_p> & mat);

/** Forms a polynomial matrix equal to the constant matrix `mat`, and returns
 * it.
 */
inline Mat<zz_pX> conv(const Mat<zz_p> & mat)
{ Mat<zz_pX> pmat; conv(pmat, mat); return pmat; }

//@} // doxygen group: Conversion from constant matrix



/** @name Conversion to/from Vec<Mat<zz_p>>
 *
 *  Two main representations are used for polynomial matrices:
 *     - a polynomial matrix stored as a matrix with polynomial entries, that
 *     is, of type `Mat<zz_pX>` (currently the default representation)
 *     - a polynomial matrix stored as a polynomial with matrix coefficients,
 *     that is, of type `Vec<Mat<zz_p>>`
 *
 *  In the latter representation, the vector has length at least `deg(mat)+1`;
 *  in particular, the zero matrix might be represented by the vector of length
 *  0, although this does not give access to the dimensions of this zero
 *  matrix. Note also the following requirement: in the second representation,
 *  all the matrix coefficients have the same row and column dimensions; this
 *  is currently assumed to hold, and not checked by the algorithms.
 *
 *  The following functions perform conversions between these types, with two
 *  variants: either the degree is deduced from the input, or a user-provided
 *  truncation order is used.
 */
//@{

/** Converts from matrix with polynomial entries `pmat` to polynomial with
 * matrix coefficients `matp`; the zero matrix is converted to the length-0
 * vector, losing information on the row and column dimensions.
 */
void conv(Vec<Mat<zz_p>> & matp, const Mat<zz_pX> & pmat);

/** Returns the conversion of a matrix with polynomial entries `pmat` into a
 * polynomial with matrix coefficients; the zero matrix is converted to the
 * length-0 vector, losing information on the row and column dimensions.
 */
inline Vec<Mat<zz_p>> conv(const Mat<zz_pX> & pmat)
{ Vec<Mat<zz_p>> matp; conv(matp, pmat); return matp; }

/** Converts from polynomial with matrix coefficients `matp` to matrix with
 * polynomial entries `pmat`; if `matp` has length 0 then `pmat` is cleared
 * (set to zero without changing its dimensions).
 */
void conv(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp);

/** Returns the conversion of a polynomial with matrix coefficients `matp` into
 * a matrix with polynomial entries; if `matp` has length 0 then it returns a
 * `Mat<zz_pX>` of dimensions 0 x 0.
 */
inline Mat<zz_pX> conv(const Vec<Mat<zz_p>> & matp)
{ Mat<zz_pX> pmat; conv(pmat, matp); return pmat; }


/** Converts from matrix with polynomial entries `pmat`, truncated at the
 * specified `order`, to polynomial with matrix coefficients `matp`. The
 * integer `order` must be nonnegative (this is currently not verified by the
 * code), and it is guaranteed that the vector `matp` has length `order`.
 */
void conv(Vec<Mat<zz_p>> & matp, const Mat<zz_pX> & pmat, const long order);

/** Returns the conversion of a matrix with polynomial entries `pmat`,
 * truncated at the specified `order`, to a polynomial with matrix
 * coefficients. The integer `order` must be nonnegative (this is currently not
 * verified by the code), and the output vector is guaranteed to have length
 * `order`.
 */
inline Vec<Mat<zz_p>> conv(const Mat<zz_pX> & pmat, const long order)
{ Vec<Mat<zz_p>> matp; conv(matp, pmat, order); return matp; }

/** Converts from polynomial with matrix coefficients `matp`, truncated at the
 * specified `order` (a nonnegative integer), to matrix with polynomial entries
 * `pmat`; if `matp` has length 0 then `pmat` is cleared (set to zero without
 * changing its dimensions).
 */
void conv(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp, const long order);

/** Returns the conversion from polynomial with matrix coefficients `matp`,
 * truncated at the specified `order` (a nonnegative integer), to matrix with
 * polynomial entries; if `matp` has length 0 then it returns a `Mat<zz_pX>` of
 * dimensions 0 x 0.
 */
inline Mat<zz_pX> conv(const Vec<Mat<zz_p>> & matp, const long order)
{ Mat<zz_pX> pmat; conv(pmat, matp, order); return pmat; }

//@} // doxygen group: Conversion to/from Vec<Mat<zz_p>>

#endif // MAT_LZZ_PX_UTILS__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
