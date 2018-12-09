#ifndef MAT_LZZ_PX_UTILS__H
#define MAT_LZZ_PX_UTILS__H

/** Some basic routines for handling `Mat<zz_pX>`.
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

/** Vectors of long's, for representing lists of degrees, list of indices, degree shifts, etc. */
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

/** Sets the matrix `pmat` to be the (square) identity matrix with `dim`
 * rows and `dim` columns.
 */
void ident(Mat<zz_pX> & pmat, long dim);

/** Build the identity matrix of the specified dimension `dim`, and return it.
 */
Mat<zz_pX> ident_mat_zz_pX(long dim);

/** Test whether `pvec` is the zero vector (whatever its dimension).
 * \return 1 if `pvec` is zero, 0 otherwise
 */
long IsZero(const Vec<zz_pX> & pvec);

/** Test whether `pmat` is the zero matrix (whatever its dimensions).
 * \return 1 if `pmat` is zero, 0 otherwise
 */
long IsZero(const Mat<zz_pX> & pmat);

/** Tests whether `pmat` is a square matrix of some dimension, and is the
 * identity matrix of that dimension.
 * \return 1 if `pmat` is the identity matrix, 0 otherwise.
 */
long IsIdent(const Mat<zz_pX> & pmat);

/** Tests whether `pmat` is the identity matrix of the specified dimension.
 * \return 1 if `pmat` is the identity matrix of dimension `dim`, 0 otherwise.
 */
long IsIdent(const Mat<zz_pX> & pmat, long dim);

//@} // doxygen group: Zero matrix and identity matrix

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MATRIX DEGREE                                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Matrix degree
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

//@} // doxygen group: Matrix degree


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
void transpose(Mat<zz_pX> & x, const Mat<zz_pX> & a);

inline Mat<zz_pX> transpose(const Mat<zz_pX> & a)
{ Mat<zz_pX> x; transpose(x, a); return x; }

/*------------------------------------------------------------*/
/* Truncate mod X^..., for all the matrix / some columns/rows */
/* output can alias input                                     */
/*------------------------------------------------------------*/

// vector versions
void trunc(Vec<zz_pX> & x, const Vec<zz_pX> & a, long n);
inline Vec<zz_pX> trunc(const Vec<zz_pX> & a, long n)
{ Vec<zz_pX> x; trunc(x, a, n); return x; }

// full matrix versions
void trunc(Mat<zz_pX> & x, const Mat<zz_pX> & a, long n);
inline Mat<zz_pX> trunc(const Mat<zz_pX> & a, long n)
{ Mat<zz_pX> x; trunc(x, a, n); return x; }

// row versions
void truncRow(Mat<zz_pX> & x, const Mat<zz_pX> & a, long r, long n);
inline Mat<zz_pX> truncRow(const Mat<zz_pX> & a, long r, long n)
{ Mat<zz_pX> x; truncRow(x, a, r, n); return x; }

// col versions
void truncCol(Mat<zz_pX> & x, const Mat<zz_pX> & a, long c, long n);
inline Mat<zz_pX> truncCol(const Mat<zz_pX> & a, long c, long n)
{ Mat<zz_pX> x; truncCol(x, a, c, n); return x; }

/* TODO: different truncation orders on columns/rows          */

/*------------------------------------------------------------*/
/* Shift operations:                                          */
/*  - LeftShift by n means multiplication by X^n              */
/*  - RightShift by n means division by X^n                   */
/*  - a negative shift reverses the direction of the shift.   */
/*------------------------------------------------------------*/

/* TODO                                                       */
/* versions with different shifting orders on rows/columns    */
/* shiftAdd, shiftSub                                         */

// left shift, vector
void LeftShift(Vec<zz_pX> & x, const Vec<zz_pX> & a, long n);
inline Vec<zz_pX> LeftShift(const Vec<zz_pX> & a, long n)
{ Vec<zz_pX> x; LeftShift(x, a, n); return x; }

// right shift, vector
void RightShift(Vec<zz_pX> & x, const Vec<zz_pX> & a, long n);
inline Vec<zz_pX> RightShift(const Vec<zz_pX> & a, long n)
{ Vec<zz_pX> x; RightShift(x, a, n); return x; }


// left shift, full matrix 
void LeftShift(Mat<zz_pX> & x, const Mat<zz_pX> & a, long n);
inline Mat<zz_pX> LeftShift(const Mat<zz_pX> & a, long n)
{ Mat<zz_pX> x; LeftShift(x, a, n); return x; }

// right shift, full matrix 
void RightShift(Mat<zz_pX> & x, const Mat<zz_pX> & a, long n);
inline Mat<zz_pX> RightShift(const Mat<zz_pX> & a, long n)
{ Mat<zz_pX> x; RightShift(x, a, n); return x; }

// left shift, single row
void LeftShiftRow(Mat<zz_pX> & x, const Mat<zz_pX> & a, const long r, long n);
inline Mat<zz_pX> LeftShiftRow(const Mat<zz_pX> & a, const long r, long n)
{ Mat<zz_pX> x; LeftShiftRow(x, a, r, n); return x; }

// right shifts, single row
void RightShiftRow(Mat<zz_pX> & x, const Mat<zz_pX> & a, const long r, long n);
inline Mat<zz_pX> RightShiftRow(const Mat<zz_pX> & a, const long r, long n)
{ Mat<zz_pX> x; RightShiftRow(x, a, r, n); return x; }

// left shift, single column
void LeftShiftCol(Mat<zz_pX> & x, const Mat<zz_pX> & a, const long c, long n);
inline Mat<zz_pX> LeftShiftCol(const Mat<zz_pX> & a, const long c, long n)
{ Mat<zz_pX> x; LeftShiftCol(x, a, c, n); return x; }

// right shifts, single column
void RightShiftCol(Mat<zz_pX> & x, const Mat<zz_pX> & a, const long c, long n);
inline Mat<zz_pX> RightShiftCol(const Mat<zz_pX> & a, const long c, long n)
{ Mat<zz_pX> x; RightShiftCol(x, a, c, n); return x; }

// operators for full matrix shift
inline Mat<zz_pX> operator<<(const Mat<zz_pX> & a, long n)
{ Mat<zz_pX> x; LeftShift(x, a, n); return x; }
inline Mat<zz_pX> operator>>(const Mat<zz_pX> & a, long n)
{ Mat<zz_pX> x; RightShift(x, a, n); return x; }

inline Mat<zz_pX> & operator<<=(Mat<zz_pX> & x, long n)
{ LeftShift(x, x, n); return x; }
inline Mat<zz_pX> & operator>>=(Mat<zz_pX> & x, long n)
{ RightShift(x, x, n); return x; }


/*------------------------------------------------------------*/
/* reverse the order of the entries in a vector               */
/* x = a[n - 1 -i], i=0..n-1, with n=length(a)                */
/*------------------------------------------------------------*/
void reverse_vector(Vec<zz_pX> & x, const Vec<zz_pX> & a);
inline Vec<zz_pX> reverse_vector(const Vec<zz_pX> & a)
{ Vec<zz_pX> x; reverse_vector(x, a); return x; }

/*------------------------------------------------------------*/
/* Reverse operations:                                        */
/* x = reverse of a[0]..a[hi] (hi >= -1);                     */
/* user provided 'hi'                                         */
/*------------------------------------------------------------*/
void reverse(Mat<zz_pX> & x, const Mat<zz_pX> & a, long hi);

inline Mat<zz_pX> reverse(const Mat<zz_pX> & a, long hi)
{ Mat<zz_pX> x; reverse(x, a, hi); return x; }

/* TODO versions with different degree on different cols/rows */
void reverse(
             Mat<zz_pX> & x, 
             const Mat<zz_pX> & a, 
             const VecLong & hi,
             const bool row_wise = true
            );

/*------------------------------------------------------------*/
/* Reverse operations:                                        */
/* x = reverse of a[0]..a[deg(a)]                             */
/*------------------------------------------------------------*/

inline void reverse(Mat<zz_pX> & x, const Mat<zz_pX> & a)
{ reverse(x, a, deg(a)); }

inline Mat<zz_pX> reverse(const Mat<zz_pX> & a)
{ Mat<zz_pX> x; reverse(x, a, deg(a)); return x; }


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* EVALUATE (SINGLE/MULTIPOINT)                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Evaluate at one point                                      */
/*------------------------------------------------------------*/
void eval(Mat<zz_p> & evmat, const Mat<zz_pX> & pmat, zz_p pt);

inline Mat<zz_p> eval(const Mat<zz_pX> & pmat, zz_p pt)
{ Mat<zz_p> evmat; eval(evmat, pmat, pt); return evmat; }


// TODO matrix-wide multipoint evaluation functions?
// (general, geometric, ..)




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* random vector of length n, degree < d                      */
/*------------------------------------------------------------*/
void random(Vec<zz_pX> & pvec, long n, long d);
inline Vec<zz_pX> random_vec_zz_pX(long n, long d)
{ Vec<zz_pX> pvec; random(pvec, n, d); return pvec; }

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d                         */
/*------------------------------------------------------------*/
void random(Mat<zz_pX> & pmat, long m, long n, long d);
inline Mat<zz_pX> random_mat_zz_pX(long n, long m, long d)
{ Mat<zz_pX> pmat; random(pmat, n, m, d); return pmat; }

/*------------------------------------------------------------*/
/* random (m, n) matrix of row degree < rdeg                  */
/*------------------------------------------------------------*/
void random_mat_zz_pX_rdeg(Mat<zz_pX> & pmat, long m, long n, VecLong rdeg);

/*------------------------------------------------------------*/
/* random (m, n) matrix of column degree < cdeg               */
/*------------------------------------------------------------*/
void random_mat_zz_pX_cdeg(Mat<zz_pX> & pmat, long m, long n, VecLong cdeg);
// TODO replace with Vec<long> ??





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
inline Mat<zz_pX> conv(const Mat<zz_p> & coeff)
{ Mat<zz_pX> mat; conv(mat, coeff); return mat; }

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
 * integer `order` must be nonnegative, and it is guaranteed that the vector
 * `matp` has length `order`.
 */
void conv(Vec<Mat<zz_p>> & matp, const Mat<zz_pX> & pmat, const long order);

/** Returns the conversion of a matrix with polynomial entries `pmat`,
 * truncated at the specified `order`, to a polynomial with matrix
 * coefficients. The integer `order` must be nonnegative, and the output
 * vector is guaranteed to have length `order`.
 */
inline Vec<Mat<zz_p>> conv(const Mat<zz_pX> & pmat, const long order)
{ Vec<Mat<zz_p>> matp; conv(matp, pmat, order); return matp; }

/** Converts from polynomial with matrix coefficients `matp`, truncated at the
 * specified `order` (nonnegative integer), to matrix with polynomial entries
 * `pmat`; if `matp` has length 0 then `pmat` is cleared (set to zero without
 * changing its dimensions).
 */
void conv(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp, const long order);

/** Returns the conversion from polynomial with matrix coefficients `matp`,
 * truncated at the specified `order` (nonnegative integer), to matrix with
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
