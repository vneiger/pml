#ifndef MAT_LZZ_PX_ARITH__H
#define MAT_LZZ_PX_ARITH__H

/** \brief Basic arithmetic for univariate polynomial matrices over `zz_p`
 *
 * \file mat_lzz_pX_arith.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-10
 *
 * This file contains the declarations of functions for performing basic
 * arithmetic with polynomial matrices in `Mat<zz_pX>`. This involves addition,
 * subtraction, negation, multiplication by a constant or by a single
 * polynomial, etc.
 *
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* vector addition                                            */
/*------------------------------------------------------------*/

/** @name Addition of polynomial vectors
 *
 * These functions all have the first two of the following parameters,
 * and some of them also have the third:
 * @param[in] a a polynomial vector
 * @param[in] b a polynomial vector
 * @param[out] c a polynomial vector
 *
 * These functions compute `a+b` and either return it, or store it in `c`.
 * They throw an error if the dimensions of `a` and `b` are not compatible.
 *
 *  \todo to be moved in a specific file for vectors
 */
//@{

/** `c = a+b`, `c` may alias `a` or `b` */
void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b);

/** `c = a+b`, `b` constant vector, `c` may alias `a` */
void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b);

/** `c = a+b`, `a` constant vector, `c` may alias `b` */
inline void add(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b)
{ add(c, b, a); }

/** Operator overloading for `a+b`, both polynomial vectors */
inline Vec<zz_pX> operator+(const Vec<zz_pX> & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> x; add(x, a, b); return x; }

/** Operator overloading for `a+b`, `b` constant vector */
inline Vec<zz_pX> operator+(const Vec<zz_pX> & a, const Vec<zz_p> & b)
{ Vec<zz_pX> x; add(x, a, b); return x; }

/** Operator overloading for `a+b`, `a` constant vector */
inline Vec<zz_pX> operator+(const Vec<zz_p> & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> x; add(x, a, b); return x; }

/** Operator overloading for `x += a`, `a` and `x` polynomial vectors */
inline Vec<zz_pX> & operator+=(Vec<zz_pX> & x, const Vec<zz_pX> & a)
{ add(x, x, a); return x; }

/** Operator overloading for `x += a`, `a` constant vector and `x` polynomial
 * vector
 */
inline Vec<zz_pX> & operator+=(Vec<zz_pX> & x, const Vec<zz_p> & a)
{ add(x, x, a); return x; }

//@} // doxygen group:  Addition of polynomial vectors

/** Computes the addition of the left shift `c = a + (b << k)`. The OUT
 * parameter `c` may alias the IN parameter `a`, but not `b`. The integer `k`
 * must be positive. */
void add_LeftShift(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b, long k);
/** Computes the addition of the left shift `c = a + (b << k)`. The OUT
 * parameter `c` may alias the IN parameter `a`, but not `b`. The integer `k`
 * must be positive. */
void add_LeftShift(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long k);

/*------------------------------------------------------------*/
/* vector subtraction                                         */
/*------------------------------------------------------------*/

/** @name Subtraction of polynomial vectors
 *
 * These functions all have the first two of the following parameters,
 * and some of them also have the third:
 * @param[in] a a polynomial vector
 * @param[in] b a polynomial vector
 * @param[out] c a polynomial vector
 *
 * These functions compute `a-b` and either return it, or store it in `c`.
 * They throw an error if the dimensions of `a` and `b` are not compatible.
 *
 *  \todo to be moved in a specific file for vectors
 */
//@{

/** `c = a-b`, `c` may alias `a` or `b` */
void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b);

/** `c = a-b`, `b` constant vector, `c` may alias `a` */
void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b);

/** `c = a-b`, `a` constant vector, `c` may alias `b` */
void sub(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b);

/** Operator overloading for `a-b`, both polynomial vectors */
inline Vec<zz_pX> operator-(const Vec<zz_pX> & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> x; sub(x, a, b); return x; }

/** Operator overloading for `a-b`, `b` constant vector */
inline Vec<zz_pX> operator-(const Vec<zz_pX> & a, const Vec<zz_p> & b)
{ Vec<zz_pX> x; sub(x, a, b); return x; }

/** Operator overloading for `a-b`, `a` constant vector */
inline Vec<zz_pX> operator-(const Vec<zz_p> & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> x; sub(x, a, b); return x; }

/** Operator overloading for `x -= a`, `a` and `x` polynomial vectors */
inline Vec<zz_pX> & operator-=(Vec<zz_pX> & x, const Vec<zz_pX> & a)
{ sub(x, x, a); return x; }

/** Operator overloading for `x -= a`, `a` constant vector and `x` polynomial
 * vector
 */
inline Vec<zz_pX> & operator-=(Vec<zz_pX> & x, const Vec<zz_p> & a)
{ sub(x, x, a); return x; }

//@} // doxygen group:  Subtraction of polynomial vectors



/*------------------------------------------------------------*/
/* matrix addition                                            */
/*------------------------------------------------------------*/

/** @name Addition of polynomial matrices
 *
 * These functions all have the first two of the following parameters,
 * and some of them also have the third:
 * @param[in] a a polynomial matrix
 * @param[in] b a polynomial matrix
 * @param[out] c a polynomial matrix
 *
 * These functions compute `a+b` and either return it, or store it in `c`.
 * They throw an error if the dimensions of `a` and `b` are not compatible.
 *
 */
//@{

/** `c = a+b`, `c` may alias `a` or `b` */
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** `c = a+b`, `b` constant matrix, `c` may alias `a` */
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

/** `c = a+b`, `a` constant matrix, `c` may alias `b` */
inline void add(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{ add(c, b, a); }

/** Operator overloading for `a+b`, both polynomial matrices */
inline Mat<zz_pX> operator+(const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; add(x, a, b); return x; }

/** Operator overloading for `a+b`, `b` constant matrix */
inline Mat<zz_pX> operator+(const Mat<zz_pX> & a, const Mat<zz_p> & b)
{ Mat<zz_pX> x; add(x, a, b); return x; }

/** Operator overloading for `a+b`, `a` constant matrix */
inline Mat<zz_pX> operator+(const Mat<zz_p> & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; add(x, a, b); return x; }

/** Operator overloading for `x += a`, `a` and `x` polynomial matrices */
inline Mat<zz_pX> & operator+=(Mat<zz_pX> & x, const Mat<zz_pX> & a)
{ add(x, x, a); return x; }

/** Operator overloading for `x += a`, `a` constant matrix and `x` polynomial
 * matrix
 */
inline Mat<zz_pX> & operator+=(Mat<zz_pX> & x, const Mat<zz_p> & a)
{ add(x, x, a); return x; }

//@} // doxygen group:  Addition of polynomial matrices


/*------------------------------------------------------------*/
/* matrix subtraction                                         */
/*------------------------------------------------------------*/
/** @name Subtraction of polynomial matrices
 *
 * These functions all have the first two of the following parameters,
 * and some of them also have the third:
 * @param[in] a a polynomial matrix
 * @param[in] b a polynomial matrix
 * @param[out] c a polynomial matrix
 *
 * These functions compute `a-b` and either return it, or store it in `c`.
 * They throw an error if the dimensions of `a` and `b` are not compatible.
 *
 */
//@{

/** `c = a-b`, `c` may alias `a` or `b` */
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** `c = a-b`, `b` constant matrix, `c` may alias `a` */
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

/** `c = a-b`, `a` constant matrix, `c` may alias `b` */
void sub(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

/** Operator overloading for `a-b`, both polynomial matrices */
inline Mat<zz_pX> operator-(const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; sub(x, a, b); return x; }

/** Operator overloading for `a-b`, `b` constant matrix */
inline Mat<zz_pX> operator-(const Mat<zz_pX> & a, const Mat<zz_p> & b)
{ Mat<zz_pX> x; sub(x, a, b); return x; }

/** Operator overloading for `a-b`, `a` constant matrix */
inline Mat<zz_pX> operator-(const Mat<zz_p> & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; sub(x, a, b); return x; }

/** Operator overloading for `x -= a`, `a` and `x` polynomial matrices */
inline Mat<zz_pX> & operator-=(Mat<zz_pX> & x, const Mat<zz_pX> & a)
{ sub(x, x, a); return x; }

/** Operator overloading for `x -= a`, `a` constant matrix and `x` polynomial
 * matrix
 */
inline Mat<zz_pX> & operator-=(Mat<zz_pX> & x, const Mat<zz_p> & a)
{ sub(x, x, a); return x; }

//@} // doxygen group: Subtraction of polynomial matrices

/*------------------------------------------------------------*/
/* negate                                                     */
/*------------------------------------------------------------*/
/** @name Negate
 *
 *  The negation of a polynomial matrix or a polynomial vector `a` is `-a`.
 */
//@{

NTL_OPEN_NNS
/** Computes the polynomial vector `x = -a` */
void negate(Vec<zz_pX> & x, const Vec<zz_pX> & a);
NTL_CLOSE_NNS

/** Computes and returns the polynomial vector `-a` */
inline Vec<zz_pX> operator-(const Vec<zz_pX> & a)
{ Vec<zz_pX> x; NTL::negate(x, a); return x; }

NTL_OPEN_NNS
/** Computes the polynomial matrix `x = -a` */
void negate(Mat<zz_pX> & x, const Mat<zz_pX> & a);
NTL_CLOSE_NNS

/** Computes and returns the polynomial matrix `-a` */
inline Mat<zz_pX> operator-(const Mat<zz_pX> & a)
{ Mat<zz_pX> x; NTL::negate(x, a); return x; }

//@} // doxygen group: Negate

/*------------------------------------------------------------*/
/* scalar and polynomial multiplication for vectors           */
/*------------------------------------------------------------*/

/** @name Multiplication of a vector by a scalar or a polynomial
 *
 *  For a given polynomial vector, the following functions implement the
 *  multiplication by a scalar from `zz_p` or by a polynomial from `zz_pX`.
 */
//@{

/** Computes `c = b*a` for a polynomial vector `a` and a scalar `b`,
 * `c` may alias `a`.
 */
void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_p & b);

/** Computes `c = a*b` for a polynomial vector `b` and a scalar `a`,
 * `c` may alias `b`.
 */
inline void mul(Vec<zz_pX> & c, const zz_p & a, const Vec<zz_pX> & b)
{ mul(c, b, a); }

/** Operator overloading for `b*a` for a polynomial vector `a` and a scalar
 * `b`.
 */
inline Vec<zz_pX> operator*(const Vec<zz_pX> & a, const zz_p & b)
{ Vec<zz_pX> x; mul(x, a, b); return x; }

/** Operator overloading for `a*b` for a polynomial vector `b` and a scalar
 * `a`.
 */
inline Vec<zz_pX> operator*(const zz_p & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> x; mul(x, a, b); return x; }

/** Computes `c = b*a` for a polynomial vector `a` and a polynomial `b`,
 * `c` may alias `a`.
 */
void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_pX & b);

/** Computes `c = a*b` for a polynomial vector `b` and a polynomial `a`,
 * `c` may alias `b`.
 */
inline void mul(Vec<zz_pX> & c, const zz_pX & a, const Vec<zz_pX> & b)
{ mul(c, b, a); }

/** Operator overloading for `b*a` for a polynomial vector `a` and a polynomial
 * `b`.
 */
inline Vec<zz_pX> operator*(const Vec<zz_pX> & a, const zz_pX & b)
{ Vec<zz_pX> x; mul(x, a, b); return x; }

/** Operator overloading for `a*b` for a polynomial vector `b` and a polynomial
 * `a`.
 */
inline Vec<zz_pX> operator*(const zz_pX & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> x; mul(x, a, b); return x; }

//@} // doxygen group: Multiplication of a vector by a scalar or a polynomial

/*------------------------------------------------------------*/
/* scalar and polynomial multiplication                       */
/*------------------------------------------------------------*/

/** @name Multiplication of a matrix by a scalar or a polynomial
 *
 *  For a given polynomial matrix, the following functions implement the
 *  multiplication by a scalar from `zz_p` or by a polynomial from `zz_pX`.
 */
//@{

/** Computes `c = b*a` for a polynomial matrix `a` and a scalar `b`,
 * `c` may alias `a`.
 */
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b);

/** Computes `c = a*b` for a polynomial matrix `b` and a scalar `a`,
 * `c` may alias `b`.
 */
inline void mul(Mat<zz_pX> & c, const zz_p & a, const Mat<zz_pX> & b)
{ mul(c, b, a); }

/** Operator overloading for `b*a` for a polynomial matrix `a` and a scalar
 * `b`.
 */
inline Mat<zz_pX> operator*(const Mat<zz_pX> & a, const zz_p & b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

/** Operator overloading for `a*b` for a polynomial matrix `b` and a scalar
 * `a`.
 */
inline Mat<zz_pX> operator*(const zz_p & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

/** Computes `c = b*a` for a polynomial matrix `a` and a polynomial `b`,
 * `c` may alias `a`.
 */
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_pX & b);

/** Computes `c = a*b` for a polynomial matrix `b` and a polynomial `a`,
 * `c` may alias `b`.
 */
inline void mul(Mat<zz_pX> & c, const zz_pX & a, const Mat<zz_pX> & b)
{ mul(c, b, a); }

/** Operator overloading for `b*a` for a polynomial matrix `a` and a polynomial
 * `b`.
 */
inline Mat<zz_pX> operator*(const Mat<zz_pX> & a, const zz_pX & b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

/** Operator overloading for `a*b` for a polynomial matrix `b` and a polynomial
 * `a`.
 */
inline Mat<zz_pX> operator*(const zz_pX & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

//@} // doxygen group: Multiplication of a matrix by a scalar or a polynomial

/*------------------------------------------------------------*/
/* multiplication by a constant matrix                        */
/*------------------------------------------------------------*/
/** @name Multiplication by a constant matrix
 *
 *  The following functions implement the multiplication of a polynomial matrix
 *  by a constant matrix. They throw an error if the dimensions are not
 *  compatible.
 *
 * \todo
 *   - try the following (more efficient?): rely on complete linearization as
 *   constant matrices. That is, expand 'a' as single (or several?) constant
 *   matrix 'cmat', and compute b*cmat, and retrieve back the entries in 'c'
 *   (preliminary work: compute cdeg(a))
 *   - make this depend on conv?
 *   - more generally, how to proceed should depend on the dimensions (for
 *   example, it is not the same story if the non-constant mat is square and if
 *   it is a single vector). Thresholds may be needed.
 *
 */
//@{

/** Computes `c = a*b` for a polynomial matrix `a` and a constant matrix `b`,
 * `c` may alias `a`
 */
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

/** Computes `c = a*b` for a constant matrix `a` and a polynomial matrix `b`,
 * `c` may alias `b`
 */
void mul(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

/** Operator overloading for `a*b` for a polynomial matrix `a` and a constant
 * matrix `b`
 */
inline Mat<zz_pX> operator*(const Mat<zz_pX> & a, const Mat<zz_p> & b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

/** Operator overloading for `a*b` for a polynomial matrix `b` and a constant
 * matrix `a`
 */
inline Mat<zz_pX> operator*(const Mat<zz_p> & a, const Mat<zz_pX> & b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

/** Operator overloading for `x *= a` for a polynomial constant matrix `a` and
 * a polynomial matrix `x`
 */
inline Mat<zz_pX> & operator*=(Mat<zz_pX> & x, const Mat<zz_p> & a)
{ mul(x, x, a); return x; }

//@} // doxygen group: Multiplication by a constant matrix

#endif /* ifndef MAT_LZZ_PX_ARITH__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
