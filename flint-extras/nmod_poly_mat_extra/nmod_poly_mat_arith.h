#ifndef NMOD_POLY_MAT_ARITH_H
#define NMOD_POLY_MAT_ARITH_H

/** \brief Basic arithmetic for univariate polynomial matrices over `zz_p`
 *
 * \file mat_lzz_pX_arith.h
 * \version 0.0
 * \date 2022-06-27
 *
 * This file contains the declarations of functions for performing basic
 * arithmetic with univariate polynomial matrices. This involves addition,
 * subtraction, negation, multiplication by a constant or by a single
 * polynomial, etc.
 *
 */

#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

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
 * These functions compute `a+b` and store it in `c`.
 * The dimensions of `a` and `b` are not checked to be compatible.
 *
 *  \todo to be moved in a specific file for vectors?
 */
//@{

/** `c = a+b`, `c` may alias `a` or `b` */
// TODO
//void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b);

/** `c = a+b`, `b` constant vector, `c` may alias `a` */
// TODO
//void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b);

/** `c = a+b`, `a` constant vector, `c` may alias `b` */
// TODO
//inline void add(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b)
//{ add(c, b, a); }

/** Computes the addition of the left shift `c = a + (b << k)`. The OUT
 * parameter `c` may alias the IN parameter `a`, but not `b`. The integer `k`
 * must be positive. */
// TODO
//void add_LeftShift(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b, long k);

/** Computes the addition of the left shift `c = a + (b << k)`. The OUT
 * parameter `c` may alias the IN parameter `a`, but not `b`. The integer `k`
 * must be positive. */
// TODO
//void add_LeftShift(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long k);

//@} // doxygen group:  Addition of polynomial vectors

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
 * These functions compute `a-b` and store it in `c`.
 * The dimensions of `a` and `b` are not checked to be compatible.
 *
 *  \todo to be moved in a specific file for vectors
 */
//@{

/** `c = a-b`, `c` may alias `a` or `b` */
// TODO
//void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b);

/** `c = a-b`, `b` constant vector, `c` may alias `a` */
// TODO
//void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b);

/** `c = a-b`, `a` constant vector, `c` may alias `b` */
// TODO
//void sub(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b);

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
 * These functions compute `a+b` and store it in `c`.
 * The dimensions of `a` and `b` are not checked to be compatible.
 *
 */
//@{

/** `c = a+b`, `c` may alias `a` or `b` */
// TODO
//void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** `c = a+b`, `b` constant matrix, `c` may alias `a` */
// TODO
//void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

/** `c = a+b`, `a` constant matrix, `c` may alias `b` */
// TODO
//inline void add(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
//{ add(c, b, a); }

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
 * These functions compute `a-b` and store it in `c`.
 * The dimensions of `a` and `b` are not checked to be compatible.
 *
 */
//@{

/** `c = a-b`, `c` may alias `a` or `b` */
// TODO
//void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** `c = a-b`, `b` constant matrix, `c` may alias `a` */
// TODO
//void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

/** `c = a-b`, `a` constant matrix, `c` may alias `b` */
// TODO
//void sub(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

//@} // doxygen group: Subtraction of polynomial matrices

/*------------------------------------------------------------*/
/* negate                                                     */
/*------------------------------------------------------------*/
/** @name Negate
 *
 *  The negation of a polynomial matrix or a polynomial vector `a` is `-a`.
 */
//@{

/** Computes the polynomial vector `x = -a` */
// TODO
//void negate(Vec<zz_pX> & x, const Vec<zz_pX> & a);

/** Computes the polynomial matrix `x = -a` */
// TODO
//void negate(Mat<zz_pX> & x, const Mat<zz_pX> & a);

//@} // doxygen group: Negate

/*------------------------------------------------------------*/
/* scalar and polynomial multiplication for vectors           */
/*------------------------------------------------------------*/

/** @name Multiplication of a vector by a scalar or a polynomial
 *
 *  For a given polynomial vector, the following functions implement the
 *  multiplication by a scalar from `nmod` or by a polynomial from `nmod_poly`.
 */
//@{

/** Computes `c = b*a` for a polynomial vector `a` and a scalar `b`,
 * `c` may alias `a`.
 */
// TODO
//void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_p & b);

/** Computes `c = a*b` for a polynomial vector `b` and a scalar `a`,
 * `c` may alias `b`.
 */
// TODO
//inline void mul(Vec<zz_pX> & c, const zz_p & a, const Vec<zz_pX> & b)
//{ mul(c, b, a); }

/** Computes `c = b*a` for a polynomial vector `a` and a polynomial `b`,
 * `c` may alias `a`.
 */
// TODO
//void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_pX & b);

/** Computes `c = a*b` for a polynomial vector `b` and a polynomial `a`,
 * `c` may alias `b`.
 */
// TODO
//inline void mul(Vec<zz_pX> & c, const zz_pX & a, const Vec<zz_pX> & b)
//{ mul(c, b, a); }

//@} // doxygen group: Multiplication of a vector by a scalar or a polynomial

/*------------------------------------------------------------*/
/* scalar and polynomial multiplication                       */
/*------------------------------------------------------------*/

/** @name Multiplication of a matrix by a scalar or a polynomial
 *
 *  For a given polynomial matrix, the following functions implement the
 *  multiplication by a scalar from `nmod` or by a polynomial from `nmod_poly`.
 */
//@{

/** Computes `c = b*a` for a polynomial matrix `a` and a scalar `b`,
 * `c` may alias `a`.
 */
// TODO
//void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b);

/** Computes `c = a*b` for a polynomial matrix `b` and a scalar `a`,
 * `c` may alias `b`.
 */
// TODO
//inline void mul(Mat<zz_pX> & c, const zz_p & a, const Mat<zz_pX> & b)
//{ mul(c, b, a); }

/** Computes `c = b*a` for a polynomial matrix `a` and a polynomial `b`,
 * `c` may alias `a`.
 */
// TODO
//void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_pX & b);

/** Computes `c = a*b` for a polynomial matrix `b` and a polynomial `a`,
 * `c` may alias `b`.
 */
// TODO
//inline void mul(Mat<zz_pX> & c, const zz_pX & a, const Mat<zz_pX> & b)
//{ mul(c, b, a); }

//@} // doxygen group: Multiplication of a matrix by a scalar or a polynomial

/*------------------------------------------------------------*/
/* multiplication by a constant matrix                        */
/*------------------------------------------------------------*/
/** @name Multiplication by a constant matrix
 *
 *  The following functions implement the multiplication of a polynomial matrix
 *  by a constant matrix. The dimensions are not checked to be compatible.
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
// TODO
//void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

/** Computes `c = a*b` for a constant matrix `a` and a polynomial matrix `b`,
 * `c` may alias `b`
 */
// TODO
//void mul(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

//@} // doxygen group: Multiplication by a constant matrix

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_ARITH_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
