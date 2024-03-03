#ifndef NMOD_POLY_MAT_ARITH_H
#define NMOD_POLY_MAT_ARITH_H

/** \brief Basic arithmetic for univariate polynomial matrices over `nmod`
 *
 * \file nmod_poly_mat_arith.h
 * \version 0.1
 * \date 2023-10-12
 *
 * This file contains the declarations of functions for performing basic
 * arithmetic with univariate polynomial matrices. This involves addition,
 * subtraction, negation, multiplication by a constant or by a single
 * polynomial, etc.
 *
 */

#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif

// TODO missing in FLINT:
// - add/sub nmod_mat with nmod_poly_mat 
// - add/sub x**k * nmod_mat to nmod_poly_mat, or more generally combined add+shift of polymats?
// - left/right multiplication by a constant matrix

/*------------------------------------------------------------*/
/* multiplication by a constant matrix                        */
/*------------------------------------------------------------*/
/** @name Multiplication by a constant matrix
 *
 *  The following functions implement the multiplication of a polynomial matrix
 *  by a constant matrix. The dimensions are not checked to be compatible.
 *
 * \todo not yet implemented, at least two strategies:
 *   - rely on complete linearization as
 *   constant matrices. That is, expand 'a' as single (or several?) constant
 *   matrix 'cmat', and compute b*cmat, and retrieve back the entries in 'c'
 *   (preliminary work: compute cdeg(a))
 *   - or convert to matpoly, multiply, convert back
 *   -> what is faster should depend on the dimensions (for example, it is not
 *   the same story if the non-constant mat is square and if it is a single
 *   vector), so thresholds may be needed.
 *
 */
//@{

/** Computes `c = a*b` for a polynomial matrix `a` and a constant matrix `b`,
 * `c` may alias `a`
 * \todo not implemented yet
 */
//void nmod_poly_mat_mul_right_nmod_mat(nmod_poly_mat_t & c, const nmod_poly_mat_t & a, const nmod_mat_t & b);

/** Computes `c = a*b` for a constant matrix `a` and a polynomial matrix `b`,
 * `c` may alias `b`
 * \todo not implemented yet
 */
//void nmod_poly_mat_mul_left_nmod_mat(nmod_poly_mat_t & c, const nmod_mat_t & a, const nmod_poly_mat_t & b);

//@} // doxygen group: Multiplication by a constant matrix

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_ARITH_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
