/**
 * \file nmod_mat_poly.h
 * \brief Matrices with univariate polynomial entries modulo word-size prime, represented as list of matrices
 * \version 0.0
 * \date 2022-06-25
 *
 * \todo doc
 *
 */

#ifndef NMOD_MAT_POLY_H
#define NMOD_MAT_POLY_H

#include <flint/perm.h>
#include <flint/nmod_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

//typedef struct
//{
//    mp_ptr coeffs;
//    slong alloc;
//    slong length;
//    nmod_t mod;
//} nmod_poly_struct;
//
//typedef nmod_poly_struct nmod_poly_t[1];
//
//
//typedef struct
//{
//    mp_limb_t * entries;
//    slong r;
//    slong c;
//    mp_limb_t ** rows;
//    nmod_t mod;
//}
//nmod_mat_struct;

/** Struct Description.
 *
 *  \todo doc
 */
typedef struct
{
    slong alloc;   /**< allocated length */
    slong length;  /**< actual length */
    slong r;       /**< number of rows */
    slong c;       /**< number of columns */
    nmod_t mod;    /**< modulus */
} nmod_mat_poly_struct;

typedef nmod_mat_poly_struct nmod_mat_poly_t[1];



/** void nmod_mat_poly_naive_mul_coef(nmod_mat_t res,
 *				       const nmod_mat_poly_t A,
 *				       const nmod_mat_poly_t B, slong k);
 *
 * A = sum^{deg_A}_{i=0} a_i x^i and B = sum^{deg_B}_{i=0} b_i x^i
 * Compute the coefficient k of the product AB
 * C = AB = sum^{deg_A + deg_B}_{i=0} c_i x^i
 * res = c_k = sum^{k}_{i=0} a_i b_{k-i}
 *
 */
void nmod_mat_poly_naive_mul_coef(nmod_mat_t res,
                                  const nmod_mat_poly_t A,
                                  const nmod_mat_poly_t B,
                                  slong k);

void nmod_mat_poly_print(const nmod_mat_poly_t A);


#ifdef __cplusplus
}
#endif

#endif /* NMOD_MAT_POLY_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
