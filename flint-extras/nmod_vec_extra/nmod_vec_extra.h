#ifndef __NMOD_VEC_EXTRA__H
#define __NMOD_VEC_EXTRA__H

/** \brief Extra functions for vectors over Z/nZ
 *
 * \file nmod_vec_extra.h
 * \version 0.0
 * \date 2023-01-23
 *
 * Some functions to deal with vectors over `nmod`.
 *
 */

#include "nmod_vec_extra.h"
#include <flint/machine_vectors.h>
#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif


// TODO augment/add documentation

/** Random */
// to be completed: random sparse? randtest small? randtest nonzero ? ...

/** Fills the entries `0`, .., `len-1` of vector with uniformly random entries.
 * Vector must already be allocated with length at least `len`. */
void _nmod_vec_rand(mp_ptr vec,
            		flint_rand_t state,
            		slong len,
            		nmod_t mod);


/*--------------------------------------------------------------*/
/* vector of n consecutive primes of exactly s bits             */
/*--------------------------------------------------------------*/
void nmod_vec_primes(mp_ptr v, slong n, mp_bitcnt_t s);

/**********************************************************************
*                            DOT PRODUCT                             *
**********************************************************************/


/* ------------------------------------------------------------ */
/* v1 and v2 have length at least len, len < 2^FLINT_BITS      */
/* all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/* all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/* computes sum(v1[i]*v2[i], 0 <= i < len)                      */
/* stores the result in 3 limbs of res                          */
/* ------------------------------------------------------------ */
void nmod_vec_integer_dot_product(mp_ptr res,
                                  mp_srcptr v1, mp_srcptr v2,
                                  ulong len, ulong max_bits1, ulong max_bits2);

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** does not assume input is reduced modulo mod.n                */
/*  ------------------------------------------------------------ */
mp_limb_t nmod_vec_dot_product(mp_srcptr v1, mp_srcptr v2,
                               ulong len, ulong max_bits1, ulong max_bits2,
                               nmod_t mod);
// note: version split16 interesting on recent laptop (gcc does some vectorization)
// limited to nbits <= ~31 (bound to be better analyzed, numterms)
mp_limb_t _nmod_vec_dot_product_2_split16(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod);
// note: version split26 interesting (beyond 30-31 bits) on recent laptop (gcc does some vectorization)
// limited to nbits <= ~52 (TODO bound to be better analyzed, numterms; potential fixes in code needed)
mp_limb_t _nmod_vec_dot_product_2_split26(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod);

/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** returns dot(a, b)                                         */
/** power_two = 2^45 mod p, pinv = 1/p                        */
/*------------------------------------------------------------*/
mp_limb_t _nmod_vec_dot_small_modulus(mp_ptr a, mp_ptr b, ulong len,
                                      mp_limb_t power_two,
                                      vec1d p, vec1d pinv);

/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** res[0] = dot(a1, b), res[1] = dot(a2, b)                  */
/** power_two = 2^45 mod p, p2 = (p,p), pinv2 = (1/p,1/p)     */
/*------------------------------------------------------------*/
void _nmod_vec_dot2_small_modulus(mp_ptr res,
                                  mp_ptr a1, mp_ptr a2, mp_ptr b, ulong len,
                                  mp_limb_t power_two,
                                  vec2d p2, vec2d pinv2);

/** Several dot products with same left operand, as in vector-matrix product.
 *
 * . u has length at least len, len < 2^FLINT_BITS
 * . all entries of u  have <= max_bits_u bits <= FLINT_BITS
 * . v points to at least len vectors v[0],..,v[len-1] each of length
 * at least k, with entries of <= max_bits_v bits <= FLINT_BITS
 * . computes uv[j] = sum(u[i]*v[i][j], 0 <= i < len) modulo mod.n,
 * for 0 <= j < k
 * . does not assume input entries are reduced modulo mod.n
 *
 * \todo variants 8_16, 16_32, 32_32
 * \todo then, thresholds needed: on two different machines, 
 * - one (recent, avx512) is most often faster with v8_8 (but sometimes v8_32 better, e.g. len=128 k=16)
 * - another one (2020, no avx512) is most often faster with v8_32 and v16_16 > v8_8
 * \todo eventually, to be compared to a good matrix-matrix product with a single row in left operand!
 */
void nmod_vec_dot_product_multi(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                ulong len, ulong k,
                                ulong max_bits_u, ulong max_bits_v,
                                nmod_t mod);
void _nmod_vec_dot_product_multi_1_v1_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                       ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_1_v8_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                       ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_1_v8_32(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                       ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_1_v16_16(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                       ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v1_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v4_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v8_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v4_32(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod);
// limited to nbits <= ~52 (TODO bound to be better analyzed, numterms; potential fixes in code needed)
void _nmod_vec_dot_product_multi_2_split26(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod);

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD_VEC_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
