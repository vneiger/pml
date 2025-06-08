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

#include <flint/flint.h>
#include <flint/machine_vectors.h>
#include <flint/nmod_types.h>
#include <flint/nmod.h> // for NMOD_RED

#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif


// TODO augment/add documentation

/** Random */
// to be completed: random sparse? randtest small? randtest nonzero ? ...

/** Fills the entries `0`, .., `len-1` of vector with uniformly random entries.
 * Vector must already be allocated with length at least `len`. */
void _nmod_vec_rand(nn_ptr vec,
            		flint_rand_t state,
            		slong len,
            		nmod_t mod);


/*--------------------------------------------------------------*/
/* vector of n consecutive primes of exactly s bits             */
/*--------------------------------------------------------------*/
void nmod_vec_primes(nn_ptr v, slong n, flint_bitcnt_t s);

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
void nmod_vec_integer_dot_product(nn_ptr res,
                                  nn_srcptr v1, nn_srcptr v2,
                                  ulong len, ulong max_bits1, ulong max_bits2);

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** does not assume input is reduced modulo mod.n                */
/*  ------------------------------------------------------------ */
ulong nmod_vec_dot_product_unbalanced(nn_srcptr v1, nn_srcptr v2,
                                      ulong len, ulong max_bits1, ulong max_bits2,
                                      nmod_t mod);


/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** res[0] = dot(a1, b), res[1] = dot(a2, b)                  */
/** power_two = 2^45 mod p, p2 = (p,p), pinv2 = (1/p,1/p)     */
/*------------------------------------------------------------*/
/**************
*  NOTES:
*  - on zen4, do not see an advantage for dot2_small_modulus
*  - on intel icelake server, small advantage (<15%) for lengths about 200 and more
**************/
#if PML_HAVE_MACHINE_VECTORS
void _nmod_vec_dot2_small_modulus(nn_ptr res, nn_ptr a1, nn_ptr a2, nn_ptr b, ulong len,
                                  ulong power_two, vec2d p2, vec2d pinv2);
#endif


/*------------------------------------------------------------*/
/* DRAFT / EXPERIMENTS                                        */
/*------------------------------------------------------------*/

// NOTES
// attempt: vectorized dot product in <= 32 bits with madd_epi16 (AVX2)
//    issue: due to madd we split at 15bits instead of 16, and are limited at 30 bits max
//    and for 30 bits or less this does not bring much compared to the avx2 _nmod_vec_dot2_split already in FLINT
// attempt: (not tried) vectorization with dp VNNI (less widely available), should bring the same issue of limiting to 30 bits


// TODO incorporate FLINT (seems interesting in all cases, see below):
// faster nmod_vec for modulus close to 32 bits (almost as fast as AVX-based version for modulus < 2**30.5 bits)
// (correctness bound related to how much we can accumulate in high part: probably same bound as for AVX already in flint?)
#if defined(__AVX2__) && FLINT_BITS == 64
ulong _nmod_vec_dot2_half_avx(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
#endif
// some timings on zen4 with modulus == 32 bits close to 2**32:
//            bit/len 50      100     200     400     600     800     1000    2000    4000    8000    16000   50000   1000000
// cf oldhalf 32hi    1.8e-08 3.3e-08 6.2e-08 1.2e-07 1.8e-07 2.4e-07 2.9e-07 5.9e-07 1.2e-06 2.4e-06 4.6e-06 1.4e-05 3.2e-04
// cf newhalf 32hi    1.1e-08 1.5e-08 2.7e-08 4.1e-08 5.9e-08 7.4e-08 9.5e-08 1.9e-07 4.4e-07 9.1e-07 1.7e-06 6.4e-06 2.0e-04
// cf new_int 32hi    1.0e-08 1.6e-08 2.4e-08 4.2e-08 6.2e-08 8.0e-08 9.5e-08 2.0e-07 3.8e-07 7.5e-07 1.5e-06 4.7e-06 9.7e-05
// cu oldhalf 32hi    2.0e-08 3.5e-08 6.3e-08 1.4e-07 2.1e-07 2.9e-07 3.4e-07 7.3e-07 1.4e-06 3.2e-06 5.4e-06 1.8e-05 3.5e-04
// cu newhalf 32hi    1.3e-08 1.9e-08 3.0e-08 9.3e-08 1.5e-07 2.3e-07 2.7e-07 6.7e-07 1.3e-06 2.8e-06 5.2e-06 1.7e-05 3.3e-04
// cu new_int 32hi    1.1e-08 1.7e-08 2.8e-08 4.6e-08 6.8e-08 9.8e-08 1.5e-07 2.9e-07 7.2e-07 1.3e-06 2.7e-06 8.5e-06 1.6e-04
// same on Groebner (icelake server)
//            bit/len 50      100     200     400     600     800     1000    2000    4000    8000    16000   50000   1000000
// cf oldhalf 32hi    3.2e-08 5.4e-08 1.0e-07 2.0e-07 2.9e-07 3.9e-07 4.9e-07 9.7e-07 2.1e-06 4.2e-06 8.4e-06 2.7e-05 5.7e-04
// cf newhalf 32hi    2.0e-08 2.9e-08 5.0e-08 8.0e-08 1.1e-07 1.3e-07 1.7e-07 3.2e-07 9.2e-07 1.8e-06 3.6e-06 1.3e-05 4.9e-04
// cf new_int 32hi    1.8e-08 2.5e-08 4.4e-08 7.3e-08 1.1e-07 1.4e-07 1.7e-07 3.3e-07 6.5e-07 1.3e-06 2.6e-06 8.7e-06 2.4e-04
// cu oldhalf 32hi    3.7e-08 6.7e-08 1.2e-07 2.3e-07 4.1e-07 7.5e-07 1.0e-06 2.2e-06 4.6e-06 9.0e-06 1.8e-05 5.5e-05 1.0e-03
// cu newhalf 32hi    2.5e-08 4.7e-08 9.4e-08 2.0e-07 3.4e-07 6.1e-07 9.0e-07 2.1e-06 4.2e-06 8.4e-06 1.6e-05 5.1e-05 9.8e-04
// cu new_int 32hi    1.8e-08 3.0e-08 5.1e-08 9.5e-08 1.4e-07 1.9e-07 2.4e-07 7.5e-07 2.1e-06 4.4e-06 8.9e-06 2.7e-05 5.2e-04
// same on machaon (skylake)
//            bit/len 50      100     200     400     600     800     1000    2000    4000    8000    16000   50000   1000000
// cf oldhalf 32hi    2.8e-08 4.2e-08 8.3e-08 1.6e-07 2.3e-07 3.3e-07 3.8e-07 7.8e-07 1.8e-06 3.4e-06 7.4e-06 2.3e-05
// cf newhalf 32hi    2.3e-08 2.8e-08 5.1e-08 7.1e-08 1.0e-07 1.4e-07 1.6e-07 3.3e-07 7.6e-07 1.6e-06 3.7e-06 1.4e-05
// cf new_int 32hi    2.5e-08 3.3e-08 4.3e-08 7.6e-08 1.0e-07 1.3e-07 1.7e-07 3.2e-07 6.6e-07 1.2e-06 2.9e-06 1.1e-05
// cu oldhalf 32hi    3.0e-08 5.5e-08 9.3e-08 2.6e-07 4.2e-07 6.3e-07 7.6e-07 1.6e-06 3.2e-06 6.6e-06 1.3e-05 4.1e-05
// cu newhalf 32hi    3.0e-08 3.8e-08 7.5e-08 1.7e-07 3.7e-07 5.1e-07 5.9e-07 1.3e-06 2.6e-06 5.3e-06 1.1e-05 3.5e-05
// cu new_int 32hi    2.2e-08 3.8e-08 4.8e-08 9.8e-08 1.5e-07 2.6e-07 3.6e-07 7.7e-07 1.7e-06 3.5e-06 7.1e-06 2.1e-05

// note: version split26_avx interesting (beyond 30-31 bits) on zen4 laptop
// --> speed-up about 2 for lengths a few 100s
// (TODO analyze correctness bounds)
ulong _nmod_vec_dot_product_split26(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
#ifdef PML_HAVE_AVX2
    ulong _nmod_vec_dot_product_split26_avx(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
#endif // PML_HAVE_AVX2


// TODO in progress: ifma attempt
#if PML_HAVE_AVX512
ulong _nmod_vec_dot_product_ifma256(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
ulong _nmod_vec_dot_product_ifma512(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
#endif  /* PML_HAVE_AVX512 */


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
 * \todo do we want to write avx512 versions... ? below, attempts at making
 * the compiler do it for us; get some improvements, but not gigantic
 *
 * for multi_1:
 * \todo variants 8_16, 16_32, 32_32
 * \todo then, thresholds needed: on two different machines, 
 * - one (recent, avx512) is most often faster with v8_8 (but sometimes v8_32 better, e.g. len=128 k=16)
 * - another one (2020, no avx512) is most often faster with v8_32 and v16_16 > v8_8
 * - split 26 is in draft version, not properly checked for overflow
 * \todo eventually, to be compared to a good matrix-matrix product with a single row in left operand!
 *
 * for multi_2:
 * - 8_8, 8_32, 16_16 look similar, and on avx512 machine, are faster than 1_8 or basic
 * - split 26 (no blocking yet; might be attempted) with avx512 is faster than the above with a factor sometimes > 2
 * - split 26 is in draft version, not properly checked for overflow
 * - again, more benchmarking and threshold needed
 *
 * for multi_3:
 * - at the moment, no attempt at blocking or other things; will depend on what happens for multi_{1,2}
 */
void nmod_vec_dot_product_multi(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                ulong len, ulong k,
                                ulong max_bits_u, ulong max_bits_v,
                                nmod_t mod);

void _nmod_vec_dot_product_multi_1_v1_8(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                        ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_1_v8_8(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                        ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_1_v8_32(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                         ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_1_v16_16(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                          ulong len, ulong k, nmod_t mod);

void _nmod_vec_dot_product_multi_2_v1_8(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                        ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v4_8(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                        ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v8_8(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                        ulong len, ulong k, nmod_t mod);
void _nmod_vec_dot_product_multi_2_v4_32(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                         ulong len, ulong k, nmod_t mod);
// limited to nbits <= ~52 (TODO bound to be better analyzed, numterms; potential fixes in code needed)
void _nmod_vec_dot_product_multi_2_split26(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
                                           ulong len, ulong k, nmod_t mod);

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD_VEC_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
