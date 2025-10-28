/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

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

/** Random */

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

/***************
*  inversion  *
***************/

/* computes res[i] = inv(v[i]) for i = 0 ... len */
/* res cannot alias inv */
void _nmod_vec_inv(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod);
void _nmod_vec_inv2(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod);

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
/* identical to flint's _nmod_vec_dot2_split,                 */
/* but handles 2 vectors at a time:                           */
/*  res[0] = dot(vec10, vec2), res[1] = dot(vec11, vec2)      */
/* FIXME constraints on input?                                */
/*------------------------------------------------------------*/
#if FLINT_BITS == 64
void _nmod_vec_2dot2_split(nn_ptr res, nn_srcptr vec10, nn_srcptr vec11, nn_srcptr vec2,
                           slong len, nmod_t mod, ulong pow2_precomp);
#endif  /* FLINT_BITS == 64 */


/*------------------------------------------------------------*/
/* DOT PRODUCT: DRAFT / EXPERIMENTS                           */
/*------------------------------------------------------------*/

// NOTES
// attempt: vectorized dot product in <= 32 bits with madd_epi16 (AVX2)
//    issue: due to madd we split at 15bits instead of 16, and are limited at 30 bits max
//    and for 30 bits or less this does not bring much compared to the avx2 _nmod_vec_dot2_split already in FLINT
// attempt: (not tried) vectorization with dp VNNI (less widely available), should bring the same issue of limiting to 30 bits

// note: version split26_avx interesting (beyond 30-31 bits) on zen4 laptop
// --> speed-up about 2 for lengths a few 100s
// (FIXME analyze correctness bounds)
ulong _nmod_vec_dot_product_split26(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
#ifdef PML_HAVE_AVX2
    ulong _nmod_vec_dot_product_split26_avx(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
#endif // PML_HAVE_AVX2

// FIXME in progress: ifma attempt
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
