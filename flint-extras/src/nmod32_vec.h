#ifndef __NMOD32_VEC__H
#define __NMOD32_VEC__H

/** \brief functions for vectors over Z/nZ, n < 2**32
 *
 * \file nmod_vec_extra.h
 *
 * Some functions to deal with vectors over `nmod32`.
 *
 * Assumes FLINT_BITS == 64. Assumes unsigned int is 32 bits.
 *
 */

#include <flint/flint.h>
#include <flint/machine_vectors.h>
#include <flint/nmod.h> // for NMOD_RED

#define HAVE_AVX512 1   // TODO handle AVX flags

#ifdef __AVX2__
#define HAS_AVX2
#endif

#ifdef HAS_AVX2  // GV


// functions below accumulate 8 terms:
//  -> modulus <= DOT2_ACC8_MAX_MODULUS
//  -> len <= DOT2_ACC8_MAX_LEN
// (see bottom of flint/src/nmod_vec/dot.c for more details)
#define DOT2_ACC8_MAX_MODULUS UWORD(1515531528)
#define DOT2_ACC8_MAX_LEN UWORD(380368697)

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned int uint; // assumes this is 32 bits on this machine
typedef uint * n32_ptr;
typedef const uint * n32_srcptr;

FLINT_INLINE
n32_ptr _nmod32_vec_init(slong len)
{
    return (n32_ptr) flint_malloc(len * sizeof(uint));
}

FLINT_INLINE
void _nmod32_vec_clear(n32_ptr vec)
{
   flint_free(vec);
}

// avx2 horizontal sum
FLINT_FORCE_INLINE ulong _mm256_hsum(__m256i a) {
    __m256i a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    __m256i sum_lo = _mm256_add_epi64(a, a_hi);
    __m128i sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (ulong) _mm_cvtsi128_si64(sum);
}

#if HAVE_AVX512
// avx512 horizontal sum
FLINT_FORCE_INLINE ulong _mm512_hsum(__m512i a) {
    return _mm512_reduce_add_epi64(a);
}
#endif

/**********************************************************************
*                            DOT PRODUCT                             *
**********************************************************************/

// WARNING: everything below assumes both
//     modulus <= DOT2_ACC8_MAX_MODULUS  (about 2**30.5)
//     len <= DOT2_ACC8_MAX_LEN          (about 2**28.5)
// behaviour is undefined if that is not satisfied

uint _nmod32_vec_dot_split(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);

uint _nmod32_vec_dot_split_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#if HAVE_AVX512
uint _nmod32_vec_dot_split_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif

// duplicate msolve's matrix-vector product (for single dot prod)
uint _nmod32_vec_dot_msolve_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, uint PRIME);

// ifma variant
// TODO unfinished, requires len multiple of 16 or 32
uint _nmod32_vec_dot_ifma_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#if HAVE_AVX512
uint _nmod32_vec_dot_ifma_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif

/**********************************************************************
*                          k-DOT PRODUCT                             *
**********************************************************************/

// similar to dot_split, computes res0 = vec1 * vec2_0 and res1 = vec1 * vec2_1
void _nmod32_vec_dot2_split(uint * res0, uint * res1,
                            n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                            slong len, nmod_t mod, ulong pow2_precomp);
void _nmod32_vec_dot2_split_avx2(uint * res0, uint * res1,
                                 n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                 slong len, nmod_t mod, ulong pow2_precomp);

// similar to dot_split, computes res[i] = vec1 * vec2_i for i = 0, 1, 2
void _nmod32_vec_dot3_split_avx2(uint * res, n32_srcptr vec1,
                                 n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                 slong len, nmod_t mod, ulong pow2_precomp);

#if HAVE_AVX512
void _nmod32_vec_dot2_split_avx512(uint * res0, uint * res1,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                   slong len, nmod_t mod, ulong pow2_precomp);

// similar to dot_split, computes res[i] = vec1 * vec2_i for i = 0, 1, 2, 3
void _nmod32_vec_dot3_split_avx512(uint * res,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                   slong len, nmod_t mod, ulong pow2_precomp);
#endif

/**********************************************************************
*                         MULTI DOT PRODUCT                          *
**********************************************************************/

/** Several dot products with fixed operand, as in matrix-vector product.
 *
 * . vec: vector of length len
 * . mat: matrix with nrows rows of length >= len, stored contiguously
 * . stride: indicates how far away mat[i,j] and mat[i+1,j] are in memory
 * . computes the matrix-vector product
 *      mv[i] = sum(mat[i,j]*vec[j], 0 <= j < len) modulo mod.n, for 0 <= i < nrows
 */

// naive via dot_split
void _nmod32_vec_mdot_split(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                            slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                 slong nrows, slong len, slong stride, nmod_t mod);

#if HAVE_AVX512
void _nmod32_vec_mdot_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                   slong nrows, slong len, slong stride, nmod_t mod);
#endif

// duplicate msolve's matrix-vector product
void _nmod32_vec_mdot_msolve_via_dot_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                          slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot_msolve_native_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                         slong nrows, slong len, slong stride, nmod_t mod);

// like dot_split but handling two rows at a time
void _nmod32_vec_mdot2_split(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                             slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot2_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot3_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod);

#if HAVE_AVX512
void _nmod32_vec_mdot2_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot3_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod);
#endif

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD32_VEC__H

#endif // GV HAS_AVX2

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
