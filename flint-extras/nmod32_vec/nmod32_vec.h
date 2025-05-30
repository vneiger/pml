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


/**********************************************************************
*                            DOT PRODUCT                             *
**********************************************************************/

#define HAVE_AVX512 0   // TODO handle AVX flags

// duplicates flint's dot2_split based on avx2

// accumulates 4 terms, limit 2**31
uint _nmod32_vec_dot_split(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);

// accumulate 4 or 8 terms
uint _nmod32_vec_dot_split_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#if HAVE_AVX512   // TODO handle AVX flags
uint _nmod32_vec_dot_split_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif

// duplicate msolve's matrix-vector product (for single dot prod)
uint _nmod32_vec_dot_msolve_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, uint PRIME);

// ifma variant
// TODO unfinished, requires len multiple of 16 or 32
uint _nmod32_vec_dot_ifma_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#if HAVE_AVX512   // TODO handle AVX flags
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
#if HAVE_AVX512   // TODO handle AVX flags
void _nmod32_vec_dot2_split_avx512(uint * res0, uint * res1,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                   slong len, nmod_t mod, ulong pow2_precomp);
#endif

// similar to dot_split, computes res_i = vec1 * vec2_i for i = 0, 1, 2, 3, ...
void _nmod32_vec_dot3_split_avx2(uint * res0, uint * res1, uint * res2,
                                 n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                 slong len, nmod_t mod, ulong pow2_precomp);
#if HAVE_AVX512   // TODO handle AVX flags
void _nmod32_vec_dot4_split_avx512(uint * res0, uint * res1, uint * res2, uint * res3,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2, n32_srcptr vec2_3,
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
#if HAVE_AVX512   // TODO handle AVX flags
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
void _nmod32_vec_mdot2_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot3_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod);
#if HAVE_AVX512   // TODO handle AVX flags
void _nmod32_vec_mdot4_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod);
#endif

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD32_VEC__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
