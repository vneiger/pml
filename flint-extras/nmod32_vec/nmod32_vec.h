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

// duplicates flint's dot2_split based on avx2
#define HAVE_AVX512 1   // TODO handle AVX flags
uint _nmod32_vec_dot2_split_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#if HAVE_AVX512   // TODO handle AVX flags
uint _nmod32_vec_dot2_split_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif

// faster nmod_vec for modulus close to 32 bits
// (correctness bound related to how much we can accumulate in high part: probably same bound as for AVX already in flint?)
ulong _nmod32_vec_dot2_half_avx2(n32_srcptr v1, n32_srcptr v2, ulong len, nmod_t mod);
ulong _nmod32_vec_dot2_half_avx512(n32_srcptr v1, n32_srcptr v2, ulong len, nmod_t mod);

// TODO in progress: ifma attempt  (may fail for large len's)
#if HAVE_AVX512   // TODO handle AVX flags
ulong _nmod32_vec_dot_product_ifma256(n32_srcptr v1, n32_srcptr v2, ulong len, nmod_t mod);
ulong _nmod32_vec_dot_product_ifma512(n32_srcptr v1, n32_srcptr v2, ulong len, nmod_t mod);
#endif
#if HAVE_AVX_IFMA   // TODO handle AVX flags
ulong _nmod_vec_dot_product_avx_ifma(n32_srcptr v1, n32_srcptr v2, ulong len, nmod_t mod);
#endif


/** Several dot products with same left operand, as in vector-matrix product.
 *
 * . u has length at least len, len < 2^FLINT_BITS
 * . all entries of u  have <= max_bits_u bits <= FLINT_BITS
 * . v points to at least len vectors v[0],..,v[len-1] each of length
 * at least k, with entries of <= max_bits_v bits <= FLINT_BITS
 * . computes uv[j] = sum(u[i]*v[i][j], 0 <= i < len) modulo mod.n,
 * for 0 <= j < k
 * . does not assume input entries are reduced modulo mod.n
 */
//void _nmod_vec_dot_product_multi_2_v1_8(nn_ptr uv, nn_srcptr u, nn_srcptr * v,
//                                        ulong len, ulong k, nmod_t mod);

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD32_VEC__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
