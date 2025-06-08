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

#include "pml.h"
#include "machine_vectors.h"

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

/**********************************************************************
*                            DOT PRODUCT                             *
**********************************************************************/

// WARNING: everything below assumes both
//     modulus <= DOT2_ACC8_MAX_MODULUS  (about 2**30.5)
//     len <= DOT2_ACC8_MAX_LEN          (about 2**28.5)
// behaviour is undefined if that is not satisfied

uint _nmod32_vec_dot_split(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);

#if PML_HAVE_AVX2
uint _nmod32_vec_dot_split_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);

// duplicate msolve's matrix-vector product (for single dot prod)
// FIXME slower -> can be removed
uint _nmod32_vec_dot_msolve_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, uint PRIME);
#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
uint _nmod32_vec_dot_split_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);

// ifma variant
// TODO unfinished, requires len multiple of 16 or 32

/* the following was written with avx_ifma flag in mind, but
 * can probably be discarded */
uint _nmod32_vec_dot_ifma_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
uint _nmod32_vec_dot_ifma_avx512(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif  /* PML_HAVE_AVX512 */

/**********************************************************************
*                      2-DOT and 3-DOT PRODUCT                       *
**********************************************************************/

// similar to dot_split, computes res0 = vec1 * vec2_0 and res1 = vec1 * vec2_1
void _nmod32_vec_dot2_split(uint * res0, uint * res1,
                            n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                            slong len, nmod_t mod, ulong pow2_precomp);


#if PML_HAVE_AVX2
void _nmod32_vec_dot2_split_avx2(uint * res0, uint * res1,
                                 n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                 slong len, nmod_t mod, ulong pow2_precomp);

// similar to dot_split, computes res[i] = vec1 * vec2_i for i = 0, 1, 2
void _nmod32_vec_dot3_split_avx2(uint * res, n32_srcptr vec1,
                                 n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                 slong len, nmod_t mod, ulong pow2_precomp);
#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
void _nmod32_vec_dot2_split_avx512(uint * res0, uint * res1,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                   slong len, nmod_t mod, ulong pow2_precomp);

// similar to dot_split, computes res[i] = vec1 * vec2_i for i = 0, 1, 2, 3
void _nmod32_vec_dot3_split_avx512(uint * res,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                   slong len, nmod_t mod, ulong pow2_precomp);
#endif  /* PML_HAVE_AVX512 */

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

void _nmod32_vec_mdot2_split(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                             slong nrows, slong len, slong stride, nmod_t mod);

#if PML_HAVE_AVX2
void _nmod32_vec_mdot_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                 slong nrows, slong len, slong stride, nmod_t mod);

// like dot_split but handling 2 or 3 rows at a time
void _nmod32_vec_mdot2_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot3_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod);

// duplicate msolve's matrix-vector product
void _nmod32_vec_mdot_msolve_via_dot_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                          slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot_msolve_native_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                         slong nrows, slong len, slong stride, nmod_t mod);
#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
void _nmod32_vec_mdot_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                   slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot2_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod);
void _nmod32_vec_mdot3_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod);
#endif  /* PML_HAVE_AVX512 */

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD32_VEC__H
