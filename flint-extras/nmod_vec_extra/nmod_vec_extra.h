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
#include "nmod_extra.h"

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




/*------------------------------------------------------------*/
/* y[i] = a*x[i] mod mod, i=0..len-1                          */
/*------------------------------------------------------------*/
static inline void _nmod_vec_scalar_mul(mp_ptr y, mp_srcptr x, slong len, mp_limb_t a, nmod_t mod)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        y[i] = nmod_mul(x[i], a, mod);
    }
}


#ifdef HAS_INT128

/*------------------------------------------------------------*/
/* y[i] = a*x[i] mod mod, i=0..len-1                          */
/* assumes int128 available                                   */
/* i_a must have been precomputed                             */
/*------------------------------------------------------------*/
static inline void _nmod_vec_int128_scalar_mul(mp_ptr y, mp_srcptr x, slong len,
                                               mp_limb_t a, mp_limb_t i_a, nmod_t mod)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        y[i] = mul_mod_precon(x[i], a, mod.n, i_a);
    }
}

#endif


#ifdef HAS_AVX2

/*------------------------------------------------------------*/
/* y[i] = a*x[i] mod mod, i=0..len-1                          */
/* works for < 32 bit primes                                  */
/* assumes avx2 available                                     */
/* i_a must have been precomputed                             */
/* x and y must be 32-byte aligned                            */
/*------------------------------------------------------------*/
static inline void _nmod_vec_avx2_32_scalar_mul(mp_hlimb_t * y, const mp_hlimb_t * x, slong len,
                                                mp_hlimb_t a, mp_hlimb_t i_a, nmod_t mod)
{
    slong i, nb;
    __m256i avx_a, avx_i_a, avx_p, avx_p_minus_1;

    nb = len/8;
    avx_a = _mm256_set1_epi32(a);
    avx_i_a = _mm256_set1_epi32(i_a);
    avx_p = _mm256_set1_epi32(mod.n);
    avx_p_minus_1 = _mm256_set1_epi32(mod.n - 1);

    for (i = 0; i < nb; i++)
    {
        __m256i avx_x;
        avx_x =  _mm256_load_si256((__m256i const*) x);
        _mm256_store_si256((__m256i*) y, mm256_mul_mod_precon(avx_x, avx_a, avx_p, avx_p_minus_1, avx_i_a));
        x += 8;
        y += 8;
    }

    nb = len - 8*nb;
    for (i = 0; i < nb; i++)
    {
        y[i] = mul_mod_precon_32(x[i], a, mod.n, i_a);
    }
}

/*------------------------------------------------------------*/
/* y[i] = a[i]*x[i] mod mod, i=0..len-1                       */
/* works for < 32 bit primes                                  */
/* assumes avx2 available                                     */
/* i_a must have been precomputed                             */
/* x and y must be 32-byte aligned                            */
/*------------------------------------------------------------*/
static inline void _nmod_vec_avx2_32_hadamard_mul(mp_hlimb_t * y, const mp_hlimb_t * x, const mp_hlimb_t * a,
                                                  const mp_hlimb_t * i_a, slong len, nmod_t mod)
{
    slong i, nb;
    __m256i avx_p, avx_p_minus_1;

    nb = len/8;
    avx_p = _mm256_set1_epi32(mod.n);
    avx_p_minus_1 = _mm256_set1_epi32(mod.n - 1);

    for (i = 0; i < nb; i++)
    {
        __m256i avx_x, avx_a, avx_i_a;
        avx_x =  _mm256_load_si256((__m256i const*) x);
        avx_a =  _mm256_load_si256((__m256i const*) a);
        avx_i_a =  _mm256_load_si256((__m256i const*) i_a);
        _mm256_store_si256((__m256i*) y, mm256_mul_mod_precon(avx_x, avx_a, avx_p, avx_p_minus_1, avx_i_a));
        a += 8;
        i_a += 8;
        x += 8;
        y += 8;
    }

    nb = len - 8*nb;
    for (i = 0; i < nb; i++)
    {
        y[i] = mul_mod_precon_32(x[i], a[i], mod.n, i_a[i]);
    }
}

#endif // HAS_AVX2


#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD_VEC_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
