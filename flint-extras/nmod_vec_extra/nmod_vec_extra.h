#ifndef __NMOD_VEC_EXTRA__H
#define __NMOD_VEC_EXTRA__H

#include "nmod_extra.h"

#ifdef __cplusplus
extern "C" {
#endif

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
        __m256i avx_x, avx_q, avx_cmp;
        /* slong j; */
        /* for (j = 0; j < 8; j++) */
        /* { */
        /*     y[j] = mul_mod_precon_32(x[j], a, mod.n, i_a); */
        /* } */

        avx_x = _mm256_load_si256((__m256i const*) x);
        avx_q = mm256_mulhi_epi32(avx_x, avx_i_a);
        avx_q = _mm256_sub_epi32(_mm256_mullo_epi32(avx_x, avx_a), _mm256_mullo_epi32(avx_q, avx_p));

        avx_cmp = _mm256_cmpgt_epi32(avx_q, avx_p_minus_1);
        avx_q = _mm256_sub_epi32(avx_q, _mm256_and_si256(avx_cmp, avx_p));
        
        _mm256_store_si256((__m256i*) y, avx_q);
        
        x += 8;
        y += 8;
    }
    
    nb = len - 8*nb;
    for (i = 0; i < nb; i++)
    {
        y[i] = mul_mod_precon_32(x[i], a, mod.n, i_a);
    }

}

#endif


#ifdef __cplusplus
}
#endif

#endif

