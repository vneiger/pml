#ifndef __NMOD_EXTRA__H
#define __NMOD_EXTRA__H

/*------------------------------------------------------------*/
/*        TODO: safe to assume this exists?                   */
/*------------------------------------------------------------*/
#include <inttypes.h>

#include <flint/flint.h>
#include <flint/nmod_vec.h>

#if (defined __SIZEOF_INT128__ && GMP_LIMB_BITS == 64)
#define HAS_INT128
#define mp_dlimb_t unsigned __int128
#endif

/*------------------------------------------------------------*/
/* TODO: find the proper way to test                          */
/*------------------------------------------------------------*/
#ifdef __AVX2__
#define HAS_AVX2
#endif

#ifdef HAS_AVX2
#include <immintrin.h>
#endif

#if (GMP_LIMB_BITS == 64)
#define mp_hlimb_t uint32_t     // half-limb
#define mp_qlimb_t uint16_t     // quarter-limb
#define mp_hlimb_signed_t int32_t      // signed half-limb
#define mp_qlimb_signed_t int16_t       // signed quarter-limb
#endif


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*        A few extra functionalities for Fp                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* in-place reduces a mod p for a in [0,2p)                   */
/*------------------------------------------------------------*/
#define CORRECT_0_2P(a, p)                      \
    do {                                        \
        a -= (a >= p) ? p : 0;                  \
    } while(0)

/*------------------------------------------------------------*/
/* in-place reduces a mod p for a in (-p,p)                   */
/*------------------------------------------------------------*/
#define CORRECT_MINUSP_P(a, p)                          \
    do {                                                \
        a += ((mp_limb_signed_t) a < 0) ? (p) : 0;      \
    } while(0)

/*------------------------------------------------------------*/
/* returns the smallest i such that 2^i >= x                  */
/*------------------------------------------------------------*/
int next_power_of_two(mp_limb_t x);

/*------------------------------------------------------------*/
/* returns 1/p mod 2^k, assuming p is odd                     */
/* ill-defined when p is even                                 */
/*------------------------------------------------------------*/
mp_limb_t inverse_mod_power_of_two(mp_limb_t p, int k);

/*------------------------------------------------------------*/
/* finds an element of order at least n                       */
/* returns 0 if not found                                     */
/*------------------------------------------------------------*/
mp_limb_t nmod_find_root(long n, nmod_t mod);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 32 bit modular multiplication using a preconditionner      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* high word of a*b, assuming words are 32 bits               */
/*------------------------------------------------------------*/
static inline mp_hlimb_t mul_hi_32(mp_hlimb_t a, mp_hlimb_t b)
{
    return ((mp_limb_t) a * (mp_limb_t) b) >> 32;
}

/*------------------------------------------------------------*/
/* returns floor (2^32*b)/p                                   */
/*------------------------------------------------------------*/
static inline mp_hlimb_t prep_mul_mod_precon_32(mp_hlimb_t b, mp_hlimb_t p)
{
    return (((mp_limb_t) b) << 32) / ((mp_limb_t) p);
}

/*------------------------------------------------------------*/
/* preconditioned product                                     */
/* returns ab mod p                                           */
/* a is in [0..2^32), b is in [0..p)                          */
/*------------------------------------------------------------*/
static inline mp_hlimb_t mul_mod_precon_32(mp_hlimb_t a, mp_hlimb_t b, mp_hlimb_t p, mp_hlimb_t i_b)
{
    mp_hlimb_t q;
    mp_hlimb_signed_t t;
    q = mul_hi_32(i_b, a);
    t = a*b - q*p - p;
    if (t < 0)
        return t + p;
    else
        return t;
}

/*------------------------------------------------------------*/
/* returns ab mod p in [0..2p)                                */
/* a is in [0..2^32), b is in [0..p)                          */
/*------------------------------------------------------------*/
static inline mp_hlimb_t mul_mod_precon_unreduced_32(mp_hlimb_t a, mp_hlimb_t b, mp_hlimb_t p, mp_hlimb_t i_b)
{
    mp_hlimb_t q;
    q = mul_hi_32(i_b, a);
    return a*b - q*p;
}

#ifdef HAS_INT128

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 64 bit modular multiplication using a preconditionner      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* high word of a*b, assuming words are 64 bits               */
/*------------------------------------------------------------*/
static inline mp_limb_t mul_hi(mp_limb_t a, mp_limb_t b)
{
    mp_dlimb_t prod;
    prod = a * (mp_dlimb_t)b;
    return prod >> 64;
}

/*------------------------------------------------------------*/
/* returns floor (2^64*b)/p                                   */
/*------------------------------------------------------------*/
static inline mp_limb_t prep_mul_mod_precon(mp_limb_t b, mp_limb_t p)
{
    return ((mp_dlimb_t) b << 64) / p;
}

/*------------------------------------------------------------*/
/* preconditioned product                                     */
/* returns ab mod p                                           */
/* a is in [0..2^64), b is in [0..p)                          */
/*------------------------------------------------------------*/
static inline mp_limb_t mul_mod_precon(mp_limb_t a, mp_limb_t b, mp_limb_t p, mp_limb_t i_b)
{
    mp_limb_t q;
    long t;

    q = mul_hi(i_b, a);
    t = a*b - q*p - p;
    if (t < 0)
        return t + p;
    else
        return t;
}

/*------------------------------------------------------------*/
/* returns ab mod p in [0..2p)                                */
/* a is in [0..2^64), b is in [0..p)                          */
/*------------------------------------------------------------*/
static inline mp_limb_t mul_mod_precon_unreduced(mp_limb_t a, mp_limb_t b, mp_limb_t p, mp_limb_t i_b)
{
    mp_limb_t q;

    q = mul_hi(i_b, a);
    return a*b - q*p;
}


#endif


#ifdef HAS_AVX2

/*------------------------------------------------------------*/
/* high 32 bits of a*b                                        */
/*------------------------------------------------------------*/
static inline __m256i mm256_mulhi_epi32(__m256i a, __m256i b)
{
    __m256i prod02, prod13;
    prod02 = _mm256_mul_epu32(a, b);
    prod13 = _mm256_mul_epu32(_mm256_shuffle_epi32(a, 0xf5), _mm256_shuffle_epi32(b, 0xf5));
    return _mm256_unpackhi_epi64(_mm256_unpacklo_epi32(prod02, prod13), _mm256_unpackhi_epi32(prod02, prod13));
}

#endif


#endif
