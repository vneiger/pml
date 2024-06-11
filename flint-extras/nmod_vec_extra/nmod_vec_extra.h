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

/*------------------------------------------------*/
/* v1 and v2 entries already reduced modulo mod.n */
/*------------------------------------------------*/
ulong nmod_vec_dot_product(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod, ulong n_limbs);

#define NMOD_VEC_DOT_PRODUCT_v1(res, i, len, expr1, expr2, mod, nlimbs) \
do                                                                   \
{                                                                    \
    res = UWORD(0);                                                  \
                                                                     \
    if (nlimbs == 1)                                                 \
    {                                                                \
        for (i = 0; i < len; i++)                                    \
            res += (expr1) * (expr2);                                \
                                                                     \
        NMOD_RED(res, res, mod);                                     \
    }                                                                \
                                                                     \
    else if (nlimbs == 2)                                            \
    {                                                                \
        ulong s0, s1;                                                \
        ulong u0 = UWORD(0);                                         \
        ulong u1 = UWORD(0);                                         \
                                                                     \
        for (i = 0; i+7 < len; )                                     \
        {                                                            \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
        }                                                            \
        for (; i < len; i++)                                         \
        {                                                            \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
        }                                                            \
                                                                     \
        NMOD2_RED2(res, u1, u0, mod);                                \
    }                                                                \
                                                                     \
    else if (nlimbs == 3)                                            \
    {                                                                \
        ulong s0, s1;                                                \
        ulong t2 = UWORD(0);                                         \
        ulong t1 = UWORD(0);                                         \
        ulong t0 = UWORD(0);                                         \
                                                                     \
        i = 0;                                                       \
        /* we can accumulate 8 terms if n == mod.n is such that */   \
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */   \
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */   \
        if (mod.n <= 6521908912666391107L)                           \
        {                                                            \
            slong u0, u1;                                            \
            for (; i+7 < len; )                                      \
            {                                                        \
                umul_ppmm(u1, u0, (expr1), (expr2));                 \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                add_sssaaaaaa(t2, t1, t0,                            \
                              t2, t1, t0,                            \
                              UWORD(0), u1, u0);                     \
            }                                                        \
                                                                     \
            u0 = UWORD(0);                                           \
            u1 = UWORD(0);                                           \
            for (; i < len; i++)                                     \
            {                                                        \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
            }                                                        \
                                                                     \
            add_sssaaaaaa(t2, t1, t0,                                \
                          t2, t1, t0,                                \
                          UWORD(0), u1, u0);                         \
        }                                                            \
        else                                                         \
        {                                                            \
            for (; i < len; i++)                                     \
            {                                                        \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, s1, s0);    \
            }                                                        \
        }                                                            \
                                                                     \
        NMOD_RED(t2, t2, mod);                                       \
        NMOD_RED3(res, t2, t1, t0, mod);                             \
    }                                                                \
} while(0);

// new general dot macro
#define DOT_SP_NB 56
#define NMOD_VEC_DOT_PRODUCT_SCALAR(res, v1, v2, len, mod, nlimbs, red_pow) \
do                                                                   \
{                                                                    \
    res = UWORD(0);                                                  \
                                                                     \
    if (nlimbs == 1)                                                 \
    {                                                                \
        ulong ixx = 0;                                               \
        for (; ixx < len; ixx++)                                     \
            res += (v1)[ixx] * (v2)[ixx];                            \
                                                                     \
        NMOD_RED(res, res, mod);                                     \
    }                                                                \
    else if (mod.n <= 1515531528 && len <= 134744072)                \
    {                                                                \
        const ulong low_bits = (1L << DOT_SP_NB) - 1;                \
        ulong dp_lo = 0;                                             \
        uint dp_hi = 0;                                              \
                                                                     \
        ulong kxx = 0;                                               \
        for (; kxx+7 < len; kxx += 8)                                \
        {                                                            \
            dp_lo += (v1)[kxx+0] * (v2)[kxx+0] +                     \
                     (v1)[kxx+1] * (v2)[kxx+1] +                     \
                     (v1)[kxx+2] * (v2)[kxx+2] +                     \
                     (v1)[kxx+3] * (v2)[kxx+3] +                     \
                     (v1)[kxx+4] * (v2)[kxx+4] +                     \
                     (v1)[kxx+5] * (v2)[kxx+5] +                     \
                     (v1)[kxx+6] * (v2)[kxx+6] +                     \
                     (v1)[kxx+7] * (v2)[kxx+7];                      \
                                                                     \
            dp_hi += dp_lo >> DOT_SP_NB;                             \
            dp_lo &= low_bits;                                       \
        }                                                            \
                                                                     \
        for (; kxx < len; kxx++)                                     \
            dp_lo += (v1)[kxx] * (v2)[kxx];                          \
                                                                     \
        NMOD_RED(res, ((ulong)red_pow * dp_hi) + dp_lo, mod);        \
    }                                                                \
    else if (nlimbs == 2)                                            \
    {                                                                \
        ulong s0, s1;                                                \
        ulong u0 = UWORD(0);                                         \
        ulong u1 = UWORD(0);                                         \
                                                                     \
        ulong kxx = 0;                                               \
        for ( ; kxx+7 < len; kxx += 8)                               \
        {                                                            \
            umul_ppmm(s1, s0, (v1)[kxx+0], (v2)[kxx+0]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+1], (v2)[kxx+1]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+2], (v2)[kxx+2]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+3], (v2)[kxx+3]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+4], (v2)[kxx+4]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+5], (v2)[kxx+5]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+6], (v2)[kxx+6]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+7], (v2)[kxx+7]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
        }                                                            \
        for (; kxx < len; kxx++)                                     \
        {                                                            \
            umul_ppmm(s1, s0, (v1)[kxx], (v2)[kxx]);                 \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
        }                                                            \
                                                                     \
        NMOD2_RED2(res, u1, u0, mod);                                \
    }                                                                \
                                                                     \
    else /* (nlimbs == 3) */                                         \
    {                                                                \
        ulong s0, s1;                                                \
        ulong t2 = UWORD(0);                                         \
        ulong t1 = UWORD(0);                                         \
        ulong t0 = UWORD(0);                                         \
                                                                     \
        ulong kxx = 0;                                               \
        /* we can accumulate 8 terms if n == mod.n is such that */   \
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */   \
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */   \
        if (mod.n <= 6521908912666391107L)                           \
        {                                                            \
            slong u0, u1;                                            \
            for ( ; kxx+7 < len; kxx += 8)                           \
            {                                                        \
                umul_ppmm(u1, u0, (v1)[kxx+0], (v2)[kxx+0]);         \
                umul_ppmm(s1, s0, (v1)[kxx+1], (v2)[kxx+1]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+2], (v2)[kxx+2]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+3], (v2)[kxx+3]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+4], (v2)[kxx+4]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+5], (v2)[kxx+5]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+6], (v2)[kxx+6]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+7], (v2)[kxx+7]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                add_sssaaaaaa(t2, t1, t0,                            \
                              t2, t1, t0,                            \
                              UWORD(0), u1, u0);                     \
            }                                                        \
                                                                     \
            u0 = UWORD(0);                                           \
            u1 = UWORD(0);                                           \
            for (; kxx < len; kxx++)                                 \
            {                                                        \
                umul_ppmm(s1, s0, (v1)[kxx], (v2)[kxx]);             \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
            }                                                        \
                                                                     \
            add_sssaaaaaa(t2, t1, t0,                                \
                          t2, t1, t0,                                \
                          UWORD(0), u1, u0);                         \
        }                                                            \
        else                                                         \
        {                                                            \
            for (ulong kxx = 0; kxx < len; kxx++)                    \
            {                                                        \
                umul_ppmm(s1, s0, (v1)[kxx], (v2)[kxx]);             \
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, s1, s0);    \
            }                                                        \
        }                                                            \
                                                                     \
        NMOD_RED(t2, t2, mod);                                       \
        NMOD_RED3(res, t2, t1, t0, mod);                             \
    }                                                                \
} while(0);

#define NMOD_VEC_DOT_PRODUCT_AVX2(res, v1, v2, len, mod, nlimbs, red_pow)                                          \
do                                                                                                                 \
{                                                                                                                  \
    res = UWORD(0);                                                                                                \
                                                                                                                   \
    if (nlimbs == 1)                                                                                               \
    {                                                                                                              \
        vec4n dp = vec4n_zero();                                                                                   \
        ulong kxx = 0;                                                                                             \
        for ( ; kxx+3 < len; kxx+=4)                                                                               \
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned((v1)+kxx), vec4n_load_unaligned((v2)+kxx)));         \
                                                                                                                   \
        res = dp[0] + dp[1] + dp[2] + dp[3];                                                                       \
                                                                                                                   \
        for (; kxx < len; kxx++)                                                                                   \
            res += (v1)[kxx] * (v2)[kxx];                                                                          \
                                                                                                                   \
        NMOD_RED(res, res, mod);                                                                                   \
    }                                                                                                              \
    else if (mod.n <= 1515531528 && len <= 134744072)                                                              \
    {                                                                                                              \
        const vec4n low_bits = vec4n_set_n((1L << DOT_SP_NB) - 1);                                                 \
        vec4n dp_lo = vec4n_zero();                                                                                \
        vec4n dp_hi = vec4n_zero();                                                                                \
                                                                                                                   \
        ulong kxx = 0;                                                                                             \
        for (; kxx+31 < len; kxx += 32)                                                                            \
        {                                                                                                          \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+ 0), vec4n_load_unaligned((v2)+kxx+ 0))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+ 4), vec4n_load_unaligned((v2)+kxx+ 4))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+ 8), vec4n_load_unaligned((v2)+kxx+ 8))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+12), vec4n_load_unaligned((v2)+kxx+12))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+16), vec4n_load_unaligned((v2)+kxx+16))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+20), vec4n_load_unaligned((v2)+kxx+20))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+24), vec4n_load_unaligned((v2)+kxx+24))); \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx+28), vec4n_load_unaligned((v2)+kxx+28))); \
                                                                                                                   \
            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SP_NB));                                     \
            dp_lo = vec4n_bit_and(dp_lo, low_bits);                                                                \
        }                                                                                                          \
                                                                                                                   \
        for (; kxx + 3 < len; kxx += 4)                                                                            \
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned((v1)+kxx), vec4n_load_unaligned((v2)+kxx)));   \
                                                                                                                   \
        dp_hi = vec4n_add(dp_hi, _mm256_srli_epi64(dp_lo, DOT_SP_NB));                                             \
        dp_lo = vec4n_bit_and(dp_lo, low_bits);                                                                    \
                                                                                                                   \
        ulong total_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3];                                                \
        const uint total_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (total_lo >> DOT_SP_NB);                 \
        total_lo &= (1L << DOT_SP_NB) - 1;                                                                         \
                                                                                                                   \
        for (; kxx < len; kxx++)                                                                                   \
            total_lo += (v1)[kxx] * (v2)[kxx];                                                                     \
                                                                                                                   \
        NMOD_RED(res, ((ulong)red_pow * total_hi) + total_lo, mod);                                                \
    }                                                                \
    else if (nlimbs == 2)                                            \
    {                                                                \
        ulong s0, s1;                                                \
        ulong u0 = UWORD(0);                                         \
        ulong u1 = UWORD(0);                                         \
                                                                     \
        ulong kxx = 0;                                               \
        for ( ; kxx+7 < len; kxx += 8)                               \
        {                                                            \
            umul_ppmm(s1, s0, (v1)[kxx+0], (v2)[kxx+0]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+1], (v2)[kxx+1]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+2], (v2)[kxx+2]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+3], (v2)[kxx+3]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+4], (v2)[kxx+4]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+5], (v2)[kxx+5]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+6], (v2)[kxx+6]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            umul_ppmm(s1, s0, (v1)[kxx+7], (v2)[kxx+7]);             \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
        }                                                            \
        for (; kxx < len; kxx++)                                     \
        {                                                            \
            umul_ppmm(s1, s0, (v1)[kxx], (v2)[kxx]);                 \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
        }                                                            \
                                                                     \
        NMOD2_RED2(res, u1, u0, mod);                                \
    }                                                                \
                                                                     \
    else /* (nlimbs == 3) */                                         \
    {                                                                \
        ulong s0, s1;                                                \
        ulong t2 = UWORD(0);                                         \
        ulong t1 = UWORD(0);                                         \
        ulong t0 = UWORD(0);                                         \
                                                                     \
        ulong kxx = 0;                                               \
        /* we can accumulate 8 terms if n == mod.n is such that */   \
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */   \
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */   \
        if (mod.n <= 6521908912666391107L)                           \
        {                                                            \
            slong u0, u1;                                            \
            for ( ; kxx+7 < len; kxx += 8)                           \
            {                                                        \
                umul_ppmm(u1, u0, (v1)[kxx+0], (v2)[kxx+0]);         \
                umul_ppmm(s1, s0, (v1)[kxx+1], (v2)[kxx+1]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+2], (v2)[kxx+2]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+3], (v2)[kxx+3]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+4], (v2)[kxx+4]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+5], (v2)[kxx+5]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+6], (v2)[kxx+6]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                umul_ppmm(s1, s0, (v1)[kxx+7], (v2)[kxx+7]);         \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                add_sssaaaaaa(t2, t1, t0,                            \
                              t2, t1, t0,                            \
                              UWORD(0), u1, u0);                     \
            }                                                        \
                                                                     \
            u0 = UWORD(0);                                           \
            u1 = UWORD(0);                                           \
            for (; kxx < len; kxx++)                                 \
            {                                                        \
                umul_ppmm(s1, s0, (v1)[kxx], (v2)[kxx]);             \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
            }                                                        \
                                                                     \
            add_sssaaaaaa(t2, t1, t0,                                \
                          t2, t1, t0,                                \
                          UWORD(0), u1, u0);                         \
        }                                                            \
        else                                                         \
        {                                                            \
            for (ulong kxx = 0; kxx < len; kxx++)                    \
            {                                                        \
                umul_ppmm(s1, s0, (v1)[kxx], (v2)[kxx]);             \
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, s1, s0);    \
            }                                                        \
        }                                                            \
                                                                     \
        NMOD_RED(t2, t2, mod);                                       \
        NMOD_RED3(res, t2, t1, t0, mod);                             \
    }                                                                \
} while(0);

#if defined(__AVX2__)
#define NMOD_VEC_DOT_PRODUCT NMOD_VEC_DOT_PRODUCT_AVX2
#else
#define NMOD_VEC_DOT_PRODUCT NMOD_VEC_DOT_PRODUCT_SCALAR
#endif

/*------------------------------------------------------------*/
/* small modulus: see implementation for constraints          */
/*------------------------------------------------------------*/

// splitting at 56 bits
ulong _nmod_vec_dot_mod32(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two);
ulong _nmod_vec_dot_mod32_avx2(nn_ptr a, nn_ptr b, ulong len, nmod_t mod, uint power_two);

FLINT_FORCE_INLINE
ulong nmod_vec_dot_mod32(nn_ptr a, nn_ptr b, ulong len, nmod_t mod)
{
    ulong power_two;
    NMOD_RED(power_two, 1L<<DOT_SP_NB, mod);
    return _nmod_vec_dot_mod32(a, b, len, mod, (uint)power_two);
}

FLINT_FORCE_INLINE
ulong nmod_vec_dot_mod32_avx2(nn_ptr a, nn_ptr b, ulong len, nmod_t mod)
{
    ulong power_two;
    NMOD_RED(power_two, 1L<<DOT_SP_NB, mod);
    return _nmod_vec_dot_mod32_avx2(a, b, len, mod, (uint)power_two);
}




/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** res[0] = dot(a1, b), res[1] = dot(a2, b)                  */
/** power_two = 2^45 mod p, p2 = (p,p), pinv2 = (1/p,1/p)     */
/*------------------------------------------------------------*/
void _nmod_vec_dot2_small_modulus(nn_ptr res,
                                  nn_ptr a1, nn_ptr a2, nn_ptr b, ulong len,
                                  ulong power_two,
                                  vec2d p2, vec2d pinv2);






/*------------------------------------------------------------*/
/* DRAFT / EXPERIMENTS                                        */
/*------------------------------------------------------------*/

ulong _nmod_vec_dot_product_1_avx2(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod);
ulong _nmod_vec_dot_product_1_avx512(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod);

// note: version split16 interesting on recent laptop (gcc does some vectorization)
// limited to nbits <= ~31 (bound to be better analyzed, numterms)
ulong _nmod_vec_dot_product_2_split16(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);
// note: version split26 interesting (beyond 30-31 bits) on recent laptop (gcc does some vectorization)
// limited to nbits <= ~52 (TODO bound to be better analyzed, numterms; potential fixes in code needed)
ulong _nmod_vec_dot_product_2_split26(nn_srcptr v1, nn_srcptr v2, ulong len, nmod_t mod);


















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
