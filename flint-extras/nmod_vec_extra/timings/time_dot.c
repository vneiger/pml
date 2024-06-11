#include <stdlib.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/machine_vectors.h>
#include <flint/nmod_vec.h>
#include <flint/nmod.h>

// small values for testing before launching test:
//#define TIME_THRES 0.002
//#define NB_ITER 10
// full test:
#define TIME_THRES 0.2
#define NB_ITER 2500

// utility
static inline
void _nmod_vec_rand(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = n_randint(state, mod.n);
}

FLINT_FORCE_INLINE vec4n vec4n_mul(vec4n u, vec4n v)
{
    return _mm256_mul_epu32(u, v);
}

FLINT_FORCE_INLINE vec4n vec4n_zero()
{
    return _mm256_setzero_si256();
}

// new general dot macro
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

#define DOT_SP_NB 56
// new general dot macro
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
        for (ulong kxx = 0; kxx < len; kxx++)                                                                      \
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

//#if defined(__AVX2__)
//#define NMOD_VEC_DOT_PRODUCT NMOD_VEC_DOT_PRODUCT_AVX2
//#else
#define NMOD_VEC_DOT_PRODUCT NMOD_VEC_DOT_PRODUCT_SCALAR
//#endif

/*------------------------------------------------------------*/
/* timing of global func against current flint                */
/*------------------------------------------------------------*/

ulong time_vs_current_flint_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    const uint red_pow = (1L << DOT_SP_NB) % mod.n;

    nn_ptr v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v1[i], state, len, mod);
    }
    nn_ptr v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v2[i], state, len, mod);
    }

    { // TEST
        ulong res_new;
        NMOD_VEC_DOT_PRODUCT(res_new, v1[0], v2[0], len, mod, n_limbs, red_pow);
        ulong res_flint = _nmod_vec_dot(v1[0], v2[0], len, mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tt;
    long nb_iter;

    nn_srcptr v1i, v2i;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            NMOD_VEC_DOT_PRODUCT(buf, v1[i], v2[i], len, mod, n_limbs, red_pow);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            NMOD_VEC_DOT_PRODUCT(buf, v1[i], v2[i], len, mod, n_limbs, red_pow);
            res += buf;
        }
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;

    t2 = 0.0; nb_iter = 0;
    while (t2 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t2 /= nb_iter;

    printf("%.1e\t%.1e\t%.1e\t", t1, t2, t2/t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        _nmod_vec_clear(v1[i]);
        _nmod_vec_clear(v2[i]);
    }

    return res;
}

ulong time_vs_current_flint_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    const uint red_pow = (1L << DOT_SP_NB) % mod.n;

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);

    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    { // TEST
        ulong res_new;
        NMOD_VEC_DOT_PRODUCT(res_new, v1, v2, len, mod, n_limbs, red_pow);
        //res_new = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        ulong res_flint = _nmod_vec_dot(v1, v2, len, mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT_PRODUCT(buf, v1, v2, len, mod, n_limbs, red_pow);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT_PRODUCT(buf, v1, v2, len, mod, n_limbs, red_pow);
            res += buf;
        }
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;

    t2 = 0.0; nb_iter = 0;
    while (t2 < TIME_THRES)
    {
        ulong buf;
        slong j;
        nn_srcptr v1s = v1;
        nn_srcptr v2s = v2;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT(buf, j, len, v1s[j], v2s[j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT(buf, j, len, v1s[j], v2s[j], mod, n_limbs);
            res += buf;
        }
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t2 /= nb_iter;

    printf("%.1e\t%.1e\t%.1e\t", t1, t2, t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}


/*--------------------------------------------------------------*/
/* main                                                         */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125);

    const slong nlens = 12;
    const slong lens[] = {2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000};

    const slong nbits = 19;
    const slong bits[] = {17, 20, 23, 26, 29, 30, 31, 32, 33, 40, 50, 55, 57, 59, 60, 61, 62, 63, 64};

    const slong nfuns = 2;
    typedef ulong (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_vs_current_flint_cf,                // 0
        time_vs_current_flint_cu,                // 1
    };

    if (argc == 1)
    {
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];

            printf("bit/len");
            for (slong i = 0; i < nlens; i++)
                printf("\t%ld\t\t", lens[i]);
            printf("\n");

            for (slong j = 0; j < nbits; j++)
            {
                const slong b = bits[j];

                printf("%ldmin\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
                printf("\n");

                printf("%ldmid\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
                printf("\n");

                printf("%ldmax\t", b);
                if (b < 64)
                    for (slong i = 0; i < nlens; i++)
                        tfun(lens[i], (UWORD(1) << b) - 1, state);
                else
                    for (slong i = 0; i < nlens; i++)
                        tfun(lens[i], UWORD_MAX, state);
                printf("\n");
            }
        }
    }
    else if (argc == 2)  // function is given
    {
        const timefun tfun = funs[atoi(argv[1])];

        printf("bit/len");
        for (slong i = 0; i < nlens; i++)
            printf("\t%ld\t\t", lens[i]);
        printf("\n");

        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];

            printf("%ldmin\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
            printf("\n");

            printf("%ldmid\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
            printf("\n");

            printf("%ldmax\t", b);
            if (b < 64)
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << b) - 1, state);
            else
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], UWORD_MAX, state);
            printf("\n");
        }
    }
    else if (argc == 3)  // function + nbits given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);

        printf("bit/len");
        for (slong i = 0; i < nlens; i++)
            printf("\t%ld\t\t", lens[i]);
        printf("\n");

        //printf("%ldmin\t", b);
        //for (slong i = 0; i < nlens; i++)
        //    tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
        //printf("\n");

        printf("%ldmid\t", b);
        for (slong i = 0; i < nlens; i++)
            tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");

        //printf("%ldmax\t", b);
        //if (b < 64)
        //    for (slong i = 0; i < nlens; i++)
        //        tfun(lens[i], (UWORD(1) << b) - 1, state);
        //else
        //    for (slong i = 0; i < nlens; i++)
        //        tfun(lens[i], UWORD_MAX, state);
        printf("\n");
    }
    else if (argc == 4)  // function + nbits + len given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);
        const slong len = atoi(argv[3]);

        printf("bit/len");
        printf("\t%ld", len);
        printf("\n");

        //printf("%ldmin\t", b);
        //tfun(len, (UWORD(1) << (b-1)) + 1, state);
        //printf("\n");

        printf("%ldmid\t", b);
        tfun(len, (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");

        //printf("%ldmax\t", b);
        //if (b < 64)
        //    tfun(len, (UWORD(1) << b) - 1, state);
        //else
        //    tfun(len, UWORD_MAX, state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
