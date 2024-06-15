#ifndef DOT_H
#define DOT_H

#include <flint/nmod_types.h>

typedef enum
{
    _DOT0 = 0,
    _DOT1 = 1,
    _DOT2_32_SPLIT = 2,
    _DOT2_32 = 3,
    _DOT2 = 4,
    _DOT3_UNROLL = 5,
    _DOT3 = 6
} dot_method_t;

typedef struct
{
    dot_method_t method;
    ulong pow2_precomp;
} dot_params_t;

// compute dot parameters; returns nlimbs
dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod);

void nmod_mat_mul_flint(nmod_mat_t, const nmod_mat_t, const nmod_mat_t);
void nmod_mat_mul_newdot(nmod_mat_t, const nmod_mat_t, const nmod_mat_t);

#define DOT_SPLIT_BITS 56
#define DOT_SPLIT_MASK UWORD(72057594037927935) // (1L << DOT_SPLIT_BITS) - 1

// new general dot macro
#define _NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs, red_pow) \
do                                                                    \
{                                                                     \
    if (nlimbs == 1)                                                  \
    {                                                                 \
        res = UWORD(0);                                               \
        for (i = 0; i < (len); i++)                                   \
            res += (expr1) * (expr2);                                 \
        NMOD_RED(res, res, mod);                                      \
    }                                                                 \
    else if (mod.n <= UWORD(1515531528) && (len) <= WORD(134744072))  \
    {                                                                 \
        ulong dp_lo = 0;                                              \
        ulong dp_hi = 0;                                              \
                                                                      \
        for (i = 0; i+7 < (len); )                                    \
        {                                                             \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
                                                                      \
            dp_hi += dp_lo >> DOT_SPLIT_BITS;                         \
            dp_lo &= DOT_SPLIT_MASK;                                  \
        }                                                             \
                                                                      \
        for ( ; i < (len); i++)                                       \
            dp_lo += (expr1) * (expr2);                               \
                                                                      \
        res = red_pow * dp_hi + dp_lo;                                \
        NMOD_RED(res, res, mod);                                      \
    }                                                                 \
    else if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))                 \
    {                                                                 \
        ulong s0zz = UWORD(0);                                        \
        ulong s1zz = UWORD(0);                                        \
        for (i = 0; i < (len); i++)                                   \
        {                                                             \
            const ulong prodzz = (expr1) * (expr2);                   \
            add_ssaaaa(s1zz, s0zz, s1zz, s0zz, 0, prodzz);            \
        }                                                             \
        NMOD2_RED2(res, s1zz, s0zz, mod);                             \
    }                                                                 \
    else if (nlimbs == 2)                                             \
    {                                                                 \
        ulong u0zz = UWORD(0);                                        \
        ulong u1zz = UWORD(0);                                        \
                                                                      \
        for (i = 0; i+7 < (len); )                                    \
        {                                                             \
            ulong s0zz, s1zz;                                         \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
        }                                                             \
        for ( ; i < (len); i++)                                       \
        {                                                             \
            ulong s0zz, s1zz;                                         \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
        }                                                             \
                                                                      \
        NMOD2_RED2(res, u1zz, u0zz, mod);                             \
    }                                                                 \
    else if (nlimbs == 3)                                             \
    {                                                                 \
        ulong t2zz = UWORD(0);                                        \
        ulong t1zz = UWORD(0);                                        \
        ulong t0zz = UWORD(0);                                        \
                                                                      \
        /* we can accumulate 8 terms if n == mod.n is such that */    \
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */    \
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */    \
        if (mod.n <= 6521908912666391107L)                            \
        {                                                             \
            for (i = 0; i+7 < (len); )                                \
            {                                                         \
                ulong s0zz, s1zz;                                     \
                ulong u0zz = UWORD(0);                                \
                ulong u1zz = UWORD(0);                                \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                add_sssaaaaaa(t2zz, t1zz, t0zz,                       \
                              t2zz, t1zz, t0zz,                       \
                              UWORD(0), u1zz, u0zz);                  \
            }                                                         \
                                                                      \
            ulong s0zz, s1zz;                                         \
            ulong u0zz = UWORD(0);                                    \
            ulong u1zz = UWORD(0);                                    \
            for ( ; i < (len); i++)                                   \
            {                                                         \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
            }                                                         \
                                                                      \
            add_sssaaaaaa(t2zz, t1zz, t0zz,                           \
                          t2zz, t1zz, t0zz,                           \
                          UWORD(0), u1zz, u0zz);                      \
        }                                                             \
        else                                                          \
        {                                                             \
            for (i = 0; i < (len); i++)                               \
            {                                                         \
                ulong s0zz, s1zz;                                     \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_sssaaaaaa(t2zz, t1zz, t0zz,                       \
                              t2zz, t1zz, t0zz,                       \
                              UWORD(0), s1zz, s0zz);                  \
            }                                                         \
        }                                                             \
                                                                      \
        NMOD_RED(t2zz, t2zz, mod);                                    \
        NMOD_RED3(res, t2zz, t1zz, t0zz, mod);                        \
    }                                                                 \
    else   /* nlimbs == 0 */                                          \
    {                                                                 \
        res = UWORD(0);                                               \
    }                                                                 \
} while(0);














// new general dot macro, version with method
// unroll triple <-  three limbs and
        /* we can accumulate 8 terms if n == mod.n is such that */    
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */    
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */    
        // if (mod.n <= 6521908912666391107L)                            

// _DOT1   (single limb)
#define _NMOD_VEC_DOT1(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    res = UWORD(0);                                                   \
    for (i = 0; i < (len); i++)                                       \
        res += (expr1) * (expr2);                                     \
    NMOD_RED(res, res, mod);                                          \
} while(0);

// _DOT2_SPLIT   (two limbs, modulus < 2**32, splitting at 56 bits, 8-unrolling)
// (constraints in nmod_vec_dot_product_final.c)
#define _NMOD_VEC_DOT2_32_SPLIT(res, i, len, expr1, expr2, mod, red_pow) \
do                                                                       \
{                                                                        \
    ulong dp_lo = 0;                                                     \
    ulong dp_hi = 0;                                                     \
                                                                         \
    for (i = 0; i+7 < (len); )                                           \
    {                                                                    \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
                                                                         \
        dp_hi += dp_lo >> DOT_SPLIT_BITS;                                \
        dp_lo &= DOT_SPLIT_MASK;                                         \
    }                                                                    \
                                                                         \
    for ( ; i < (len); i++)                                              \
        dp_lo += (expr1) * (expr2);                                      \
                                                                         \
    res = red_pow * dp_hi + dp_lo;                                       \
    NMOD_RED(res, res, mod);                                             \
} while(0);

// _DOT2_SPLIT4   (two limbs, modulus < 2**32, splitting at 56 bits, 4-unrolling)
// (constraints in nmod_vec_dot_product_final.c)
#define _NMOD_VEC_DOT2_32_SPLIT4(res, i, len, expr1, expr2, mod, red_pow) \
do                                                                       \
{                                                                        \
    ulong dp_lo = 0;                                                     \
    ulong dp_hi = 0;                                                     \
                                                                         \
    for (i = 0; i+4 < (len); )                                           \
    {                                                                    \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
        dp_lo += (expr1) * (expr2); i++;                                 \
                                                                         \
        dp_hi += dp_lo >> DOT_SPLIT_BITS;                                \
        dp_lo &= DOT_SPLIT_MASK;                                         \
    }                                                                    \
                                                                         \
    for ( ; i < (len); i++)                                              \
        dp_lo += (expr1) * (expr2);                                      \
                                                                         \
    res = red_pow * dp_hi + dp_lo;                                       \
    NMOD_RED(res, res, mod);                                             \
} while(0);


// _DOT2_32   (two limbs, modulus < 2**32)
#define _NMOD_VEC_DOT2_32(res, i, len, expr1, expr2, mod)             \
do                                                                    \
{                                                                     \
    ulong s0zz = UWORD(0);                                            \
    ulong s1zz = UWORD(0);                                            \
    for (i = 0; i < (len); i++)                                       \
    {                                                                 \
        const ulong prodzz = (expr1) * (expr2);                       \
        add_ssaaaa(s1zz, s0zz, s1zz, s0zz, 0, prodzz);                \
    }                                                                 \
    NMOD2_RED2(res, s1zz, s0zz, mod);                                 \
} while(0);

// _DOT2   (two limbs, general)
#define _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    ulong u0zz = UWORD(0);                                            \
    ulong u1zz = UWORD(0);                                            \
                                                                      \
    for (i = 0; i+7 < (len); )                                        \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
    }                                                                 \
    for ( ; i < (len); i++)                                           \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
    }                                                                 \
                                                                      \
    NMOD2_RED2(res, u1zz, u0zz, mod);                                 \
} while(0);

// _DOT3   (three limbs, general)
#define _NMOD_VEC_DOT3(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    ulong t2zz = UWORD(0);                                            \
    ulong t1zz = UWORD(0);                                            \
    ulong t0zz = UWORD(0);                                            \
    for (i = 0; i < (len); i++)                                       \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_sssaaaaaa(t2zz, t1zz, t0zz,                               \
                        t2zz, t1zz, t0zz,                             \
                        UWORD(0), s1zz, s0zz);                        \
    }                                                                 \
                                                                      \
    NMOD_RED(t2zz, t2zz, mod);                                        \
    NMOD_RED3(res, t2zz, t1zz, t0zz, mod);                            \
} while(0);

// _DOT3_UNROLL   (three limbs, 8-unroll, 8*(mod.n - 1)**2 < 2**128)
#define _NMOD_VEC_DOT3_UNROLL(res, i, len, expr1, expr2, mod)         \
do                                                                    \
{                                                                     \
    ulong t2zz = UWORD(0);                                            \
    ulong t1zz = UWORD(0);                                            \
    ulong t0zz = UWORD(0);                                            \
                                                                      \
    for (i = 0; i+7 < (len); )                                        \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        ulong u0zz = UWORD(0);                                        \
        ulong u1zz = UWORD(0);                                        \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        add_sssaaaaaa(t2zz, t1zz, t0zz,                               \
                        t2zz, t1zz, t0zz,                             \
                        UWORD(0), u1zz, u0zz);                        \
    }                                                                 \
                                                                      \
    ulong s0zz, s1zz;                                                 \
    ulong u0zz = UWORD(0);                                            \
    ulong u1zz = UWORD(0);                                            \
    for ( ; i < (len); i++)                                           \
    {                                                                 \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
    }                                                                 \
                                                                      \
    add_sssaaaaaa(t2zz, t1zz, t0zz,                                   \
                    t2zz, t1zz, t0zz,                                 \
                    UWORD(0), u1zz, u0zz);                            \
                                                                      \
    NMOD_RED(t2zz, t2zz, mod);                                        \
    NMOD_RED3(res, t2zz, t1zz, t0zz, mod);                            \
} while(0);

#define _NMOD_VEC_DOT_NEW(res, i, len, expr1, expr2, mod, params) \
do                                                                         \
{                                                                          \
    res = UWORD(0);                                                        \
    if (params.method == _DOT1)                                                   \
        _NMOD_VEC_DOT1(res, i, len, expr1, expr2, mod)                     \
    else if (params.method == _DOT2_32_SPLIT)                                     \
        /* _NMOD_VEC_DOT2_32_SPLIT(res, i, len, expr1, expr2, mod, red_pow) */   \
        _NMOD_VEC_DOT2_32(res, i, len, expr1, expr2, mod)                  \
    else if (params.method == _DOT2_32)                                           \
        _NMOD_VEC_DOT2_32(res, i, len, expr1, expr2, mod)                  \
    else if (params.method == _DOT2)                                              \
        _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                     \
    else if (params.method == _DOT3_UNROLL)                                       \
        _NMOD_VEC_DOT3_UNROLL(res, i, len, expr1, expr2, mod)              \
    else if (params.method == _DOT3)                                              \
        _NMOD_VEC_DOT3(res, i, len, expr1, expr2, mod)                     \
} while(0);



#define FLINT_NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs) \
  do                                             \
  {                                              \
    ulong s0, s1, s2, t0, t1;                    \
    s0 = s1 = s2 = UWORD(0);                     \
    switch (nlimbs)                              \
    {                                            \
      case 1:                                    \
        for (i = 0; i < (len); i++)              \
          s0 += (expr1) * (expr2);               \
        NMOD_RED(s0, s0, mod);                   \
        break;                                   \
      case 2:                                    \
        if (mod.n <= (UWORD(1) << (FLINT_BITS / 2))) \
        {                                        \
          for (i = 0; i < (len); i++)            \
          {                                      \
            t0 = (expr1) * (expr2);              \
            add_ssaaaa(s1, s0, s1, s0, 0, t0);   \
          }                                      \
        }                                        \
        else if ((len) < 8)                      \
        {                                        \
          for (i = 0; i < len; i++)              \
          {                                      \
            umul_ppmm(t1, t0, (expr1), (expr2)); \
            add_ssaaaa(s1, s0, s1, s0, t1, t0);  \
          }                                      \
        }                                        \
        else                                     \
        {                                        \
          ulong v0, v1, u0, u1;                  \
          i = 0;                                 \
          if ((len) & 1)                         \
            umul_ppmm(v1, v0, (expr1), (expr2)); \
          else                                   \
            v0 = v1 = 0;                         \
          for (i = (len) & 1; i < (len); i++)    \
          {                                      \
            umul_ppmm(t1, t0, (expr1), (expr2)); \
            add_ssaaaa(s1, s0, s1, s0, t1, t0);  \
            i++;                                 \
            umul_ppmm(u1, u0, (expr1), (expr2)); \
            add_ssaaaa(v1, v0, v1, v0, u1, u0);  \
          }                                      \
          add_ssaaaa(s1, s0, s1, s0, v1, v0);    \
        }                                        \
        NMOD2_RED2(s0, s1, s0, mod);             \
        break;                                   \
      default:                                   \
        for (i = 0; i < (len); i++)              \
        {                                        \
          umul_ppmm(t1, t0, (expr1), (expr2));   \
          add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0); \
        }                                        \
        NMOD_RED(s2, s2, mod);                   \
        NMOD_RED3(s0, s2, s1, s0, mod);          \
        break;                                   \
    }                                            \
    res = s0;                                    \
  } while (0);





#endif
