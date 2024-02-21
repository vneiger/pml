//#include <flint/longlong.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

#define nbits_adhoc 32
#define __ll_B_adhoc (1 << (nbits_adhoc / 2))
#define __ll_lowhi_parts(tlo,thi,t)       \
      tlo = (__uint32_t) t;               \
      thi = ((tlo) >> (nbits_adhoc / 2)); \
      tlo = (tlo) & (__ll_B_adhoc - 1);

/* ------------------------------------------------------------ */
/* number of limbs needed for a dot product of length len       */
/* all entries 1st vector have <= max_bits1 bits <= FLINT_BITS  */
/* all entries 2nd vector have <= max_bits2 bits <= FLINT_BITS  */
/* returns 0, 1, 2, 3                                           */
/* ------------------------------------------------------------ */
static inline
ulong _nmod_vec_dot_bound_limbs_unbalanced(ulong len, ulong max_bits1, ulong max_bits2)
{
    const mp_limb_t a1 = (max_bits1 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits1) - 1;
    const mp_limb_t a2 = (max_bits2 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits2) - 1;

    mp_limb_t t2, t1, t0, u1, u0;
    umul_ppmm(t1, t0, a1, a2);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0)
        return 3;
    if (t1 != 0)
        return 2;
    return (t0 != 0);
}



/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 1 limb to store the result before reduction             */
/*  ------------------------------------------------------------ */
// Feb 2024: tried other block sizes (8, 16, 32), and also
// storing intermediate 8/16/32 products in vector before summing;
// it seems this simplest variant remains the most efficient (gcc is better
// at vectorizing it, both on slightly older machine and recent avx512 machine)
static inline
mp_limb_t _nmod_vec_dot_product_1(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t res = UWORD(0);

    for (ulong i = 0; i < len; i++)
        res += v1[i] * v2[i];

    NMOD_RED(res, res, mod);
    return res;
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 2 limbs to store the result before reduction            */
/*  ------------------------------------------------------------ */
static inline
mp_limb_t _nmod_vec_dot_product_2(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t s0, s1;
    mp_limb_t u0 = UWORD(0);
    mp_limb_t u1 = UWORD(0);

    ulong i = 0;
    for (; i+7 < len; i += 8)
    {
        umul_ppmm(s1, s0, v1[i+0], v2[i+0]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+1], v2[i+1]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+2], v2[i+2]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+3], v2[i+3]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+4], v2[i+4]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+5], v2[i+5]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+6], v2[i+6]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+7], v2[i+7]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }

    for (; i < len; i++)
    {
        umul_ppmm(s1, s0, v1[i], v2[i]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }

    mp_limb_t res;
    NMOD2_RED2(res, u1, u0, mod);
    return res;
}

mp_limb_t _nmod_vec_dot_product_2_vuint32(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    __uint32_t v1hi, v1lo, v2hi, v2lo;
    ulong ulo = UWORD(0);
    ulong umi = UWORD(0);
    ulong uhi = UWORD(0);
    for (ulong i = 0; i < len; i++)
    {
        __ll_lowhi_parts(v1lo, v1hi, v1[i]);
        __ll_lowhi_parts(v2lo, v2hi, v2[i]);
        ulo += v1lo * v2lo;
        umi += v1lo * v2hi + v1hi * v2lo;
        uhi += v1hi * v2hi;
    }

    // result:
    // ulo + 2**x umi + 2**(2x) uhi,  x = nbits_adhoc / 2
    // umi : split from 0...63-x and 64-x to 64
    // uhi : split from 0...63-2x and 63-2x to 64
    ulong of = (umi >> (64 - nbits_adhoc/2)) + (uhi >> (64 - nbits_adhoc));
    ulong uf = ((umi & ((1UL << (64 - nbits_adhoc/2)) - 1)) << (nbits_adhoc/2))
             + ((uhi & ((1UL << (64 - nbits_adhoc)) - 1)) << (nbits_adhoc));
    mp_limb_t res;
    NMOD2_RED2(res, of, uf + ulo, mod);
    return res;
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 3 limbs to store the result before reduction            */
/*  ------------------------------------------------------------ */
static inline
mp_limb_t _nmod_vec_dot_product_3(mp_srcptr v1, mp_srcptr v2, ulong len, ulong max_bits1, ulong max_bits2, nmod_t mod)
{
    /* number of products we can do before overflow */
    const ulong log_nt = 2*FLINT_BITS - (max_bits1 + max_bits2);
    const ulong num_terms = (log_nt < FLINT_BITS) ? (UWORD(1) << log_nt) : (UWORD_MAX);

    mp_limb_t s0, s1, u0, u1;
    mp_limb_t t2 = UWORD(0);
    mp_limb_t t1 = UWORD(0);
    mp_limb_t t0 = UWORD(0);

    ulong i = 0;
    if (num_terms >= 8)
        for (; i+7 < len; i += 8)
        {
            umul_ppmm(u1, u0, v1[i+0], v2[i+0]);
            umul_ppmm(s1, s0, v1[i+1], v2[i+1]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+2], v2[i+2]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+3], v2[i+3]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+4], v2[i+4]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+5], v2[i+5]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+6], v2[i+6]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[i+7], v2[i+7]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
        }
    else
        for (; i+num_terms < len; i += num_terms)
        {
            umul_ppmm(u1, u0, v1[i], v2[i]);
            for (ulong ii = 1; ii < num_terms; ii++)
            {
                umul_ppmm(s1, s0, v1[i+ii], v2[i+ii]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
            }
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
        }

    u0 = UWORD(0);
    u1 = UWORD(0);
    for (; i < len; i++)
    {
        umul_ppmm(s1, s0, v1[i], v2[i]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }

    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
    NMOD_RED(t2, t2, mod);

    mp_limb_t res;
    NMOD_RED3(res, t2, t1, t0, mod);
    return res;
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len <= 2^FLINT_BITS      */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** does not assume input is reduced modulo mod.n                */
/*  ------------------------------------------------------------ */
mp_limb_t nmod_vec_dot_product(mp_srcptr v1, mp_srcptr v2, ulong len, ulong max_bits1, ulong max_bits2, nmod_t mod)
{
    const ulong n_limbs = _nmod_vec_dot_bound_limbs_unbalanced(len, max_bits1, max_bits2);

    if (n_limbs == 2)
        return _nmod_vec_dot_product_2(v1, v2, len, mod);
    if (n_limbs == 3)
        return _nmod_vec_dot_product_3(v1, v2, len, max_bits1, max_bits2, mod);
    if (n_limbs == 1)
        return _nmod_vec_dot_product_1(v1, v2, len, mod);
    return UWORD(0);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
