#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


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
// TODO try other block sizes (16, 32)
static inline
mp_limb_t _nmod_vec_dot_product_1(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t res = UWORD(0);

    ulong i = 0;
    for (; i+7 < len; i += 8)
    {
        res +=   v1[i+0] * v2[i+0]
               + v1[i+1] * v2[i+1]
               + v1[i+2] * v2[i+2]
               + v1[i+3] * v2[i+3]
               + v1[i+4] * v2[i+4]
               + v1[i+5] * v2[i+5]
               + v1[i+6] * v2[i+6]
               + v1[i+7] * v2[i+7];
    }

    for (; i < len; i++)
        res += v1[i] * v2[i];

    NMOD_RED(res, res, mod);
    return res;
}

mp_limb_t _nmod_vec_dot_product_1_v0(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t res = UWORD(0);

    for (ulong i = 0; i < len; i++)
        res += v1[i] * v2[i];

    NMOD_RED(res, res, mod);
    return res;
}

mp_limb_t _nmod_vec_dot_product_1_vec8(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t vecres[8];
    mp_limb_t res = UWORD(0);

    ulong i = 0;
    for (; i+7 < len; i += 8)
    {
        vecres[0] = v1[i+0] * v2[i+0];
        vecres[1] = v1[i+1] * v2[i+1];
        vecres[2] = v1[i+2] * v2[i+2];
        vecres[3] = v1[i+3] * v2[i+3];
        vecres[4] = v1[i+4] * v2[i+4];
        vecres[5] = v1[i+5] * v2[i+5];
        vecres[6] = v1[i+6] * v2[i+6];
        vecres[7] = v1[i+7] * v2[i+7];
        res += vecres[0] + vecres[1] + vecres[2] + vecres[3] + vecres[4] + vecres[5] + vecres[6] + vecres[7];
    }

    for (; i < len; i++)
        res += v1[i] * v2[i];

    NMOD_RED(res, res, mod);
    return res;
}

mp_limb_t _nmod_vec_dot_product_1_vec32(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t vecres[32];
    mp_limb_t res = UWORD(0);

    ulong i = 0;
    for (; i+31 < len; i += 32)
    {
        vecres[0] = v1[i+0] * v2[i+0];
        vecres[1] = v1[i+1] * v2[i+1];
        vecres[2] = v1[i+2] * v2[i+2];
        vecres[3] = v1[i+3] * v2[i+3];
        vecres[4] = v1[i+4] * v2[i+4];
        vecres[5] = v1[i+5] * v2[i+5];
        vecres[6] = v1[i+6] * v2[i+6];
        vecres[7] = v1[i+7] * v2[i+7];
        vecres[8] = v1[i+8] * v2[i+8];
        vecres[9] = v1[i+9] * v2[i+9];
        vecres[10] = v1[i+10] * v2[i+10];
        vecres[11] = v1[i+11] * v2[i+11];
        vecres[12] = v1[i+12] * v2[i+12];
        vecres[13] = v1[i+13] * v2[i+13];
        vecres[14] = v1[i+14] * v2[i+14];
        vecres[15] = v1[i+15] * v2[i+15];
        vecres[16] = v1[i+16] * v2[i+16];
        vecres[17] = v1[i+17] * v2[i+17];
        vecres[18] = v1[i+18] * v2[i+18];
        vecres[19] = v1[i+19] * v2[i+19];
        vecres[20] = v1[i+20] * v2[i+20];
        vecres[21] = v1[i+21] * v2[i+21];
        vecres[22] = v1[i+22] * v2[i+22];
        vecres[23] = v1[i+23] * v2[i+23];
        vecres[24] = v1[i+24] * v2[i+24];
        vecres[25] = v1[i+25] * v2[i+25];
        vecres[26] = v1[i+26] * v2[i+26];
        vecres[27] = v1[i+27] * v2[i+27];
        vecres[28] = v1[i+28] * v2[i+28];
        vecres[29] = v1[i+29] * v2[i+29];
        vecres[30] = v1[i+30] * v2[i+30];
        vecres[31] = v1[i+31] * v2[i+31];
        res += vecres[0] + vecres[1] + vecres[2] + vecres[3] + vecres[4] + vecres[5] + vecres[6] + vecres[7] + vecres[8] + vecres[9] +
               vecres[10] + vecres[11] + vecres[12] + vecres[13] + vecres[14] + vecres[15] + vecres[16] + vecres[17] + vecres[18] + vecres[19] +
               vecres[20] + vecres[21] + vecres[22] + vecres[23] + vecres[24] + vecres[25] + vecres[26] + vecres[27] + vecres[28] + vecres[29] +
               vecres[30] + vecres[31];
    }

    for (; i < len; i++)
        res += v1[i] * v2[i];

    NMOD_RED(res, res, mod);
    return res;
}

mp_limb_t _nmod_vec_dot_product_1_v32(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t res = UWORD(0);

    ulong i = 0;
    for (; i+31 < len; i += 32)
    {
        res +=   v1[i+ 0] * v2[i+ 0]
               + v1[i+ 1] * v2[i+ 1]
               + v1[i+ 2] * v2[i+ 2]
               + v1[i+ 3] * v2[i+ 3]
               + v1[i+ 4] * v2[i+ 4]
               + v1[i+ 5] * v2[i+ 5]
               + v1[i+ 6] * v2[i+ 6]
               + v1[i+ 7] * v2[i+ 7]
               + v1[i+ 8] * v2[i+ 8]
               + v1[i+ 9] * v2[i+ 9]
               + v1[i+10] * v2[i+10]
               + v1[i+11] * v2[i+11]
               + v1[i+12] * v2[i+12]
               + v1[i+13] * v2[i+13]
               + v1[i+14] * v2[i+14]
               + v1[i+15] * v2[i+15]
               + v1[i+16] * v2[i+16]
               + v1[i+17] * v2[i+17]
               + v1[i+18] * v2[i+18]
               + v1[i+19] * v2[i+19]
               + v1[i+20] * v2[i+20]
               + v1[i+21] * v2[i+21]
               + v1[i+22] * v2[i+22]
               + v1[i+23] * v2[i+23]
               + v1[i+24] * v2[i+24]
               + v1[i+25] * v2[i+25]
               + v1[i+26] * v2[i+26]
               + v1[i+27] * v2[i+27]
               + v1[i+28] * v2[i+28]
               + v1[i+29] * v2[i+29]
               + v1[i+30] * v2[i+30]
               + v1[i+31] * v2[i+31];
    }

    for (; i < len; i++)
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


mp_limb_t _nmod_vec_dot_product_2_v16(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t s0, s1;
    mp_limb_t u0 = UWORD(0);
    mp_limb_t u1 = UWORD(0);

    ulong i = 0;
    for (; i+15 < len; i += 16)
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
        umul_ppmm(s1, s0, v1[i+8], v2[i+8]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+9], v2[i+9]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+10], v2[i+10]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+11], v2[i+11]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+12], v2[i+12]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+13], v2[i+13]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+14], v2[i+14]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[i+15], v2[i+15]);
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
