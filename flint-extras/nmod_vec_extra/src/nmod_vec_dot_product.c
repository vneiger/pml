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


/**********************************************************************
*                         SINGLE DOT PRODUCT                         *
**********************************************************************/


/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len < 2^FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 1 limb to store the result before reduction             */
/*  ------------------------------------------------------------ */
static inline
mp_limb_t _nmod_vec_dot_product_1(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t res = UWORD(0);

    while (8 < len)
    {
        res +=   v1[0] * v2[0]
               + v1[1] * v2[1]
               + v1[2] * v2[2]
               + v1[3] * v2[3]
               + v1[4] * v2[4]
               + v1[5] * v2[5]
               + v1[6] * v2[6]
               + v1[7] * v2[7];
        v1 += 8;
        v2 += 8;
        len -= 8;
    }

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

    while (8 < len)
    {
        umul_ppmm(s1, s0, v1[0], v2[0]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[1], v2[1]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[2], v2[2]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[3], v2[3]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[4], v2[4]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[5], v2[5]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[6], v2[6]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[7], v2[7]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        v1 += 8;
        v2 += 8;
        len -= 8;
    }

    for (ulong i = 0; i < len; i++)
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

    if (num_terms >= 8)
        while (8 < len)
        {
            umul_ppmm(u1, u0, v1[0], v2[0]);
            umul_ppmm(s1, s0, v1[1], v2[1]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[2], v2[2]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[3], v2[3]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[4], v2[4]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[5], v2[5]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[6], v2[6]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[7], v2[7]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
            v1 += 8;
            v2 += 8;
            len -= 8;
        }
    else
        while (num_terms < len)
        {
            umul_ppmm(u1, u0, v1[0], v2[0]);
            for (ulong i = 1; i < num_terms; i++)
            {
                umul_ppmm(s1, s0, v1[i], v2[i]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
            }
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, UWORD(0), u1, u0);
            v1 += num_terms;
            v2 += num_terms;
            len -= num_terms;
        }


    u0 = UWORD(0);
    u1 = UWORD(0);
    for (ulong i = 0; i < len; i++)
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



/**********************************************************************
*                      MULTI DOT PRODUCT  (vec*mat)                  *
**********************************************************************/

static inline
void _nmod_vec_dot_product_multi_1(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    _nmod_vec_zero(uv, k);
    for (ulong i = 0; i < len; i++)
    {
        ulong j = 0;
        for (; j+8 < k; j += 8)
        {
            uv[j+0] += u[i] * v[i][j+0];
            uv[j+1] += u[i] * v[i][j+1];
            uv[j+2] += u[i] * v[i][j+2];
            uv[j+3] += u[i] * v[i][j+3];
            uv[j+4] += u[i] * v[i][j+4];
            uv[j+5] += u[i] * v[i][j+5];
            uv[j+6] += u[i] * v[i][j+6];
            uv[j+7] += u[i] * v[i][j+7];
        }
        for (; j < k; j++)
            uv[j] += u[i] * v[i][j];
    }

    for (ulong j = 0; j < k; j++)
        NMOD_RED(uv[j], uv[j], mod);
}

static inline
void _nmod_vec_dot_product_multi_2(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
}

static inline
void _nmod_vec_dot_product_multi_3(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k,
                                   ulong max_bits_u, ulong max_bits_v,
                                   nmod_t mod)
{
}





void nmod_vec_dot_product_multi(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                ulong len, ulong k,
                                ulong max_bits_u, ulong max_bits_v,
                                nmod_t mod)
{
    const ulong n_limbs = _nmod_vec_dot_bound_limbs_unbalanced(len, max_bits_u, max_bits_v);

    if (n_limbs == 2)
    {
        printf("NOT IMPLEMENTED YET\n");
        _nmod_vec_dot_product_multi_2(uv, u, v, len, k, mod);
        return;
    }
    if (n_limbs == 3)
    {
        printf("NOT IMPLEMENTED YET\n");
        _nmod_vec_dot_product_multi_3(uv, u, v, len, max_bits_u, max_bits_v, k, mod);
        return;
    }
    if (n_limbs == 1)
    {
        _nmod_vec_dot_product_multi_1(uv, u, v, len, k, mod);
        return;
    }
    // n_limbs == 0
    _nmod_vec_zero(uv, k);
    return;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
