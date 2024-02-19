#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

/* ------------------------------------------------------------ */
/* v1 and v2 have length at least len, len < 2^FLINT_BITS      */
/* all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/* all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/* computes sum(v1[i]*v2[i], 0 <= i < len)                      */
/* stores the result in 3 limbs of res                          */
/* ------------------------------------------------------------ */
void nmod_vec_integer_dot_product(mp_ptr res, mp_srcptr v1, mp_srcptr v2, ulong len, ulong max_bits1, ulong max_bits2)
{
    mp_limb_t s0, s1, u0, u1, t0, t1, t2;

    /* number of products we can do before overflow */
    const ulong log_nt = 2*FLINT_BITS - (max_bits1 + max_bits2);
    const ulong num_terms = (log_nt < FLINT_BITS) ? (UWORD(1) << log_nt) : (UWORD_MAX);

    t2 = UWORD(0);
    t1 = UWORD(0);
    t0 = UWORD(0);

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
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, u1, u0);
            v1 += num_terms;
            v2 += num_terms;
            len -= num_terms;
        }

    u1 = UWORD(0);
    u0 = UWORD(0);
    for (ulong i = 0; i < len; i++)
    {
        umul_ppmm(s1, s0, v1[i], v2[i]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }

    add_sssaaaaaa(res[2], res[1], res[0], t2, t1, t0, 0, u1, u0);

    return;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
