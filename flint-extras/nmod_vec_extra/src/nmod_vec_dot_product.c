#include <flint/flint.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


/* ------------------------------------------------------------ */
/* number of limbs needed for a dot product of length len       */
/* all entries 1st vector have <= max_bits1 bits <= FLINT_BITS  */
/* all entries 2nd vector have <= max_bits2 bits <= FLINT_BITS  */
/* returns either 2 or 3                                        */
/* ------------------------------------------------------------ */
static inline 
ulong _nmod_vec_dot_bound_limbs_unbalanced(ulong len, ulong max_bits1, ulong max_bits2)
{
    mp_limb_t t2, t1, t0, u1, u0;

    const mp_limb_t a1 = (max_bits1 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits1) - 1;
    const mp_limb_t a2 = (max_bits2 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits2) - 1;
    
    umul_ppmm(t1, t0, a1, a2);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0)
        return 3;
    else
        return 2;
}


/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len <= 2^FLINT_BITS      */
/** all entries of v1 have <= max_bits1 bits <= FLINT_BITS       */
/** all entries of v2 have <= max_bits2 bits <= FLINT_BITS       */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 3 limbs to store the result before reduction            */
/*  ------------------------------------------------------------ */
static inline
mp_limb_t _nmod_vec_dot_product_3(mp_srcptr v1, mp_srcptr v2, ulong len, ulong max_bits1, ulong max_bits2, nmod_t mod)
{
    mp_limb_t s0, s1, u0, u1, t0, t1, t2, res;
    ulong i;
    /* number of products we can do before overflow */
    ulong num_terms;
    
    num_terms = 1 << (2*FLINT_BITS - (max_bits1 + max_bits2));
    i = 0;
    
    t2 = 0;
    t1 = 0;
    t0 = 0;

    if (num_terms >= 10)
        while (10 < len)
        {
            umul_ppmm(s1, s0, v1[0], v2[0]);
            add_ssaaaa(u1, u0, 0, 0, s1, s0);
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
            umul_ppmm(s1, s0, v1[8], v2[8]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            umul_ppmm(s1, s0, v1[9], v2[9]);
            add_ssaaaa(u1, u0, u1, u0, s1, s0);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, u1, u0);
            v1 += 10;
            v2 += 10;
            len -= 10;
        }
    else 
        while (num_terms < len)
        {
            u0 = UWORD(0);
            u1 = UWORD(0);
            for (i = 0; i < num_terms; i++)
            {
                umul_ppmm(s1, s0, v1[i], v2[i]);
                add_ssaaaa(u1, u0, u1, u0, s1, s0);
            }
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, u1, u0);
            v1 += num_terms;
            v2 += num_terms;
            len -= num_terms;
        }
    
    
    u0 = UWORD(0);
    u1 = UWORD(0);
    for (i = 0; i < len; i++)
    {
	umul_ppmm(s1, s0, v1[i], v2[i]);
	add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }

    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, u1, u0);
    NMOD_RED(t2, t2, mod);
    NMOD_RED3(res, t2, t1, t0, mod);
    return res;
}

/*  ------------------------------------------------------------ */
/** v1 and v2 have length at least len, len <= 2^FLINT_BITS      */
/** computes sum(v1[i]*v2[i], 0 <= i < len) modulo mod.n         */
/** uses 2 limbs to store the result before reduction            */
/*  ------------------------------------------------------------ */
static inline
mp_limb_t _nmod_vec_dot_product_2(mp_srcptr v1, mp_srcptr v2, ulong len, nmod_t mod)
{
    mp_limb_t s0, s1, u0, u1, res;
    ulong i;

    i = 0;
    u0 = 0;
    u1 = 0;
    
    while (10 < len)
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
        umul_ppmm(s1, s0, v1[8], v2[8]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        umul_ppmm(s1, s0, v1[9], v2[9]);
        add_ssaaaa(u1, u0, u1, u0, s1, s0);
        v1 += 10;
        v2 += 10;
        len -= 10;
    }

    for (i = 0; i < len; i++)
    {
	umul_ppmm(s1, s0, v1[i], v2[i]);
	add_ssaaaa(u1, u0, u1, u0, s1, s0);
    }
    
    NMOD2_RED2(res, u1, u0, mod);
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
    ulong n_limbs;

    n_limbs = _nmod_vec_dot_bound_limbs_unbalanced(len, max_bits1, max_bits2);

    if (n_limbs == 2)
        return _nmod_vec_dot_product_2(v1, v2, len, mod);
    else
        return _nmod_vec_dot_product_3(v1, v2, len, max_bits1, max_bits2, mod);
}
