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

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
