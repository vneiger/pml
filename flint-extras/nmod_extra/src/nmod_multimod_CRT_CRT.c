#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/fft_small.h>
#include <flint/templates.h>
#include <flint/machine_vectors.h>
#include <flint/crt_helpers.h>

#include "nmod_extra.h"

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* helper functions for large modulus                           */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

/* ------------------------------------------------------------ */
/* out[i] = residues0[i] mod p, i < nb                          */
/* used for p >= 2^50                                           */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_1(mp_ptr out, 
                               mp_ptr residues0,
                               ulong nb,
                               nmod_t mod)
{
    ulong i;
    
    for (i = 0; i < nb; i++)
        NMOD_RED(out[i], residues0[i], mod);
}

/* ------------------------------------------------------------ */
/* out[i] = CRT(residues0[i],residues1[i]) mod p, i < nb        */
/* assumes p >= 2^50                                            */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_2(mp_ptr out,
                               mp_ptr residues0, mp_ptr residues1,
                               ulong nb,
                               mp_limb_t inv0, mp_limb_t inv1, // inv0 = 1/p1 mod p0, inv1 = 1/p0 mod p1
                               mp_ptr data,
                               nmod_t mod, nmod_t mod0, nmod_t mod1)
{
    ulong i;
    ulong t[2], r[2];

    for (i = 0; i < nb; i++)
    {
        _big_mul_2_1(r, t, data + 0*2, nmod_mul(inv0, residues0[i], mod0));
        _big_addmul_2_1(r, t, data + 1*2, nmod_mul(inv1, residues1[i], mod1));
        _reduce_big_sum_2(r, t, data + 2*2);
         NMOD2_RED2(out[i], r[1], r[0], mod);
    }
}

/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < 3) mod p, i < nb            */
/* assumes p >= 2^50                                            */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_3(mp_ptr out,
                               mp_ptr residues0, mp_ptr residues1, mp_ptr residues2,
                               ulong nb,
                               mp_limb_t inv0, mp_limb_t inv1, mp_limb_t inv2,
                               mp_ptr data,
                               nmod_t mod, nmod_t mod0, nmod_t mod1, nmod_t mod2)
{
    ulong i;
    ulong t[3], r[3];

    for (i = 0; i < nb; i++)
    {
        _big_mul_3_2(r, t, data + 0*3, nmod_mul(inv0, residues0[i], mod0));
        _big_addmul_3_2(r, t, data + 1*3, nmod_mul(inv1, residues1[i], mod1));
        _big_addmul_3_2(r, t, data + 2*3, nmod_mul(inv2, residues2[i], mod2));
        _reduce_big_sum_3(r, t, data + 3*3);
        NMOD_RED3(out[i], r[2], r[1], r[0], mod);
    }
}

/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < 4) mod p, i < nb            */
/* assumes p >= 2^50                                            */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_4(mp_ptr out,
                               mp_ptr residues0, mp_ptr residues1, mp_ptr residues2, mp_ptr residues3,
                               ulong nb,
                               mp_limb_t inv0, mp_limb_t inv1, mp_limb_t inv2, mp_limb_t inv3,
                               mp_ptr data,
                               nmod_t mod, nmod_t mod0, nmod_t mod1, nmod_t mod2, nmod_t mod3)
{
    ulong i;
    ulong t[4], r[4];

    for (i = 0; i < nb; i++)
    {
        _big_mul_4_3(r, t, data + 0*4, nmod_mul(inv0, residues0[i], mod0));
        _big_addmul_4_3(r, t, data + 1*4, nmod_mul(inv1, residues1[i], mod1));
        _big_addmul_4_3(r, t, data + 2*4, nmod_mul(inv2, residues2[i], mod2));
        _big_addmul_4_3(r, t, data + 3*4, nmod_mul(inv3, residues3[i], mod3));
        _reduce_big_sum_4(r, t, data + 4*4);
        NMOD_RED3(out[i], r[2], r[1], r[0], mod);
    }
}

/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < k) mod p, i < nb            */
/* k = 1,2,3,4                                                  */
/* used for p >= 2^50                                           */
/* ------------------------------------------------------------ */
static void nmod_large_modulus_CRT(mp_ptr out, mp_ptr *residues, ulong nb, nmod_multimod_CRT_t C)
{
    switch(C->num_primes)
    {
        case 1:
            _crt_1(out, residues[0], nb, C->mod);
            break;
        case 2:
            _crt_2(out, residues[0], residues[1], nb,
                   C->inverse_cofactors[0], C->inverse_cofactors[1],
                   C->data, C->mod, C->mod_primes[0], C->mod_primes[1]);
            break;
        case 3:
            _crt_3(out, residues[0], residues[1], residues[2], nb,
                   C->inverse_cofactors[0], C->inverse_cofactors[1], C->inverse_cofactors[2], 
                   C->data, C->mod, C->mod_primes[0], C->mod_primes[1], C->mod_primes[2]);
            break;
        case 4:
            _crt_4(out, residues[0], residues[1], residues[2], residues[3], nb,
                   C->inverse_cofactors[0], C->inverse_cofactors[1], C->inverse_cofactors[2], C->inverse_cofactors[3], 
                   C->data, C->mod, C->mod_primes[0], C->mod_primes[1], C->mod_primes[2], C->mod_primes[3]);
            break;
    }
}



/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* helper functions for small modulus                           */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

/* ------------------------------------------------------------ */
/* out[i] = residues0[i] mod p, i < nb                          */
/* used for p < 2^50                                            */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_1_small(mp_ptr out,
                                     mp_ptr residues0, 
                                     ulong nb,
                                     mp_limb_t p, double pinv)
{
    ulong i;

    i = 0;
    if (nb >= 4)
    {
        vec4d p4, pinv4;
        p4 = vec4d_set_d(p);
        pinv4 = vec4d_set_d(pinv);
        for (; i + 4 < nb; i += 4)
            vec4d_store_unaligned_mp_ptr(out + i,
                                         vec4d_reduce_to_0n(vec4d_load_unaligned_mp_ptr(residues0 + i), p4, pinv4));
    }
    for (; i < nb; i++)
        out[i] = (vec1n) vec1d_reduce_to_0n(residues0[i], p, pinv);
}


/* ------------------------------------------------------------ */
/* out[i] = CRT(residues0[i],residues1[i]) mod p, i < nb        */
/* assumes p < 2^50                                             */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_2_small(mp_ptr out,
                                     mp_ptr residues0, mp_ptr residues1,
                                     ulong nb,
                                     mp_limb_t p, double pinv, 
                                     mp_limb_t p1, double p1inv,
                                     double invp0_p1, // invp0_p1 = 1/p0 mod p1
                                     double p0_redp) // p0_redp = p0 mod p
{
    ulong i;
    vec1d alpha, alphap, alphapp;
    vec1n a, m0;
    
    i = 0;

    if (nb >= 4)
    {
        vec4d invp0_p14, p14, p1inv4, p4, pinv4, p0_redp4;
        invp0_p14 = vec4d_set_d(invp0_p1);
        p14 = vec4d_set_d(p1);
        p1inv4 = vec4d_set_d(p1inv);
        p4 = vec4d_set_d(p);
        pinv4 = vec4d_set_d(pinv);
        p0_redp4 = vec4d_set_d(p0_redp);

        for (; i + 4 < nb; i += 4)
        {
            vec4d alpha04, alpha4, alphap4, alphapp4, m04, a4;

            alpha04 = vec4d_load_unaligned_mp_ptr(residues0 + i);
            alpha4 = vec4d_sub(vec4d_load_unaligned_mp_ptr(residues1 + i), alpha04);
            alphap4 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha4, invp0_p14, p14, p1inv4), p14);   // alphap4 in [0..p1)
            alphapp4 = vec4d_reduce_to_0n(alphap4, p4, pinv4); // alphapp4 in [0..p)
            m04 = vec4d_reduce_to_0n(alpha04, p4, pinv4); // m04 in [0..p)
            a4 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alphapp4, p0_redp4, p4, pinv4), p4); // a in [0..p)

            vec4d_store_unaligned_mp_ptr(out + i, vec4d_addmod(m04, a4, p4));
        }
    }
        
    for (; i < nb; i++)
    {
        alpha = (vec1d) residues1[i] - (vec1d) residues0[i];   // (-p1..p1)
        /* mulmod needs product to be in (-2p1^2..2p1^2), so we are good */
        alphap = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha, invp0_p1, p1, p1inv), p1);   // alphap in [0..p1)
        alphapp = vec1d_reduce_to_0n(alphap, p, pinv); // alphapp in [0..p)
        m0 = (vec1n) vec1d_reduce_to_0n(residues0[i], p, pinv); // m0 in [0..p)
        /* mulmod needs product to be in (-2p^2..2p^2), result in (-p..p) */
        a = (vec1n) vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alphapp, p0_redp, p, pinv), p); // a in [0..p)
        out[i] = vec1n_addmod(m0, a, p);
    }
}

/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < 3) mod p, i < nb            */
/* assumes p < 2^50                                             */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_3_small(mp_ptr out,
                                     mp_ptr residues0, mp_ptr residues1, mp_ptr residues2,
                                     ulong nb,
                                     mp_limb_t p, double pinv, // pinv for reduction
                                     mp_limb_t p1, double p1inv, // p1inv for reduction
                                     mp_limb_t p2, double p2inv, // p2inv for reduction
                                     double p0_redpi, // p0 reduced modulo all others pi = p0 as a double
                                     double invp0_p1, // invp0_p1 = 1/p0 mod p1
                                     double p0_redp,  // p0 reduced modulo p
                                     double p0p1_red, double invp0p1_p2)
{
    ulong i;
    i = 0;

    // check: do we need p1 as mp_limb_t?
    if (nb >= 4)
    {
        vec4d invp0_p14, p14, p1inv4, p24, p2inv4, p4, pinv4, p0_redp4, p0_redpi4, invp0p1_p24, p0p1_redp4;
        invp0_p14 = vec4d_set_d(invp0_p1);
        p14 = vec4d_set_d((vec1d) p1);
        p1inv4 = vec4d_set_d(p1inv);
        p24 = vec4d_set_d((vec1d) p2);
        p2inv4 = vec4d_set_d(p2inv);
        p4 = vec4d_set_d(p);
        pinv4 = vec4d_set_d(pinv);
        p0_redp4 = vec4d_set_d(p0_redp);
        p0_redpi4 = vec4d_set_d(p0_redpi);
        invp0p1_p24 = vec4d_set_d(invp0p1_p2);
        p0p1_redp4 = vec4d_set_d(p0p1_red);

        for (; i + 4 < nb; i += 4)
        {
            vec4d alpha04, alpha14, alpha1_p0_red24, alpha24;

                
            /* find coefficients alpha0, alpha1, alpha2 s.t. c = alpha0 + p0 alpha1 + p0 p1 alpha2 */
            alpha04 = vec4d_load_unaligned_mp_ptr(residues0 + i); // [0..p0), so [0..p1)
            alpha14 = vec4d_sub(vec4d_load_unaligned_mp_ptr(residues1 + i), alpha04); // (-p1..p1)
            
            alpha14 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha14, invp0_p14, p14, p1inv4), p14); // [0..p1) so [0..p2)
            alpha1_p0_red24 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha14, p0_redpi4, p24, p2inv4), p24); // [0..p2)
            /* residues0[i] in [0..p0) so in [0..p2) */
            alpha1_p0_red24 = vec4d_add(alpha04, alpha1_p0_red24); // [0..2p2)
            alpha24 = vec4d_sub(vec4d_load_unaligned_mp_ptr(residues2 + i), alpha1_p0_red24); // (-2p2..p2)
            /* sufficient: alpha2 in (-2p2..2p2) */
            alpha24 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha24, invp0p1_p24, p24, p2inv4), p24); // [0..p2)
            
            /* reduce everything mod p */
            alpha04 = vec4d_reduce_to_0n(alpha04, p4, pinv4); // [0..p)
            alpha14 = vec4d_reduce_to_pm1no(alpha14, p4, pinv4); // (-p..p)
            alpha24 = vec4d_reduce_to_pm1no(alpha24, p4, pinv4); // (-p..p)

            vec4d_store_unaligned_mp_ptr(out + i,
                                         vec4d_addmod(alpha04,
                                                      vec4d_addmod(vec4d_reduce_pm1no_to_0n(
                                                                       vec4d_mulmod(alpha14, p0_redp4, p4, pinv4), p4),
                                                                   vec4d_reduce_pm1no_to_0n(
                                                                       vec4d_mulmod(alpha24, p0p1_redp4, p4, pinv4), p4),
                                                                   p4),
                                                      p4)
                );
        }

    }

    for (; i < nb; i++)
    {
        vec1d alpha0, alpha1, alpha1_p0_red2, alpha2;
        
        /* find coefficients alpha0, alpha1, alpha2 s.t. c = alpha0 + p0 alpha1 + p0 p1 alpha2 */
        alpha0 = (vec1d) residues0[i]; // [0..p0), so [0..p1)
        alpha1 = (vec1d) residues1[i] - alpha0; // (-p1..p1)

        alpha1 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, invp0_p1, (vec1d) p1, p1inv), p1); // [0..p1) so [0..p2)
        alpha1_p0_red2 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, p0_redpi, p2, p2inv), p2); // [0..p2)
        /* residues0[i] in [0..p0) so in [0..p2) */
        alpha1_p0_red2 = vec1d_add(alpha0, alpha1_p0_red2); // [0..2p2)
        alpha2 = vec1d_sub(residues2[i], alpha1_p0_red2); // (-2p2..p2)
        /* sufficient: alpha2 in (-2p2..2p2) */
        alpha2 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha2, invp0p1_p2, p2, p2inv), p2); // [0..p2)
                
        /* reduce everything mod p */
        alpha0 = vec1d_reduce_to_0n(alpha0, p, pinv); // [0..p)
        alpha1 = vec1d_reduce_to_pm1no(alpha1, p, pinv); // (-p..p)
        alpha2 = vec1d_reduce_to_pm1no(alpha2, p, pinv); // (-p..p)

        out[i] = (vec1n) vec1d_addmod(alpha0,
                                      vec1d_addmod(vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, p0_redp, p, pinv), p),
                                                   vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha2, p0p1_red, p, pinv), p),
                                                   p),
                                      p);
    }
}

/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < 4) mod p, i < nb            */
/* k = 1,2,3,4                                                  */
/* assumes p < 2^50                                             */
/* ------------------------------------------------------------ */
FLINT_FORCE_INLINE void _crt_4_small(mp_ptr out,
                                     mp_ptr residues0, mp_ptr residues1, mp_ptr residues2, mp_ptr residues3,
                                     ulong nb,
                                     mp_limb_t p, double pinv,
                                     mp_limb_t p1, double p1inv,
                                     mp_limb_t p2, double p2inv,
                                     mp_limb_t p3, double p3inv,
                                     double p0_redpi,
                                     double invp0_p1, 
                                     double p0_redp, double p0p1_red, double p0p1p2_red,
                                     double invp0p1_p2,
                                     double p0p1_red3, double invp0p1p2_p3)
{
    ulong i;
    i = 0;
    
    /* check: do we need p1 as mp_limb_t? */
    if (nb >= 4)
    {
        vec4d invp0_p14, p14, p1inv4, p24, p2inv4, p4, p34, p3inv4, pinv4, p0_redp4, p0_redpi4, invp0p1_p24, p0p1_redp4, p0p1_red34, invp0p1p2_p34, p0p1p2_red4;
        invp0_p14 = vec4d_set_d(invp0_p1);
        p14 = vec4d_set_d((vec1d) p1);
        p1inv4 = vec4d_set_d(p1inv);
        p24 = vec4d_set_d((vec1d) p2);
        p2inv4 = vec4d_set_d(p2inv);
        p34 = vec4d_set_d((vec1d) p3);
        p3inv4 = vec4d_set_d(p3inv);
        p4 = vec4d_set_d(p);
        pinv4 = vec4d_set_d(pinv);
        p0_redp4 = vec4d_set_d(p0_redp);
        p0_redpi4 = vec4d_set_d(p0_redpi);
        invp0p1_p24 = vec4d_set_d(invp0p1_p2);
        p0p1_redp4 = vec4d_set_d(p0p1_red);
        p0p1_red34 = vec4d_set_d(p0p1_red3);
        invp0p1p2_p34 = vec4d_set_d(invp0p1p2_p3);
        p0p1p2_red4 = vec4d_set_d(p0p1p2_red);
        
        for (; i + 4 < nb; i += 4)
        {
            vec4d alpha04, alpha14, alpha1_p0_red24, alpha24, alpha34, alpha1_p0_red34, alpha2_p0p1_red34;
            alpha04 = vec4d_load_unaligned_mp_ptr(residues0 + i); // [0..p0), so [0..p1)
            alpha14 = vec4d_sub(vec4d_load_unaligned_mp_ptr(residues1 + i), alpha04);   // (-p1..p1)
            
            alpha14 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha14, invp0_p14, p14, p1inv4), p14);  // [0..p1) so [0..p2)
            alpha1_p0_red24 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha14, p0_redpi4, p24, p2inv4), p24); // [0..p2)
            alpha1_p0_red24 = vec4d_add(alpha04, alpha1_p0_red24); // [0..2p2)
            alpha24 = vec4d_sub(vec4d_load_unaligned_mp_ptr(residues2 + i), alpha1_p0_red24); // (-2p2..p2)
            alpha24 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha24, invp0p1_p24, p24, p2inv4), p24); // [0..p2)
            
            /* residues0[i] in [0..p0) so in [0..p3) */
            alpha1_p0_red34 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha14, p0_redpi4, p34, p3inv4), p34); // [0..p3)
            alpha2_p0p1_red34 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha24, p0p1_red34, p34, p3inv4), p34); // [0..p3)
            alpha34 = vec4d_add(vec4d_addmod(alpha04, alpha1_p0_red34, p34), alpha2_p0p1_red34); // inner sum [0..p3], result [0..2p3)

            alpha34 = vec4d_sub(vec4d_load_unaligned_mp_ptr(residues3 + i), alpha34); // (-2p3..p3)
            alpha34 = vec4d_reduce_pm1no_to_0n(vec4d_mulmod(alpha34, invp0p1p2_p34, p34, p3inv4), p34); // [0..p3)
            
            /* reduce everything mod p */
            alpha04 = vec4d_reduce_to_0n(alpha04, p4, pinv4); // [0..p)
            alpha14 = vec4d_reduce_to_pm1no(alpha14, p4, pinv4); // (-p..p)
            alpha24 = vec4d_reduce_to_pm1no(alpha24, p4, pinv4); // (-p..p)
            alpha34 = vec4d_reduce_to_pm1no(alpha34, p4, pinv4); // (-p..p)

            vec4d_store_unaligned_mp_ptr(out + i,
                                         vec4d_addmod(alpha04,
                                                      vec4d_addmod(
                                                          vec4d_addmod(
                                                              vec4d_reduce_pm1no_to_0n(
                                                                  vec4d_mulmod(alpha14, p0_redp4, p4, pinv4), p4),
                                                              vec4d_reduce_pm1no_to_0n(
                                                                  vec4d_mulmod(alpha24, p0p1_redp4, p4, pinv4), p4),
                                                              p4),
                                                          vec4d_reduce_pm1no_to_0n(
                                                              vec4d_mulmod(alpha34, p0p1p2_red4, p4, pinv4), p4),
                                                          p4),
                                                      p4));

        }
    }
    
    for (; i < nb; i++)
    {
        vec1d alpha0, alpha1, alpha1_p0_red2, alpha2, alpha3, alpha1_p0_red3, alpha2_p0p1_red3;
        
        /* find coefficients alpha0, alpha1, alpha2, alpha4 s.t. c = alpha0 + p0 alpha1 + p0 p1 alpha2 + p0 p1 p2 alpha3 */
        // alpha0 = c0
        // alpha1 = (c1 - alpha0)/p0 mod p1
        // alpha2 = (c2 - (alpha0 + alpha1*p0))/p0p1 mod p2
        // alpha3 = (c3 - (alpha0 + alpha1*p0 + alpha2*p0*p1))/p0p1p2 mod p3

        alpha0 = (vec1d) residues0[i]; // [0..p0), so [0..p1)
        alpha1 = (vec1d) residues1[i] - alpha0;   // (-p1..p1)

        alpha1 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, invp0_p1, (vec1d) p1, p1inv), p1);   // [0..p1) so [0..p2)
        alpha1_p0_red2 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, p0_redpi, p2, p2inv), p2); // [0..p2)
        // residues0[i] in [0..p0) so in [0..p2)
        alpha1_p0_red2 = vec1d_add(alpha0, alpha1_p0_red2); // [0..2p2)
        alpha2 = vec1d_sub(residues2[i], alpha1_p0_red2); // (-2p2..p2)
        // sufficient: alpha2 in (-2p2..2p2)
        alpha2 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha2, invp0p1_p2, p2, p2inv), p2); // [0..p2)

        // residues0[i] in [0..p0) so in [0..p3)
        alpha1_p0_red3 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, p0_redpi, p3, p3inv), p3); // [0..p3)
        alpha2_p0p1_red3 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha2, p0p1_red3, p3, p3inv), p3); // [0..p3)
        alpha3 = vec1d_add(vec1d_addmod(alpha0, alpha1_p0_red3, p3), alpha2_p0p1_red3); // inner sum [0..p3], result [0..2p3)

        alpha3 = vec1d_sub((vec1d) residues3[i], alpha3); // (-2p3..p3)
        alpha3 = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha3, invp0p1p2_p3, p3, p3inv), p3); // [0..p3)
        
        // reduce everything mod p
        alpha0 = vec1d_reduce_to_0n(alpha0, p, pinv); // [0..p)
        alpha1 = vec1d_reduce_to_pm1no(alpha1, p, pinv); // (-p..p)
        alpha2 = vec1d_reduce_to_pm1no(alpha2, p, pinv); // (-p..p)
        alpha3 = vec1d_reduce_to_pm1no(alpha3, p, pinv); // (-p..p)

        out[i] = (vec1n) vec1d_addmod(alpha0,
                                      vec1d_addmod(
                                          vec1d_addmod(
                                              vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha1, p0_redp, p, pinv), p),
                                              vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha2, p0p1_red, p, pinv), p),
                                              p),
                                          vec1d_reduce_pm1no_to_0n(vec1d_mulmod(alpha3, p0p1p2_red, p, pinv), p),
                                          p),
                                      p);
    }
}


/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < k) mod p, i < nb            */
/* k = 1,2,3,4                                                  */
/* assumes p < 2^50                                             */
/* ------------------------------------------------------------ */
static void nmod_small_modulus_CRT(mp_ptr out, mp_ptr *residues, ulong nb, nmod_multimod_CRT_t C)
{
    switch(C->num_primes)
    {
        case 1:
            _crt_1_small(out, residues[0], nb, C->p, C->pinv);
            break;
        case 2:
            _crt_2_small(out, residues[0], residues[1], nb, C->p, C->pinv, C->primes[1], C->primes_inv[1],
                         C->invp0_p1, C->p0_red);
            break;
        case 3:
            _crt_3_small(out, residues[0], residues[1], residues[2], nb,
                         C->p, C->pinv, C->primes[1], C->primes_inv[1], C->primes[2], C->primes_inv[2],
                         (vec1d) C->primes[0],
                         C->invp0_p1, C->p0_red, C->p0p1_red,  C->invp0p1_p2);
            break;
        case 4:
            _crt_4_small(out, residues[0], residues[1], residues[2], residues[3], nb,
                         C->p, C->pinv,
                         C->primes[1], C->primes_inv[1], C->primes[2], C->primes_inv[2],
                         C->primes[3], C->primes_inv[3],
                         (vec1d) C->primes[0], C->invp0_p1,
                         C->p0_red, C->p0p1_red, C->p0p1p2_red,
                         C->invp0p1_p2,
                         C->p0p1_red3, C->invp0p1p2_p3);
            break;
        default:
            break;
    }
}


/* ------------------------------------------------------------ */
/* out[i] = CRT(residues[j][i], j < k) mod p, i < nb            */
/* k = 1,2,3,4                                                  */
/* ------------------------------------------------------------ */
void nmod_multimod_CRT_CRT(mp_ptr out, mp_ptr *residues, ulong nb, nmod_multimod_CRT_t C)
{
    if (C->p < (1L << 50)) // small modulus: use SIMD floating-point representation 
        nmod_small_modulus_CRT(out, residues, nb, C);
    else
        nmod_large_modulus_CRT(out, residues, nb, C);
}
