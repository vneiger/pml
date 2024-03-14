#include <flint/nmod.h>
#include "nmod_poly_fft.h"

/** fft dif evaluation, in place                                 
 * - len == 2**order, order >= 0
 * - input p has length at least len, seen as polynomial
 * - output p contains, in 0...len-1, the evaluations of trunc(p, len)
 *       at w**k, for 0 <= k < len, in bit reverse order
 * TODO : ASSUMES order suitable for w,                      
 *          and tables for 'order' already precomputed in F  
 *
 *  also assumes we can add a few terms before overflowing
 */
void _nmod_poly_dif_inplace_radix2_rec(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(p[0], p[1], F->mod);
    }
    else if (order == 2)
    {
        // for I == w**(2**(order-2)) square root of -1,
        //                            [1  1  1  1]
        // p <- [p[0] p[1] p[2] p[3]] [1 -1  I -I]
        //                            [1  1 -1 -1]
        //                            [1 -1 -I  I]
        // p[0] = (p[0] + p[2]) + (p[1] + p[3]);
        // p[1] = (p[0] + p[2]) - (p[1] + p[3]);
        // p[2] = (p[0] - p[2]) + I*(p[1] - p[3]);
        // p[3] = (p[0] - p[2]) - I*(p[1] - p[3]);
        ulong tmp02 = p[0] + p[2]; // p[0] + p[2], [0,2n-2]
        p[2] = p[0] + F->mod.n - p[2];    // p[0] - p[2], [1,2n-1]
        ulong tmp13 = p[1] + F->mod.n - p[3]; // p[1] - p[3], [1,2n-1]
        p[1] = p[1] + p[3];  // p[1] + p[3], [0,2n-2]

        p[0] = tmp02 + p[1];  // p[0] now ok,  [0, 4n-4]
        // TODO keep track of range and do separate functions with / without correcting excess
        if (p[0] >= 2*F->mod.n)
            p[0] -= 2*F->mod.n;
        if (p[0] >= F->mod.n)
            p[0] -= F->mod.n;
        p[1] = tmp02 + 2*F->mod.n - p[1];  // p[1] now ok,  [2,4n-2]
        if (p[1] >= 2*F->mod.n)
            p[1] -= 2*F->mod.n;
        if (p[1] >= F->mod.n)
            p[1] -= F->mod.n;

        // TODO compare also NMOD_MUL_PRENOM with n_mulmod2_preinv
        // (but maybe not for here: the former requires operands < modulus)
        // compute I*(p[1] - p[3]),  [0, n-1]
        tmp13 = n_mulmod2_preinv(F->tab_w[0][1], tmp13, F->mod.n, F->mod.ninv);
        p[3] = p[2] + F->mod.n - tmp13; // p[3] now ok, [2, 3n-1]
        if (p[3] >= 2*F->mod.n)
            p[3] -= 2*F->mod.n;
        else if (p[3] >= F->mod.n)
            p[3] -= F->mod.n;
        p[2] = p[2] + tmp13; // p[2] now ok, [1, 3n-2]
        if (p[2] >= 2*F->mod.n)
            p[2] -= 2*F->mod.n;
        else if (p[2] >= F->mod.n)
            p[2] -= F->mod.n;
    }
    else
    {
        const mp_ptr phi = p+(len/2);
        for (ulong k = 0; k < len/2; k++)
        {
            const ulong tmp = p[k];
            //p[k] = p[k] + phi[k];
            //phi[k] = tmp + F->mod.n - phi[k];
            p[k] = nmod_add(p[k], phi[k], F->mod);
            phi[k] = nmod_sub(tmp, phi[k], F->mod);
            phi[k] = n_mulmod2_preinv(phi[k], F->tab_w[order-2][k], F->mod.n, F->mod.ninv);
        }
        _nmod_poly_dif_inplace_radix2_rec(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec(phi, len/2, order-1, F);
    }
//    {
//        // apply radix-4, then call recursively
//        // for I == w**(2**(order-2)) square root of -1,
//        // block-wise with blocks of length len/4,
//        //                    [1  1  1  1]
//        // p <- [p0 p1 p2 p3] [1 -1  I -I]
//        //                    [1  1 -1 -1]
//        //                    [1 -1 -I  I]
//        // p0 = (p0 + p2) + (p1 + p3)
//        // p1[k] = F->tab_w[order-2][k] * ((p0[k] + p2[k]) - (p1[k] + p3[k]))
//        // p2 = (p0 - p2) + I*(p1 - p3) TODO
//        // p3 = (p0 - p2) - I*(p1 - p3) TODO
//        mp_ptr p0 = p;
//        mp_ptr p1 = p0 + (len/4);
//        mp_ptr p2 = p1 + (len/4);
//        mp_ptr p3 = p2 + (len/4);
//        for (ulong k = 0; k < len/4; k++)
//        {
//            ulong tmp02 = p0[k] + p2[k]; // p0 + p2, [0,2n-2]
//            p2[k] = p0[k] + F->mod.n - p2[k];    // p0 - p2, [1,2n-1]
//            ulong tmp13 = p1[k] + F->mod.n - p3[k]; // p1 - p3, [1,2n-1]
//            p1[k] = p1[k] + p3[k];  // p1 + p3, [0,2n-2]
//
//            p0[k] = tmp02 + p1[k];  // p0 now ok,  [0, 4n-4]
//                                  // TODO keep track of range and do separate functions with / without correcting excess
//            if (p0[k] >= 2*F->mod.n)
//                p0[k] -= 2*F->mod.n;
//            if (p0[k] >= F->mod.n)
//                p0[k] -= F->mod.n;
//            p1[k] = tmp02 + 2*F->mod.n - p1[k];  // p1 now ok,  [2,4n-2]
//            if (p1[k] >= 2*F->mod.n)
//                p1[k] -= 2*F->mod.n;
//            if (p1[k] >= F->mod.n)
//                p1[k] -= F->mod.n;
//
//            // TODO use precomputations on I
//            // TODO compare also NMOD_MUL_PRENOM with n_mulmod2_preinv
//            // (but maybe not for here: the former requires operands < modulus)
//            // compute I*(p1 - p3),  [0, n-1]
//            tmp13 = n_mulmod2_preinv(F->tab_w[0][1], tmp13, F->mod.n, F->mod.ninv);
//            p3[k] = p2[k] + F->mod.n - tmp13; // p3 now ok, [2, 3n-1]
//            if (p3[k] >= 2*F->mod.n)
//                p3[k] -= 2*F->mod.n;
//            else if (p3[k] >= F->mod.n)
//                p3[k] -= F->mod.n;
//            p2[k] = p2[k] + tmp13; // p2 now ok, [1, 3n-2]
//            if (p2[k] >= 2*F->mod.n)
//                p2[k] -= 2*F->mod.n;
//            else if (p2[k] >= F->mod.n)
//                p2[k] -= F->mod.n;
//        }
//    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
