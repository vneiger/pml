#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include "nmod_poly_fft.h"

/***************************
*  DIF radix 2 recursive  *
***************************/

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
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        mp_limb_t p7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else
    {
        const mp_ptr phi = p+(len/2);
        for (ulong k = 0; k < len/2; k++)
        {
            const ulong tmp = p[k];
            p[k] = nmod_add(p[k], phi[k], F->mod);
            phi[k] = nmod_sub(tmp, phi[k], F->mod);
            NMOD_MUL_PRENORM(phi[k], phi[k] << F->mod.norm, F->tab_w[order-2][k], F->mod);
        }
        _nmod_poly_dif_inplace_radix2_rec(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec(phi, len/2, order-1, F);
    }

}

void _nmod_poly_dif_inplace_radix2_rec_v2(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        mp_limb_t p7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else
    {
        // here order >= 3, len >= 8, try some unrolling
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/2);
        const mp_ptr ww = F->tab_w[order-2];
        for (ulong k = 0; k < len/2; k+=4)
        {
            const mp_limb_t u0 = p0[k+0];
            const mp_limb_t u1 = p0[k+1];
            const mp_limb_t u2 = p0[k+2];
            const mp_limb_t u3 = p0[k+3];
            const mp_limb_t v0 = p1[k+0];
            const mp_limb_t v1 = p1[k+1];
            const mp_limb_t v2 = p1[k+2];
            const mp_limb_t v3 = p1[k+3];
            p0[k+0] = nmod_add(u0, v0, F->mod);
            p0[k+1] = nmod_add(u1, v1, F->mod);
            p0[k+2] = nmod_add(u2, v2, F->mod);
            p0[k+3] = nmod_add(u3, v3, F->mod);
            NMOD_MUL_PRENORM(p1[k+0], nmod_sub(u0, v0, F->mod) << F->mod.norm, ww[k+0], F->mod);
            NMOD_MUL_PRENORM(p1[k+1], nmod_sub(u1, v1, F->mod) << F->mod.norm, ww[k+1], F->mod);
            NMOD_MUL_PRENORM(p1[k+2], nmod_sub(u2, v2, F->mod) << F->mod.norm, ww[k+2], F->mod);
            NMOD_MUL_PRENORM(p1[k+3], nmod_sub(u3, v3, F->mod) << F->mod.norm, ww[k+3], F->mod);
        }
        _nmod_poly_dif_inplace_radix2_rec_v2(p0, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_v2(p1, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_rec_v3(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        mp_limb_t p7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else
    {
        // here order >= 3, len >= 8, try some unrolling
        // for multiplication:
        mp_limb_t q0xx, q1xx, rxx, p_hixx, p_loxx;
        const mp_limb_t nxx = F->mod.n << F->mod.norm;
        // second half of array
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/2);
        const mp_ptr ww = F->tab_w[order-2];
        // temporary space
        for (ulong k = 0; k < len/2; k+=4)
        {
            const mp_limb_t u0 = p0[k+0];
            const mp_limb_t u1 = p0[k+1];
            const mp_limb_t u2 = p0[k+2];
            const mp_limb_t u3 = p0[k+3];
            const mp_limb_t v0 = p1[k+0];
            const mp_limb_t v1 = p1[k+1];
            const mp_limb_t v2 = p1[k+2];
            const mp_limb_t v3 = p1[k+3];
            p0[k+0] = nmod_add(u0, v0, F->mod);
            p0[k+1] = nmod_add(u1, v1, F->mod);
            p0[k+2] = nmod_add(u2, v2, F->mod);
            p0[k+3] = nmod_add(u3, v3, F->mod);

            umul_ppmm(p_hixx, p_loxx, nmod_sub(u0, v0, F->mod), ww[k+0] << F->mod.norm);
            umul_ppmm(q1xx, q0xx, F->mod.ninv, p_hixx);
            add_ssaaaa(q1xx, q0xx, q1xx, q0xx, p_hixx, p_loxx);
            rxx = (p_loxx - (q1xx + 1) * nxx);
            if (rxx > q0xx)
                rxx += nxx;
            p1[k+0] = (rxx < nxx ? rxx : rxx - nxx) >> F->mod.norm;

            umul_ppmm(p_hixx, p_loxx, nmod_sub(u1, v1, F->mod), ww[k+1] << F->mod.norm);
            umul_ppmm(q1xx, q0xx, F->mod.ninv, p_hixx);
            add_ssaaaa(q1xx, q0xx, q1xx, q0xx, p_hixx, p_loxx);
            rxx = (p_loxx - (q1xx + 1) * nxx);
            if (rxx > q0xx)
                rxx += nxx;
            p1[k+1] = (rxx < nxx ? rxx : rxx - nxx) >> F->mod.norm;

            umul_ppmm(p_hixx, p_loxx, nmod_sub(u2, v2, F->mod), ww[k+2] << F->mod.norm);
            umul_ppmm(q1xx, q0xx, F->mod.ninv, p_hixx);
            add_ssaaaa(q1xx, q0xx, q1xx, q0xx, p_hixx, p_loxx);
            rxx = (p_loxx - (q1xx + 1) * nxx);
            if (rxx > q0xx)
                rxx += nxx;
            p1[k+2] = (rxx < nxx ? rxx : rxx - nxx) >> F->mod.norm;

            umul_ppmm(p_hixx, p_loxx, nmod_sub(u3, v3, F->mod), ww[k+3] << F->mod.norm);
            umul_ppmm(q1xx, q0xx, F->mod.ninv, p_hixx);
            add_ssaaaa(q1xx, q0xx, q1xx, q0xx, p_hixx, p_loxx);
            rxx = (p_loxx - (q1xx + 1) * nxx);
            if (rxx > q0xx)
                rxx += nxx;
            p1[k+3] = (rxx < nxx ? rxx : rxx - nxx) >> F->mod.norm;
        }
        _nmod_poly_dif_inplace_radix2_rec_v3(p0, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_v3(p1, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_rec_v4(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        const mp_limb_t p7 = n_mulmod_shoup(F->tab_w[0][1], nmod_sub(p1, p3, F->mod), F->tab_w_pre[0][1], F->mod.n);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else
    {
        // here order >= 3, len >= 8, try some unrolling
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/2);
        const mp_ptr ww = F->tab_w[order-2];
        const mp_ptr wwpre = F->tab_w_pre[order-2];
        for (ulong k = 0; k < len/2; k+=4)
        {
            const mp_limb_t u0 = p0[k+0];
            const mp_limb_t u1 = p0[k+1];
            const mp_limb_t u2 = p0[k+2];
            const mp_limb_t u3 = p0[k+3];
            const mp_limb_t v0 = p1[k+0];
            const mp_limb_t v1 = p1[k+1];
            const mp_limb_t v2 = p1[k+2];
            const mp_limb_t v3 = p1[k+3];
            p0[k+0] = nmod_add(u0, v0, F->mod);
            p0[k+1] = nmod_add(u1, v1, F->mod);
            p0[k+2] = nmod_add(u2, v2, F->mod);
            p0[k+3] = nmod_add(u3, v3, F->mod);
            p1[k+0] = n_mulmod_shoup(ww[k+0], nmod_sub(u0, v0, F->mod), wwpre[k+0], F->mod.n);
            p1[k+1] = n_mulmod_shoup(ww[k+1], nmod_sub(u1, v1, F->mod), wwpre[k+1], F->mod.n);
            p1[k+2] = n_mulmod_shoup(ww[k+2], nmod_sub(u2, v2, F->mod), wwpre[k+2], F->mod.n);
            p1[k+3] = n_mulmod_shoup(ww[k+3], nmod_sub(u3, v3, F->mod), wwpre[k+3], F->mod.n);
        }
        _nmod_poly_dif_inplace_radix2_rec_v4(p0, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_v4(p1, len/2, order-1, F);
    }

}

void _nmod_poly_dif_inplace_radix2_rec_bench(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = p[0] + p[1];
        p[1] = tmp - p[1];
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = p0 + p2;
        const mp_limb_t p5 = p0 - p2;
        const mp_limb_t p6 = p1 + p3;
        const mp_limb_t p7 = F->tab_w[0][1] * (p1 - p3);
        p[0] = p4 + p6;
        p[1] = p4 - p6;
        p[2] = p5 + p7;
        p[3] = p5 - p7;
    }
    else
    {
        const mp_ptr ww = F->tab_w[order-2];
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/2);
        // here order >= 3, len >= 8
        for (ulong k = 0; k < len/2; k+=4)
        {
            const mp_limb_t u0 = p0[k+0];
            const mp_limb_t u1 = p0[k+1];
            const mp_limb_t u2 = p0[k+2];
            const mp_limb_t u3 = p0[k+3];
            const mp_limb_t v0 = p1[k+0];
            const mp_limb_t v1 = p1[k+1];
            const mp_limb_t v2 = p1[k+2];
            const mp_limb_t v3 = p1[k+3];
            p[k+0] = u0 + v0;
            p[k+1] = u1 + v1;
            p[k+2] = u2 + v2;
            p[k+3] = u3 + v3;
            p1[k+0] = ww[k] * (u0 - v0);
            p1[k+1] = ww[k+1] * (u1 - v1);
            p1[k+2] = ww[k+2] * (u2 - v2);
            p1[k+3] = ww[k+3] * (u3 - v3);
        }
        _nmod_poly_dif_inplace_radix2_rec_bench(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_bench(p+len/2, len/2, order-1, F);
    }
}


/***************************
*  DIF radix 2 iterative  *
***************************/

// iterative version, using NMOD_MUL_PRENORM
void _nmod_poly_dif_inplace_radix2_iter(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
    {
        const mp_ptr ww = F->tab_w[ell-2];
        for (ulong k = 0; k < len; k+=llen)
        {
            const mp_ptr p0 = p+k;
            const mp_ptr p1 = p+(k+llen/2);
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                const mp_limb_t u0 = p0[kk+0];
                const mp_limb_t u1 = p0[kk+1];
                const mp_limb_t u2 = p0[kk+2];
                const mp_limb_t u3 = p0[kk+3];
                const mp_limb_t v0 = p1[kk+0];
                const mp_limb_t v1 = p1[kk+1];
                const mp_limb_t v2 = p1[kk+2];
                const mp_limb_t v3 = p1[kk+3];
                p0[kk+0] = nmod_add(u0, v0, F->mod);
                p0[kk+1] = nmod_add(u1, v1, F->mod);
                p0[kk+2] = nmod_add(u2, v2, F->mod);
                p0[kk+3] = nmod_add(u3, v3, F->mod);
                NMOD_MUL_PRENORM(p1[kk+0], nmod_sub(u0, v0, F->mod) << F->mod.norm, ww[kk+0], F->mod);
                NMOD_MUL_PRENORM(p1[kk+1], nmod_sub(u1, v1, F->mod) << F->mod.norm, ww[kk+1], F->mod);
                NMOD_MUL_PRENORM(p1[kk+2], nmod_sub(u2, v2, F->mod) << F->mod.norm, ww[kk+2], F->mod);
                NMOD_MUL_PRENORM(p1[kk+3], nmod_sub(u3, v3, F->mod) << F->mod.norm, ww[kk+3], F->mod);
            }
        }
    }
    // perform last two FFT layers
    for (ulong k = 0; k < len; k+=4)
    {
        const mp_limb_t p0 = p[k+0];
        const mp_limb_t p1 = p[k+1];
        const mp_limb_t p2 = p[k+2];
        const mp_limb_t p3 = p[k+3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        const mp_limb_t tmp7 = nmod_sub(p1, p3, F->mod);
        mp_limb_t p7;
        NMOD_MUL_PRENORM(p7, tmp7 << F->mod.norm, F->tab_w[0][1], F->mod);
        p[k+0] = nmod_add(p4, p6, F->mod);
        p[k+1] = nmod_sub(p4, p6, F->mod);
        p[k+2] = nmod_add(p5, p7, F->mod);
        p[k+3] = nmod_sub(p5, p7, F->mod);
    }
    //for (ulong k = 0; k < len; k+=2)
    //{
    //    const mp_limb_t tmp = p[k+0];
    //    p[k+0] = nmod_add(p[k+0], p[k+1], F->mod);
    //    p[k+1] = nmod_sub(tmp, p[k+1], F->mod);
    //}
}

// iterative version, using n_mulmod_shoup
void _nmod_poly_dif_inplace_radix2_iter_v2(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
    {
        const mp_ptr ww = F->tab_w[ell-2];
        const mp_ptr wwpre = F->tab_w_pre[ell-2];
        for (ulong k = 0; k < len; k+=llen)
        {
            const mp_ptr p0 = p+k;
            const mp_ptr p1 = p+(k+llen/2);
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                const mp_limb_t u0 = p0[kk+0];
                const mp_limb_t u1 = p0[kk+1];
                const mp_limb_t u2 = p0[kk+2];
                const mp_limb_t u3 = p0[kk+3];
                const mp_limb_t v0 = p1[kk+0];
                const mp_limb_t v1 = p1[kk+1];
                const mp_limb_t v2 = p1[kk+2];
                const mp_limb_t v3 = p1[kk+3];
                p0[kk+0] = nmod_add(u0, v0, F->mod);
                p0[kk+1] = nmod_add(u1, v1, F->mod);
                p0[kk+2] = nmod_add(u2, v2, F->mod);
                p0[kk+3] = nmod_add(u3, v3, F->mod);
                p1[kk+0] = n_mulmod_shoup(ww[kk+0], nmod_sub(u0, v0, F->mod), wwpre[kk+0], F->mod.n);
                p1[kk+1] = n_mulmod_shoup(ww[kk+1], nmod_sub(u1, v1, F->mod), wwpre[kk+1], F->mod.n);
                p1[kk+2] = n_mulmod_shoup(ww[kk+2], nmod_sub(u2, v2, F->mod), wwpre[kk+2], F->mod.n);
                p1[kk+3] = n_mulmod_shoup(ww[kk+3], nmod_sub(u3, v3, F->mod), wwpre[kk+3], F->mod.n);
            }
        }
    }
    // perform last two FFT layers
    for (ulong k = 0; k < len; k+=4)
    {
        const mp_limb_t p0 = p[k+0];
        const mp_limb_t p1 = p[k+1];
        const mp_limb_t p2 = p[k+2];
        const mp_limb_t p3 = p[k+3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        const mp_limb_t tmp7 = nmod_sub(p1, p3, F->mod);
        const mp_limb_t p7 = n_mulmod_shoup(F->tab_w[0][1], tmp7, F->tab_w_pre[0][1], F->mod.n);
        p[k+0] = nmod_add(p4, p6, F->mod);
        p[k+1] = nmod_sub(p4, p6, F->mod);
        p[k+2] = nmod_add(p5, p7, F->mod);
        p[k+3] = nmod_sub(p5, p7, F->mod);
    }
    //for (ulong k = 0; k < len; k+=2)
    //{
    //    const mp_limb_t tmp = p[k+0];
    //    p[k+0] = nmod_add(p[k+0], p[k+1], F->mod);
    //    p[k+1] = nmod_sub(tmp, p[k+1], F->mod);
    //}
}

void _nmod_poly_dif_inplace_radix2_iter_bench(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
    {
        const mp_ptr ww = F->tab_w[ell-2];
        for (ulong k = 0; k < len; k+=llen)
        {
            const mp_ptr p0 = p+k;
            const mp_ptr p1 = p+(k+llen/2);
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                const mp_limb_t u0 = p0[kk+0];
                const mp_limb_t u1 = p0[kk+1];
                const mp_limb_t u2 = p0[kk+2];
                const mp_limb_t u3 = p0[kk+3];
                const mp_limb_t v0 = p1[kk+0];
                const mp_limb_t v1 = p1[kk+1];
                const mp_limb_t v2 = p1[kk+2];
                const mp_limb_t v3 = p1[kk+3];
                p0[kk+0] = u0 + v0;
                p0[kk+1] = u1 + v1;
                p0[kk+2] = u2 + v2;
                p0[kk+3] = u3 + v3;
                p1[kk+0] = ww[kk] * (u0 - v0);
                p1[kk+1] = ww[kk+1] * (u1 - v1);
                p1[kk+2] = ww[kk+2] * (u2 - v2);
                p1[kk+3] = ww[kk+3] * (u3 - v3);
            }
        }
    }
    // perform last two FFT layers
    for (ulong k = 0; k < len; k+=4)
    {
        const mp_limb_t p0 = p[k+0];
        const mp_limb_t p1 = p[k+1];
        const mp_limb_t p2 = p[k+2];
        const mp_limb_t p3 = p[k+3];
        const mp_limb_t p4 = p0 + p2;
        const mp_limb_t p5 = p0 - p2;
        const mp_limb_t p6 = p1 + p3;
        const mp_limb_t p7 = F->tab_w[0][1] * (p1 - p3);
        p[k+0] = p4 + p6;
        p[k+1] = p4 - p6;
        p[k+2] = p5 + p7;
        p[k+3] = p5 - p7;
    }
    //for (ulong k = 0; k < len; k+=2)
    //{
    //    const mp_limb_t tmp = p[k+0];
    //    p[k+0] = p[k+0] + p[k+1];
    //    p[k+1] = tmp - p[k+1];
    //}
}

/***************************
*  DIF radix 4 recursive  *
***************************/

// radix 4, rec, using prenorm
void _nmod_poly_dif_inplace_radix4_rec(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        mp_limb_t p7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else
    {
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/4);
        const mp_ptr p2 = p1+(len/4);
        const mp_ptr p3 = p2+(len/4);
        for (ulong k = 0; k < len/4; k++)
        {
            // for I == w**(2**(order-2)) square root of -1,
            // splitting p in 4 blocks of len/4,
            //                    [1  1  1  1]
            // p <- [p0 p1 p2 p3] [1 -1  I -I]
            //                    [1  1 -1 -1]
            //                    [1 -1 -I  I]
            // p0 = (p0 + p2) + (p1 + p3);
            // p1 = (p0 + p2) - (p1 + p3), multiplied by powers of w**2;
            // p2 = (p0 - p2) + I*(p1 - p3);
            // p3 = (p0 - p2) - I*(p1 - p3), multiplied by powers of w;
            // TODO comment to explain below, mod x^{len/2 - 1}, mod x^{len/2 + 1}, mod x^{len/2 - I}, mod x^{len/2 + I}, 
            const mp_limb_t u0 = p0[k];
            const mp_limb_t u1 = p1[k];
            const mp_limb_t u2 = p2[k];
            const mp_limb_t u3 = p3[k];
            const mp_limb_t u4 = nmod_add(u0, u2, F->mod);
            const mp_limb_t u5 = nmod_sub(u0, u2, F->mod);
            const mp_limb_t u6 = nmod_add(u1, u3, F->mod);
            mp_limb_t u7;
            NMOD_MUL_PRENORM(u7, nmod_sub(u1, u3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
            p0[k] = nmod_add(u4, u6, F->mod);
            NMOD_MUL_PRENORM(p1[k], nmod_sub(u4, u6, F->mod) << F->mod.norm, F->tab_w[order-2][2*k], F->mod);
            NMOD_MUL_PRENORM(p2[k], nmod_add(u5, u7, F->mod) << F->mod.norm, F->tab_w[order-2][k], F->mod);
            if (3*k < len/2)
                NMOD_MUL_PRENORM(p3[k], nmod_sub(u5, u7, F->mod) << F->mod.norm, F->tab_w[order-2][3*k], F->mod);
            else
                NMOD_MUL_PRENORM(p3[k], nmod_sub(u5, u7, F->mod) << F->mod.norm, F->mod.n - F->tab_w[order-2][3*k-len/2], F->mod);
        }
        _nmod_poly_dif_inplace_radix4_rec(p0, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec(p1, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec(p2, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec(p3, len/4, order-2, F);
    }
}

// radix 4, rec, using prenorm, unrolled
void _nmod_poly_dif_inplace_radix4_rec_v2(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        mp_limb_t tmp = p[0];
        p[0] = nmod_add(p[0], p[1], F->mod);
        p[1] = nmod_sub(tmp, p[1], F->mod);
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
        const mp_limb_t p0 = p[0];
        const mp_limb_t p1 = p[1];
        const mp_limb_t p2 = p[2];
        const mp_limb_t p3 = p[3];
        const mp_limb_t p4 = nmod_add(p0, p2, F->mod);
        const mp_limb_t p5 = nmod_sub(p0, p2, F->mod);
        const mp_limb_t p6 = nmod_add(p1, p3, F->mod);
        mp_limb_t p7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else
    {
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/4);
        const mp_ptr p2 = p1+(len/4);
        const mp_ptr p3 = p2+(len/4);
        for (ulong k = 0; k < len/4; k++)
        {
            // for I == w**(2**(order-2)) square root of -1,
            // splitting p in 4 blocks of len/4,
            //                    [1  1  1  1]
            // p <- [p0 p1 p2 p3] [1 -1  I -I]
            //                    [1  1 -1 -1]
            //                    [1 -1 -I  I]
            // p0 = (p0 + p2) + (p1 + p3);
            // p1 = (p0 + p2) - (p1 + p3), multiplied by powers of w**2;
            // p2 = (p0 - p2) + I*(p1 - p3);
            // p3 = (p0 - p2) - I*(p1 - p3), multiplied by powers of w;
            // TODO comment to explain below, mod x^{len/2 - 1}, mod x^{len/2 + 1}, mod x^{len/2 - I}, mod x^{len/2 + I}, 
            const mp_limb_t u0 = p0[k];
            const mp_limb_t u1 = p1[k];
            const mp_limb_t u2 = p2[k];
            const mp_limb_t u3 = p3[k];
            const mp_limb_t u4 = nmod_add(u0, u2, F->mod);
            const mp_limb_t u5 = nmod_sub(u0, u2, F->mod);
            const mp_limb_t u6 = nmod_add(u1, u3, F->mod);
            mp_limb_t u7;
            NMOD_MUL_PRENORM(u7, nmod_sub(u1, u3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
            p0[k] = nmod_add(u4, u6, F->mod);
            NMOD_MUL_PRENORM(p1[k], nmod_sub(u4, u6, F->mod) << F->mod.norm, F->tab_w[order-2][2*k], F->mod);
            NMOD_MUL_PRENORM(p2[k], nmod_add(u5, u7, F->mod) << F->mod.norm, F->tab_w[order-2][k], F->mod);
            if (3*k < len/2)
                NMOD_MUL_PRENORM(p3[k], nmod_sub(u5, u7, F->mod) << F->mod.norm, F->tab_w[order-2][3*k], F->mod);
            else
                NMOD_MUL_PRENORM(p3[k], nmod_sub(u5, u7, F->mod) << F->mod.norm, F->mod.n - F->tab_w[order-2][3*k-len/2], F->mod);
        }
        _nmod_poly_dif_inplace_radix4_rec(p0, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec(p1, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec(p2, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec(p3, len/4, order-2, F);
    }
}

/***************************
*  DIF radix 4 iterative  *
***************************/



/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
