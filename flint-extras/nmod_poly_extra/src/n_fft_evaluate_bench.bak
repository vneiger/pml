#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include <flint/n_fft.h>

/*************************************************************************
*  the below functions are here for bench: they don't do any modular
*  reduction so somehow estimate what we lose by doing them
*
*  Very incomplete:
*  - many variants are missing (reduction tree, iterative, radix 4,...)
*  - would require more base cases up to length ~ 32
*************************************************************************/


/*------------------------------------------------------------*/
/* 2-point butterfly                                          */
/*------------------------------------------------------------*/

/** 2-points Discrete Fourier Transform:
 * returns (p(1), p(-1)) for p(x) = a+b*x
 *                    [1   1]
 * [a  b]  <-  [a  b] [1  -1]
 * */
#define DFT2_NMOD(a,b,mod)                \
    do {                                  \
        const ulong tmp = (a);        \
        (a) = nmod_add((a), (b), (mod));  \
        (b) = nmod_sub(tmp, (b), (mod));  \
    } while(0)

// lazy1: input [0..n), output [0..2n)
#define DFT2_LAZY1(a,b,n)           \
    do {                            \
        const ulong tmp = (a);  \
        (a) = (a) + (b);            \
        (b) = tmp + (n) - (b);      \
    } while(0)

// lazy2: input [0..2n), output [0..4n)
// n2 is 2*n
#define DFT2_LAZY2(a,b,n2)          \
    do {                            \
        const ulong tmp = (a);  \
        (a) = (a) + (b);            \
        (b) = tmp + (n2) - (b);     \
    } while(0)


// lazy2 red1:
// input [0..2n) x [0..2n)
// output [0..2n) x [0..4n)
// n2 is 2*n
#define DFT2_LAZY2_RED1(a,b,n2)     \
    do {                            \
        const ulong tmp = (b);  \
        (b) = (a) + (n2) - tmp;     \
        (a) = (a) + tmp;            \
        if ((a) >= (n2))            \
            (a) -= (n2);            \
    } while(0)

// lazy2 red:
// input [0..2n), output [0..2n)
// n2 is 2*n
#define DFT2_LAZY2_RED(a,b,n2)      \
    do {                            \
        const ulong tmp = (b);  \
        (b) = (a) + (n2) - tmp;     \
        if ((b) >= (n2))            \
            (b) -= (n2);            \
        (a) += tmp;                 \
        if ((a) >= (n2))            \
            (a) -= (n2);            \
    } while(0)


// lazy3 red:
// input [0..4n), output [0..4n)
// n4 is 4*n
#define DFT2_LAZY3_RED(a,b,n4)      \
    do {                            \
        const ulong tmp = (a);  \
        (a) = (a) + (b);            \
        if ((a) >= (n4))            \
            (a) -= (n4);            \
        (b) = tmp + n4 - (b);       \
        if ((b) >= (n4))            \
            (b) -= (n4);            \
    } while(0)

/*------------------------------------------------------------*/
/* 4-point DFT (DIF style)                                    */
/*------------------------------------------------------------*/

/*****************************************
*  TEMPORARIES: BENCH WITHOUT REDUCTION  *
*****************************************/

/** 4-points DFT (DIF style):
 * (Decimation In Frequency returns bit reversed order evaluations)
 * returns (p(1), p(-1), p(I), p(-I)) where p(x) = a + b*x + c*x**2 + d*x**3
 * and I is typically a square root of -1 (this property is not exploited)
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 * */

// 2-points DFT, bench version without mod:
#define DFT2_BENCH(a,b)         \
    do {                        \
        ulong tmp = (a);    \
        (a) = (a) + (b);        \
        (b) = tmp - (b);        \
    } while(0)

// 4-points DFT (DIF style):
#define DFT4_DIF_BENCH(a,b,c,d,I)            \
    do {                                      \
        const ulong p0 = (a);             \
        const ulong p1 = (b);             \
        const ulong p2 = (c);             \
        const ulong p3 = (d);             \
        const ulong p4 = p0 + p2;         \
        const ulong p5 = p0 - p2;         \
        const ulong p6 = p1 + p3;         \
        const ulong p7 = (I) * (p1 - p3); \
        (a) = p4 + p6;                        \
        (b) = p4 - p6;                        \
        (c) = p5 + p7;                        \
        (d) = p5 - p7;                        \
    } while(0)

/***************************
*  DIF radix 2 recursive  *
***************************/

void _nmod_poly_dif_inplace_radix2_rec_bench(nn_ptr p, ulong len, ulong order, n_fft_ctx_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_BENCH(p[0], p[1]);
    else if (order == 2)
        DFT4_DIF_BENCH(p[0], p[1], p[2], p[3], F->tab_w[2]);
    else
    {
        const nn_ptr ww = F->tab_w + 2*(order-2); // !!!order-2 not checked!!!
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+(len/2);
        // here order >= 3, len >= 8
        for (ulong k = 0; k < len/2; k+=4)
        {
            const ulong u0 = p0[k+0];
            const ulong u1 = p0[k+1];
            const ulong u2 = p0[k+2];
            const ulong u3 = p0[k+3];
            const ulong v0 = p1[k+0];
            const ulong v1 = p1[k+1];
            const ulong v2 = p1[k+2];
            const ulong v3 = p1[k+3];
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

// BENCH with Shoup's stuff but without "if corrections"
//void _nmod_poly_dif_inplace_radix2_rec_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
//{
//    // order == 0: nothing to do
//    if (order == 1)
//        DFT2_BENCH(p[0], p[1]);
//    else if (order == 2)
//        DFT4_DIF_BENCH(p[0], p[1], p[2], p[3], F->tab_w[0][1]);
//    else
//    {
//        const nn_ptr ww = F->tab_w[order-2];
//        const nn_ptr p0 = p;
//        const nn_ptr p1 = p+(len/2);
//        // here order >= 3, len >= 8
//        for (ulong k = 0; k < len/2; k+=4)
//        {
//            const ulong u0 = p0[k+0];
//            const ulong u1 = p0[k+1];
//            const ulong u2 = p0[k+2];
//            const ulong u3 = p0[k+3];
//            const ulong v0 = p1[k+0];
//            const ulong v1 = p1[k+1];
//            const ulong v2 = p1[k+2];
//            const ulong v3 = p1[k+3];
//            p[k+0] = u0 + v0;
//            p[k+1] = u1 + v1;
//            p[k+2] = u2 + v2;
//            p[k+3] = u3 + v3;
//
//            ulong p_hi, p_lo;
//            umul_ppmm(p_hi, p_lo, ww[k], u0-v0);
//            p1[k+0] = ww[k] * p_hi - u0 * F->mod.n;
//
//            umul_ppmm(p_hi, p_lo, ww[k+1], u1-v1);
//            p1[k+1] = ww[k+1] * p_hi - u1 * F->mod.n;
//
//            umul_ppmm(p_hi, p_lo, ww[k+2], u0-v0);
//            p1[k+2] = ww[k+2] * p_hi - u2 * F->mod.n;
//
//            umul_ppmm(p_hi, p_lo, ww[k+3], u0-v0);
//            p1[k+3] = ww[k+3] * p_hi - u3 * F->mod.n;
//        }
//        _nmod_poly_dif_inplace_radix2_rec_bench(p, len/2, order-1, F);
//        _nmod_poly_dif_inplace_radix2_rec_bench(p+len/2, len/2, order-1, F);
//    }
//}

/***************************
*  DIF radix 2 iterative  *
***************************/

void _nmod_poly_dif_inplace_radix2_iter_bench(nn_ptr p, ulong len, ulong order, n_fft_ctx_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
    {
        const nn_ptr ww = F->tab_w + 2*(ell-2);  // !!! ell - 2 not checked !!!
        for (ulong k = 0; k < len; k+=llen)
        {
            const nn_ptr p0 = p+k;
            const nn_ptr p1 = p+(k+llen/2);
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                const ulong u0 = p0[kk+0];
                const ulong u1 = p0[kk+1];
                const ulong u2 = p0[kk+2];
                const ulong u3 = p0[kk+3];
                const ulong v0 = p1[kk+0];
                const ulong v1 = p1[kk+1];
                const ulong v2 = p1[kk+2];
                const ulong v3 = p1[kk+3];
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
        DFT4_DIF_BENCH(p[k+0], p[k+1], p[k+2], p[k+3], F->tab_w[2]);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
