#include <flint/nmod.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_fft.h"

/*----------------*/
/* 2-point DFT    */
/*----------------*/
/** 2-points Discrete Fourier Transform:
 * returns (p(1), p(-1)) for p(x) = a+b*x
 *                    [1   1]
 * [a  b]  <-  [a  b] [1  -1]
 * */

// lazy1: input [0..n), output [0..2n)
#define DFT2_LAZY1(a,b,n)           \
    do {                            \
        const mp_limb_t tmp = (a);  \
        (a) = (a) + (b);            \
        (b) = tmp + (n) - (b);      \
    } while(0)

// lazy2: input [0..2n), output [0..4n)
// n2 is 2*n
#define DFT2_LAZY2(a,b,n2)          \
    do {                            \
        const mp_limb_t tmp = (a);  \
        (a) = (a) + (b);            \
        (b) = tmp + (n2) - (b);     \
    } while(0)

// lazy2 red1:
// input [0..2n) x [0..2n)
// output [0..2n) x [0..4n)
// n2 is 2*n
#define DFT2_LAZY2_RED1(a,b,n2)     \
    do {                            \
        const mp_limb_t tmp = (b);  \
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
        const mp_limb_t tmp = (b);  \
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
        const mp_limb_t tmp = (a);  \
        (a) = (a) + (b);            \
        if ((a) >= (n4))            \
            (a) -= (n4);            \
        (b) = tmp + n4 - (b);       \
        if ((b) >= (n4))            \
            (b) -= (n4);            \
    } while(0)



/*----------------*/
/* 4-point DFT    */
/*----------------*/

/** 4-points DFT:
 * (Decimation In Frequency returns bit reversed order evaluations)
 * returns (p(1), p(-1), p(I), p(-I)) where p(x) = a + b*x + c*x**2 + d*x**3
 * and I is typically a square root of -1 (this property is not exploited)
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 * */

// lazy1: input in [0..n) --> output [0..4*n)
#define DFT4_DIF_SHOUP_LAZY1(a,b,c,d,I,Ipre,n,n2)                  \
    do {                                                           \
        const mp_limb_t p0 = (a);                                  \
        const mp_limb_t p1 = (b);                                  \
        const mp_limb_t p2 = (c);                                  \
        const mp_limb_t p3 = (d);                                  \
        const mp_limb_t p4 = p0 + p2;              /* < 2*n */     \
        const mp_limb_t p5 = p0 + (n) - p2;        /* < 2*n */     \
        const mp_limb_t p6 = p1 + p3;              /* < 2*n */     \
        const mp_limb_t p7 =                       /* < 2*n */     \
             n_mulmod_shoup_lazy((I), p1 + (n) - p3, (Ipre), (n)); \
        (a) = p4 + p6;                             /* < 4*n */     \
        (b) = p4 + (n2) - p6;                     /* < 4*n */      \
        (c) = p5 + p7;                             /* < 4*n */     \
        (d) = p5 + (n2) - p7;                     /* < 4*n */      \
    } while(0)

// lazy2 red: input in [0..2*n) --> output [0..4*n)
// n2 is 2*n
#define DFT4_DIF_SHOUP_LAZY2_RED(a,b,c,d,I,Ipre,n,n2)                  \
    do {                                                               \
        const mp_limb_t p0 = (a);                                      \
        const mp_limb_t p1 = (b);                                      \
        const mp_limb_t p2 = (c);                                      \
        const mp_limb_t p3 = (d);                                      \
        mp_limb_t p4 = p0 + p2;                     /* < 4*n */        \
        if (p4 >= (n2))                                                \
            p4 -= (n2);                             /* < 2*n */        \
        mp_limb_t p5 = p0 + (n2) - p2;              /* < 4*n */        \
        if (p5 >= (n2))                                                \
            p5 -= (n2);                             /* < 2*n */        \
        mp_limb_t p6 = p1 + p3;                     /* < 4*n */        \
        if (p6 >= (n2))                                                \
            p6 -= (n2);                             /* < 2*n */        \
        const mp_limb_t p7 =                        /* < 2*n */        \
             n_mulmod_shoup_lazy((I), p1 + (n2) - p3, (Ipre), (n));    \
        (a) = p4 + p6;                              /* < 4*n */        \
        (b) = p4 + (n2) - p6;                       /* < 4*n */        \
        (c) = p5 + p7;                              /* < 4*n */        \
        (d) = p5 + (n2) - p7;                       /* < 4*n */        \
    } while(0)


/*----------------*/
/* 8-point DFT    */
/*----------------*/

// returns a*b % n  in [0..2*n)
FLINT_FORCE_INLINE ulong n_mulmod_shoup_lazy(ulong a, ulong b, ulong apre, ulong n)
{
    ulong p_hi, p_lo;
    umul_ppmm(p_hi, p_lo, apre, b);
    return a * b - p_hi * n;
}

// reduction tree 8-point lazy DFT
// lazy red: input in [0..2*n) --> output in [0..4*n)
FLINT_FORCE_INLINE void dft8_red_lazy(mp_ptr p, nmod_fft_t F)
{
    ulong p_hi, p_lo;
    mp_limb_t u0 = p[0];
    mp_limb_t u1 = p[1];
    mp_limb_t u2 = p[2];
    mp_limb_t u3 = p[3];
    mp_limb_t v0 = p[4];
    mp_limb_t v1 = p[5];
    mp_limb_t v2 = p[6];
    mp_limb_t v3 = p[7];

    // mod x**4 - 1 | x**4 + 1
    mp_limb_t p0 = u0 + v0;  // [0..4n)
    mp_limb_t p1 = u1 + v1;  // [0..4n)
    mp_limb_t p2 = u2 + v2;  // [0..4n)
    mp_limb_t p3 = u3 + v3;  // [0..4n)
    u0 += F->modn2 - v0;  // [0..4n)
    u1 += F->modn2 - v1;  // [0..4n)
    u2 += F->modn2 - v2;  // [0..4n)
    u3 += F->modn2 - v3;  // [0..4n)

    // left, mod x**2 - 1 | x**2 + 1
    v0 = p0 + p2;             // [0..8n)
    v1 = p1 + p3;             // [0..8n)
    v2 = p0 + F->modn4 - p2;  // [0..8n)
    v3 = p1 + F->modn4 - p3;  // [0..8n)

    // left-left, mod x-1 | x+1
    p0 = v0 + v1;               // [0..16n)
    p1 = v0 + 2*F->modn4 - v1;  // [0..16n)
    if (p0 > 2*F->modn4)
        p0 -= 2*F->modn4;
    if (p0 > F->modn4)
        p0 -= F->modn4;        // [0..4n)
    if (p1 > 2*F->modn4)
        p1 -= 2*F->modn4;
    if (p1 > F->modn4)
        p1 -= F->modn4;        // [0..4n)
    p[0] = p0;
    p[1] = p1;

    // left-right, mod x-I | x+I
    umul_ppmm(p_hi, p_lo, F->Ipre, v3);
    v3 = F->I * v3 - p_hi * F->mod.n;  // [0..2n)
    if (v2 >= F->modn4)
        v2 -= F->modn4;
    if (v2 >= F->modn2)
        v2 -= F->modn2;         // [0..2n)
    p[2] = v2 + v3;             // [0..4n)
    p[3] = v2 + F->modn2 - v3;  // [0..4n)

    // right, mod x**2 - I | x**2 + I
    umul_ppmm(p_hi, p_lo, F->Ipre, u2);
    u2 = F->I * u2 - p_hi * F->mod.n;  // [0..2n)
    umul_ppmm(p_hi, p_lo, F->Ipre, u3);
    u3 = F->I * u3 - p_hi * F->mod.n;  // [0..2n)
    if (u0 >= F->modn2)
        u0 -= F->modn2;         // [0..2n)
    if (u1 >= F->modn2)
        u1 -= F->modn2;         // [0..2n)
    v0 = u0 + u2;  // [0..4n)
    v1 = u1 + u3;  // [0..4n)
    v2 = u0 + F->modn2 - u2;  // [0..4n)
    v3 = u1 + F->modn2 - u3;  // [0..4n)

    // right-left, mod x - J | x + J
    umul_ppmm(p_hi, p_lo, F->Jpre, v1);
    v1 = F->J * v1 - p_hi * F->mod.n;  // [0..2n)
    if (v0 >= F->modn2)
        v0 -= F->modn2;         // [0..2n)
    p[4] = v0 + v1;
    p[5] = v0 + F->modn2 - v1;

    // right-right, mod x - I*J | x + I*J
    umul_ppmm(p_hi, p_lo, F->IJpre, v3);
    v3 = F->IJ * v3 - p_hi * F->mod.n;  // [0..2n)
    if (v2 >= F->modn2)
        v2 -= F->modn2;         // [0..2n)
    p[6] = v2 + v3;
    p[7] = v2 + F->modn2 - v3;
}





// input [0..2*n), output [0..4*n)
void _nmod_fft_dif_rec2_lazy(mp_ptr p, ulong len, ulong order, nmod_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_LAZY3_RED(p[0], p[1], F->modn4);
    else if (order == 2)
        DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->Ipre, F->mod.n, F->modn2);
    else if (order == 3)
        dft8_red_lazy(p, F);
    else
    {
        for (ulong k = 0; k < len/2; k+=4)
        {
            // in: p[k], p[len/2+k] in [0..2n)
            // out: p[k], p[len/2+k] in [0..2n)
            ulong p_hi, p_lo, tmp;

            tmp = p[k+0] + F->modn2 - p[len/2+k+0];
            p[k+0] = p[k+0] + p[len/2+k+0];
            if (p[k+0] >= F->modn2)
                p[k+0] -= F->modn2;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][k+0], tmp);
            p[len/2+k+0] = F->tab_w[order-2][k+0] * tmp - p_hi * F->mod.n;

            tmp = p[k+1] + F->modn2 - p[len/2+k+1];
            p[k+1] = p[k+1] + p[len/2+k+1];
            if (p[k+1] >= F->modn2)
                p[k+1] -= F->modn2;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][k+1], tmp);
            p[len/2+k+1] = F->tab_w[order-2][k+1] * tmp - p_hi * F->mod.n;

            tmp = p[k+2] + F->modn2 - p[len/2+k+2];
            p[k+2] = p[k+2] + p[len/2+k+2];
            if (p[k+2] >= F->modn2)
                p[k+2] -= F->modn2;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][k+2], tmp);
            p[len/2+k+2] = F->tab_w[order-2][k+2] * tmp - p_hi * F->mod.n;

            tmp = p[k+3] + F->modn2 - p[len/2+k+3];
            p[k+3] = p[k+3] + p[len/2+k+3];
            if (p[k+3] >= F->modn2)
                p[k+3] -= F->modn2;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][k+3], tmp);
            p[len/2+k+3] = F->tab_w[order-2][k+3] * tmp - p_hi * F->mod.n;
        }
        _nmod_fft_dif_rec2_lazy(p, len/2, order-1, F);
        _nmod_fft_dif_rec2_lazy(p+len/2, len/2, order-1, F);
    }
}

void _nmod_fft_dif_iter2_lazy(mp_ptr p, ulong len, ulong order, nmod_fft_t F)
{
    // perform FFT layers up to order 3
    ulong llen = len;
    for (ulong ell = order; ell > 3; ell--, llen>>=1)
        for (ulong k = 0; k < len; k+=llen)
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                // in: p[k+kk], p[llen/2+k+kk] in [0..2n)
                // out: p[k+kk], p[llen/2+k+kk] in [0..2n)
                ulong p_hi, p_lo, tmp;

                tmp = p[k+kk+0] + F->modn2 - p[llen/2+k+kk+0];
                p[k+kk+0] = p[k+kk+0] + p[llen/2+k+kk+0];
                if (p[k+kk+0] >= F->modn2)
                    p[k+kk+0] -= F->modn2;
                umul_ppmm(p_hi, p_lo, F->tab_w_pre[ell-2][kk+0], tmp);
                p[llen/2+k+kk+0] = F->tab_w[ell-2][kk+0] * tmp - p_hi * F->mod.n;

                tmp = p[k+kk+1] + F->modn2 - p[llen/2+k+kk+1];
                p[k+kk+1] = p[k+kk+1] + p[llen/2+k+kk+1];
                if (p[k+kk+1] >= F->modn2)
                    p[k+kk+1] -= F->modn2;
                umul_ppmm(p_hi, p_lo, F->tab_w_pre[ell-2][kk+1], tmp);
                p[llen/2+k+kk+1] = F->tab_w[ell-2][kk+1] * tmp - p_hi * F->mod.n;

                tmp = p[k+kk+2] + F->modn2 - p[llen/2+k+kk+2];
                p[k+kk+2] = p[k+kk+2] + p[llen/2+k+kk+2];
                if (p[k+kk+2] >= F->modn2)
                    p[k+kk+2] -= F->modn2;
                umul_ppmm(p_hi, p_lo, F->tab_w_pre[ell-2][kk+2], tmp);
                p[llen/2+k+kk+2] = F->tab_w[ell-2][kk+2] * tmp - p_hi * F->mod.n;

                tmp = p[k+kk+3] + F->modn2 - p[llen/2+k+kk+3];
                p[k+kk+3] = p[k+kk+3] + p[llen/2+k+kk+3];
                if (p[k+kk+3] >= F->modn2)
                    p[k+kk+3] -= F->modn2;
                umul_ppmm(p_hi, p_lo, F->tab_w_pre[ell-2][kk+3], tmp);
                p[llen/2+k+kk+3] = F->tab_w[ell-2][kk+3] * tmp - p_hi * F->mod.n;
            }
    // perform last two FFT layers
    for (ulong k = 0; k < len; k+=8)
        dft8_red_lazy(p+k, F);
    //for (ulong k = 0; k < len; k+=4)
    //    DFT4_DIF_SHOUP_LAZY2_RED(p[k+0], p[k+1], p[k+2], p[k+3], F->tab_w[0][1], F->tab_w_pre[0][1], F->mod.n, F->modn2);
}

void _nmod_fft_red_rec2_lazy(mp_ptr p, ulong len, ulong order, ulong node, nmod_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        if (node==0)  // w == 1
            DFT2_LAZY3_RED(p[0], p[1], F->modn4);
        else
        {
            // in [0..2n), out [0..4n)
            const mp_limb_t u = p[0];
            const mp_limb_t v = n_mulmod_shoup_lazy(F->tab_w[1][node], p[1], F->tab_w_pre[1][node], F->mod.n);
            p[0] = u + v;
            p[1] = u + F->modn2 - v;
        }
    }
    else if (order == 2)
    {
        // in [0..2n), out [0..4n)
        if (node==0)  // w == 1
            DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->tab_w[1][1], F->tab_w_pre[1][1], F->mod.n, F->modn2);
        else
        {
            const mp_limb_t u0 = p[0];
            const mp_limb_t u1 = p[1];
            mp_limb_t v0 = p[2];
            mp_limb_t v1 = p[3];

            mp_limb_t w = F->tab_w[1][node];
            mp_limb_t wpre = F->tab_w_pre[1][node];
            v0 = n_mulmod_shoup_lazy(w, v0, wpre, F->mod.n);
            v1 = n_mulmod_shoup_lazy(w, v1, wpre, F->mod.n);
            mp_limb_t p0 = u0 + v0;  // [0..4n)
            if (p0 >= F->modn2)
                p0 -= F->modn2;      // [0..2n)
            mp_limb_t p1 = u1 + v1;  // [0..4n)
            if (p1 >= F->modn2)
                p1 -= F->modn2;       // [0..2n)
            v0 = u0 + F->modn2 - v0;  // [0..4n)
            if (v0 >= F->modn2)
                v0 -= F->modn2;       // [0..2n)
            v1 = u1 + F->modn2 - v1;  // [0..4n)
            if (v1 >= F->modn2)
                v1 -= F->modn2;       // [0..2n)

            w = F->tab_w[1][2*node];
            wpre = F->tab_w_pre[1][2*node];
            p1 = n_mulmod_shoup_lazy(w, p1, wpre, F->mod.n);
            w = F->tab_w[1][2*node+1];
            wpre = F->tab_w_pre[1][2*node+1];
            v1 = n_mulmod_shoup_lazy(w, v1, wpre, F->mod.n);
            p[0] = p0 + p1;             // [0..4n)
            p[1] = p0 + F->modn2 - p1;  // [0..4n)
            p[2] = v0 + v1;             // [0..4n)
            p[3] = v0 + F->modn2 - v1;  // [0..4n)
        }
    }
    else if (order == 3)
    {
        // in [0..2n), out [0..4n)
        if (node==0)  // w == 1
        {
            mp_limb_t u0 = p[0];
            mp_limb_t u1 = p[1];
            mp_limb_t u2 = p[2];
            mp_limb_t u3 = p[3];
            mp_limb_t v0 = p[4];
            mp_limb_t v1 = p[5];
            mp_limb_t v2 = p[6];
            mp_limb_t v3 = p[7];

            // mod x**4 - 1 | x**4 + 1
            mp_limb_t p0 = u0 + v0;  // [0..4n)
            mp_limb_t p1 = u1 + v1;  // [0..4n)
            mp_limb_t p2 = u2 + v2;  // [0..4n)
            mp_limb_t p3 = u3 + v3;  // [0..4n)
            u0 += F->modn2 - v0;  // [0..4n)
            u1 += F->modn2 - v1;  // [0..4n)
            u2 += F->modn2 - v2;  // [0..4n)
            u3 += F->modn2 - v3;  // [0..4n)

            // left, mod x**2 - 1 | x**2 + 1
            v0 = p0 + p2;             // [0..8n)
            v1 = p1 + p3;             // [0..8n)
            v2 = p0 + F->modn4 - p2;  // [0..8n)
            v3 = p1 + F->modn4 - p3;  // [0..8n)

            // left-left, mod x-1 | x+1
            p0 = v0 + v1;               // [0..16n)
            p1 = v0 + 2*F->modn4 - v1;  // [0..16n)
            if (p0 > 2*F->modn4)
                p0 -= 2*F->modn4;
            if (p0 > F->modn4)
                p0 -= F->modn4;        // [0..4n)
            if (p1 > 2*F->modn4)
                p1 -= 2*F->modn4;
            if (p1 > F->modn4)
                p1 -= F->modn4;        // [0..4n)
            p[0] = p0;
            p[1] = p1;

            // left-right, mod x-I | x+I
            v3 = n_mulmod_shoup_lazy(F->tab_w[1][1], v3, F->tab_w_pre[1][1], F->mod.n);
            if (v2 >= F->modn4)
                v2 -= F->modn4;
            if (v2 >= F->modn2)
                v2 -= F->modn2;         // [0..2n)
            p[2] = v2 + v3;             // [0..4n)
            p[3] = v2 + F->modn2 - v3;  // [0..4n)

            // right, mod x**2 - I | x**2 + I
            u2 = n_mulmod_shoup_lazy(F->tab_w[1][1], u2, F->tab_w_pre[1][1], F->mod.n);
            u3 = n_mulmod_shoup_lazy(F->tab_w[1][1], u3, F->tab_w_pre[1][1], F->mod.n);
            if (u0 >= F->modn2)
                u0 -= F->modn2;         // [0..2n)
            if (u1 >= F->modn2)
                u1 -= F->modn2;         // [0..2n)
            v0 = u0 + u2;  // [0..4n)
            v1 = u1 + u3;  // [0..4n)
            v2 = u0 + F->modn2 - u2;  // [0..4n)
            v3 = u1 + F->modn2 - u3;  // [0..4n)

            // right-left, mod x - J | x + J
            v1 = n_mulmod_shoup_lazy(F->tab_w[1][2], v1, F->tab_w_pre[1][2], F->mod.n);
            if (v0 >= F->modn2)
                v0 -= F->modn2;         // [0..2n)
            p[4] = v0 + v1;
            p[5] = v0 + F->modn2 - v1;

            // right-right, mod x - I*J | x + I*J
            v3 = n_mulmod_shoup_lazy(F->tab_w[1][3], v3, F->tab_w_pre[1][3], F->mod.n);
            if (v2 >= F->modn2)
                v2 -= F->modn2;         // [0..2n)
            p[6] = v2 + v3;
            p[7] = v2 + F->modn2 - v3;
        }
        else
        {
            mp_limb_t u0 = p[0];
            mp_limb_t u1 = p[1];
            mp_limb_t u2 = p[2];
            mp_limb_t u3 = p[3];
            mp_limb_t v0 = p[4];
            mp_limb_t v1 = p[5];
            mp_limb_t v2 = p[6];
            mp_limb_t v3 = p[7];

            // mod x**4 - w | x**4 + w
            mp_limb_t w = F->tab_w[1][node];
            mp_limb_t wpre = F->tab_w_pre[1][node];
            v0 = n_mulmod_shoup_lazy(w, v0, wpre, F->mod.n);
            v1 = n_mulmod_shoup_lazy(w, v1, wpre, F->mod.n);
            v2 = n_mulmod_shoup_lazy(w, v2, wpre, F->mod.n);
            v3 = n_mulmod_shoup_lazy(w, v3, wpre, F->mod.n);
            mp_limb_t p0 = u0 + v0;   // [0..4n)
            mp_limb_t p1 = u1 + v1;   // [0..4n)
            mp_limb_t p2 = u2 + v2;   // [0..4n)
            mp_limb_t p3 = u3 + v3;   // [0..4n)
            u0 = u0 + F->modn2 - v0;  // [0..4n)
            u1 = u1 + F->modn2 - v1;  // [0..4n)
            u2 = u2 + F->modn2 - v2;  // [0..4n)
            u3 = u3 + F->modn2 - v3;  // [0..4n)

            // left, mod x**2 - sqrt(w) | x**2 + sqrt(w)
            w = F->tab_w[1][2*node];
            wpre = F->tab_w_pre[1][2*node];
            p2 = n_mulmod_shoup_lazy(w, p2, wpre, F->mod.n);
            p3 = n_mulmod_shoup_lazy(w, p3, wpre, F->mod.n);
            if (p0 >= F->modn2)
                p0 -= F->modn2;  // [0..2n)
            if (p1 >= F->modn2)
                p1 -= F->modn2;  // [0..2n)
            v0 = p0 + p2;             // [0..4n)
            v1 = p1 + p3;             // [0..4n)
            v2 = p0 + F->modn2 - p2;  // [0..4n)
            v3 = p1 + F->modn2 - p3;  // [0..4n)

            // left-left, mod x - fort(w) | x + fort(w)
            v1 = n_mulmod_shoup_lazy(F->tab_w[1][4*node], v1, F->tab_w_pre[1][4*node], F->mod.n);
            if (v0 >= F->modn2)
                v0 -= F->modn2;  // [0..2n)
            p0 = v0 + v1;             // [0..4n)
            p1 = v0 + F->modn2 - v1;  // [0..4n)

            // left-right, mod x - I*fort(w) | x+ I*fort(w)
            v3 = n_mulmod_shoup_lazy(F->tab_w[1][4*node+1], v3, F->tab_w_pre[1][4*node+1], F->mod.n);
            if (v2 >= F->modn2)
                v2 -= F->modn2;  // [0..2n)
            p2 = v2 + v3;              // [0..4n)
            p3 = v2 + F->modn2 - v3;   // [0..4n)

            // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
            w = F->tab_w[1][2*node+1];
            wpre = F->tab_w_pre[1][2*node+1];
            u2 = n_mulmod_shoup_lazy(w, u2, wpre, F->mod.n);
            u3 = n_mulmod_shoup_lazy(w, u3, wpre, F->mod.n);
            if (u0 >= F->modn2)
                u0 -= F->modn2;  // [0..2n)
            if (u1 >= F->modn2)
                u1 -= F->modn2;  // [0..2n)
            v0 = u0 + u2;             // [0..4n)
            v1 = u1 + u3;             // [0..4n)
            v2 = u0 + F->modn2 - u2;  // [0..4n)
            v3 = u1 + F->modn2 - u3;  // [0..4n)

            // right-left, mod x - J*fort(w) | x + J*fort(w)
            v1 = n_mulmod_shoup_lazy(F->tab_w[1][4*node+2], v1, F->tab_w_pre[1][4*node+2], F->mod.n);
            if (v0 >= F->modn2)
                v0 -= F->modn2;  // [0..2n)
            u0 = v0 + v1;             // [0..4n)
            u1 = v0 + F->modn2 - v1;  // [0..4n)

            // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
            v3 = n_mulmod_shoup_lazy(F->tab_w[1][4*node+3], v3, F->tab_w_pre[1][4*node+3], F->mod.n);
            if (v2 >= F->modn2)
                v2 -= F->modn2;  // [0..2n)
            u2 = v2 + v3;             // [0..4n)
            u3 = v2 + F->modn2 - v3;  // [0..4n)

            p[0] = p0;
            p[1] = p1;
            p[2] = p2;
            p[3] = p3;
            p[4] = u0;
            p[5] = u1;
            p[6] = u2;
            p[7] = u3;
        }
    }
    else
    {
        // in: [0..4n), out: [0..4n)
        if (node==0)  // w == 1
            for (ulong k = 0; k < len/2; k++)
                DFT2_LAZY2_RED1(p[k+0], p[len/2+k+0], F->modn2);
        else
        {
            const mp_limb_t w = F->tab_w[1][node];
            const mp_limb_t wpre = F->tab_w_pre[1][node];
            for (ulong k = 0; k < len/2; k+=4)
            {
                ulong p_hi, p_lo, u, v;

                u = p[k+0];
                if (u >= F->modn2)
                    u -= F->modn2;  // [0..2n)
                v = p[len/2+k+0];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod.n;
                p[k+0] = u + v;                   // [0..4n)
                p[len/2+k+0] = u + F->modn2 - v;  // [0..4n)

                u = p[k+1];
                if (u >= F->modn2)
                    u -= F->modn2;  // [0..2n)
                v = p[len/2+k+1];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod.n;
                p[k+1] = u + v;  // [0..4n)
                p[len/2+k+1] = u + F->modn2 - v;  // [0..4n)

                u = p[k+2];
                if (u >= F->modn2)
                    u -= F->modn2;  // [0..2n)
                v = p[len/2+k+2];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod.n;
                p[k+2] = u + v;  // [0..4n)
                p[len/2+k+2] = u + F->modn2 - v;  // [0..4n)

                u = p[k+3];
                if (u >= F->modn2)
                    u -= F->modn2;  // [0..2n)
                v = p[len/2+k+3];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod.n;
                p[k+3] = u + v;  // [0..4n)
                p[len/2+k+3] = u + F->modn2 - v;  // [0..4n)
            }
        }
        _nmod_fft_red_rec2_lazy(p, len/2, order-1, 2*node, F);
        _nmod_fft_red_rec2_lazy(p+len/2, len/2, order-1, 2*node+1, F);
    }
}

// in [0..2n) out [0..4n)
void _nmod_poly_dif_rec4_lazy(mp_ptr p, ulong len, ulong order, nmod_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        // in [0..2n) out [0..4n)
        DFT2_LAZY3_RED(p[0], p[1], F->modn2);
    else if (order == 2)
        // in [0..2n) out [0..4n)
        DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->tab_w[0][1], F->tab_w_pre[0][1], F->mod.n, F->modn2);
    else if (order == 3)
    {
        dft8_red_lazy(p, F);
        //dft8_dif_lazy(p, F);
    }
    else
    {
        // in [0..2n) out [0..2n)
        const mp_ptr p0 = p;
        const mp_ptr p1 = p+(len/4);
        const mp_ptr p2 = p1+(len/4);
        const mp_ptr p3 = p2+(len/4);
        for (ulong k = 0; k < len/4; k++)
        {
            const mp_limb_t u0 = p0[k];
            const mp_limb_t u1 = p1[k];
            const mp_limb_t u2 = p2[k];
            const mp_limb_t u3 = p3[k];
            mp_limb_t u4 = u0 + u2;  // [0..4n)
            mp_limb_t u5 = u0 + F->modn2 - u2;  // [0..4n)
            mp_limb_t u6 = u1 + u3;  // [0..4n)
            mp_limb_t u7 = u1 + F->modn2 - u3;

            ulong p_hi, p_lo;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[0][1], u7);
            u7 = F->tab_w[0][1] * u7 - p_hi * F->mod.n;  // [0..2n)

            p_lo = u4 + u6;  // [0..8n)
            if (p_lo >= F->modn4)
                p_lo -= F->modn4;
            if (p_lo >= F->modn2)
                p_lo -= F->modn2;
            p0[k] = p_lo;  // [0..2n)

            u4 = u4 + F->modn4 - u6;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][2*k], u4);
            p1[k] = F->tab_w[order-2][2*k] * u4 - p_hi * F->mod.n;  // [0..2n)

            u6 = u5 + u7;
            umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][k], u6);
            p2[k] = F->tab_w[order-2][k] * u6 - p_hi * F->mod.n;  // [0..2n)

            if (3*k < len/2)
            {
                u5 = u5 + F->modn2 - u7;
                umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][3*k], u5);
                p3[k] = F->tab_w[order-2][3*k] * u5 - p_hi * F->mod.n;  // [0..2n)
            }
            else
            {
                u5 = u7 + F->modn4 - u5;
                umul_ppmm(p_hi, p_lo, F->tab_w_pre[order-2][3*k-len/2], u5);
                p3[k] = F->tab_w[order-2][3*k-len/2] * u5 - p_hi * F->mod.n;  // [0..2n)
            }
        }
        _nmod_poly_dif_rec4_lazy(p0, len/4, order-2, F);
        _nmod_poly_dif_rec4_lazy(p1, len/4, order-2, F);
        _nmod_poly_dif_rec4_lazy(p2, len/4, order-2, F);
        _nmod_poly_dif_rec4_lazy(p3, len/4, order-2, F);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
