#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <immintrin.h>

#include "nmod_poly_fft.h"

/*--------------------------------*/
/* Modular multiplication helpers */
/*--------------------------------*/

// Shoup's modular multiplication with precomputation, lazy
// let r = (a*b) % n, this returns r or r+n
// result stored in res; a_pr is the precomputation for n, p_hi and p_lo are temporaries
#define N_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n, p_hi, p_lo) \
    do {                                                      \
        umul_ppmm(p_hi, p_lo, (a_pr), (b));                   \
        res = (a) * (b) - p_hi * (n);                         \
    } while(0)

/*----------------*/
/* 2-point DFT    */
/*----------------*/
/** 2-points Discrete Fourier Transform:
 * returns (p(1), p(-1)) for p(x) = a+b*x
 *                    [1   1]
 * [a  b]  <-  [a  b] [1  -1]
 * */

// lazy2 red1:
// input [0..2n) x [0..2n)
// output [0..2n) x [0..4n)
// n2 is 2*n, tmp is a temporary (ulong)
#define DFT2_LAZY2_RED1(a, b, n2, tmp) \
    do {                               \
        tmp = (b);                     \
        (b) = (a) + (n2) - tmp;        \
        (a) = (a) + tmp;               \
        if ((a) >= (n2))               \
            (a) -= (n2);               \
    } while(0)

// lazy2 red1 out:
// input [0..2n) x [0..2n)
// output [0..2n) x [0..4n)
// b unchanged, a-b stored in tmp
// n2 is 2*n
#define DFT2_LAZY2_RED1_OUT(a, b, n2, tmp) \
    do {                                   \
        tmp = (a) + (n2) - (b);            \
        (a) = (a) + (b);                   \
        if ((a) >= (n2))                   \
            (a) -= (n2);                   \
    } while(0)

// Gentleman-Sande butterfly, in-place, lazy
// in:  a,b in [0..2n) x [0..2n)
// out: a,b in [0..2n) x [0..2n)
#define BUTTERFLY_GS_LAZY(a, b, n, n2, w, w_pr, p_hi, p_lo, tmp) \
    do {                                                         \
        DFT2_LAZY2_RED1_OUT(a, b, n2, tmp);                      \
        N_MULMOD_PRECOMP_LAZY(b, w, tmp, w_pr, n, p_hi, p_lo);   \
    } while(0)

// lazy3 red:
// input [0..4n), output [0..4n)
// n4 is 4*n
// currently only used in base cases
#define DFT2_LAZY3_RED(a, b, n4) \
    do {                         \
        const ulong tmp = (a);   \
        (a) = (a) + (b);         \
        if ((a) >= (n4))         \
            (a) -= (n4);         \
        (b) = tmp + n4 - (b);    \
        if ((b) >= (n4))         \
            (b) -= (n4);         \
    } while(0)




// OTHER DFT2 NOT USED CURRENTLY

// lazy1: input [0..n), output [0..2n)
#define DFT2_LAZY1(a,b,n)           \
    do {                            \
        const ulong u = (a);    \
        const ulong v = (b);    \
        (a) = u + v;                \
        (b) = u + (n) - v;          \
    } while(0)

// lazy2: input [0..2n), output [0..4n)
// n2 is 2*n
#define DFT2_LAZY2(a,b,n2)          \
    do {                            \
        const ulong tmp = (a);  \
        (a) = (a) + (b);            \
        (b) = tmp + (n2) - (b);     \
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

// lazy2 red: input in [0..2*n) --> output [0..4*n)
// n2 is 2*n
#define DFT4_DIF_SHOUP_LAZY2_RED(a,b,c,d,I,Ipre,n,n2)                  \
    do {                                                               \
        const ulong p0 = (a);                                          \
        const ulong p1 = (b);                                          \
        const ulong p2 = (c);                                          \
        const ulong p3 = (d);                                          \
        ulong p4 = p0 + p2;                         /* < 4*n */        \
        if (p4 >= (n2))                                                \
            p4 -= (n2);                             /* < 2*n */        \
        ulong p5 = p0 + (n2) - p2;                  /* < 4*n */        \
        if (p5 >= (n2))                                                \
            p5 -= (n2);                             /* < 2*n */        \
        ulong p6 = p1 + p3;                         /* < 4*n */        \
        if (p6 >= (n2))                                                \
            p6 -= (n2);                             /* < 2*n */        \
        ulong p7, p_hi, p_lo;                                          \
        N_MULMOD_PRECOMP_LAZY(p7, (I), p1 + (n2) - p3, (Ipre), (n),    \
                              p_hi, p_lo);                             \
        (a) = p4 + p6;                              /* < 4*n */        \
        (b) = p4 + (n2) - p6;                       /* < 4*n */        \
        (c) = p5 + p7;                              /* < 4*n */        \
        (d) = p5 + (n2) - p7;                       /* < 4*n */        \
    } while(0)

// in [0..2n) out [0..2n)
#define BUTTERFLY4_GS_LAZY(a,b,c,d,                                    \
                           I,I_pr,w1,w1_pr,w2,w2_pr,w3,w3_pr,          \
                           n,n2,n4,p_hi,p_lo)                          \
    do {                                                               \
            const ulong u0 = (a);                                      \
            const ulong u1 = (b);                                      \
            const ulong u2 = (c);                                      \
            const ulong u3 = (d);                                      \
                                                                       \
            ulong u4 = u0 + u2;            /* [0..4n) */               \
            ulong u5 = u0 + n2 - u2;       /* [0..4n) */               \
            ulong u6 = u1 + u3;            /* [0..4n) */               \
            ulong u7 = u1 + n2 - u3;       /* [0..4n) */               \
                                                                       \
            N_MULMOD_PRECOMP_LAZY(u7, I, u7, I_pr, n, p_hi, p_lo);     \
                                                                       \
            p_lo = u4 + u6;                /* [0..8n) */               \
            if (p_lo >= n4)                                            \
                p_lo -= n4;                                            \
            if (p_lo >= n2)                                            \
                p_lo -= n2;                                            \
            (a) = p_lo;                    /* [0..2n) */               \
                                                                       \
            u4 = u4 + n4 - u6;                                         \
            N_MULMOD_PRECOMP_LAZY((b), w2, u4, w2_pr, n, p_hi, p_lo);  \
            u6 = u5 + u7;                                              \
            N_MULMOD_PRECOMP_LAZY((c), w1, u6, w1_pr, n, p_hi, p_lo);  \
            u5 = u5 + n2 - u7;                                         \
            N_MULMOD_PRECOMP_LAZY((d), w3, u5, w3_pr, n, p_hi, p_lo);  \
    } while(0)


// NOT USED
// lazy1: input in [0..n) --> output [0..4*n)
#define DFT4_DIF_SHOUP_LAZY1(a,b,c,d,I,Ipre,n,n2)                  \
    do {                                                           \
        const ulong p0 = (a);                                      \
        const ulong p1 = (b);                                      \
        const ulong p2 = (c);                                      \
        const ulong p3 = (d);                                      \
        const ulong p4 = p0 + p2;                  /* < 2*n */     \
        const ulong p5 = p0 + (n) - p2;            /* < 2*n */     \
        const ulong p6 = p1 + p3;                  /* < 2*n */     \
        const ulong p7 =                           /* < 2*n */     \
             n_mulmod_shoup_lazy((I), p1 + (n) - p3, (Ipre), (n)); \
        (a) = p4 + p6;                             /* < 4*n */     \
        (b) = p4 + (n2) - p6;                      /* < 4*n */     \
        (c) = p5 + p7;                             /* < 4*n */     \
        (d) = p5 + (n2) - p7;                      /* < 4*n */     \
    } while(0)


/*----------------*/
/* 8-point DFT    */
/*----------------*/

// reduction tree 8-point lazy DFT
// lazy red: input in [0..2*n) --> output in [0..4*n)
// FIXME goes through [0..16n)... easily improved to [0..8n)
//       (just reduce v0,v1 to [0..4n) immediately --> check no speed difference)
FLINT_FORCE_INLINE void dft8_red_lazy(nn_ptr p, nmod_fft_ctx_t F)
{
    ulong p_hi, p_lo;
    ulong u0 = p[0];
    ulong u1 = p[1];
    ulong u2 = p[2];
    ulong u3 = p[3];
    ulong v0 = p[4];
    ulong v1 = p[5];
    ulong v2 = p[6];
    ulong v3 = p[7];

    // mod x**4 - 1 | x**4 + 1
    ulong p0 = u0 + v0;  // [0..4n)
    ulong p1 = u1 + v1;  // [0..4n)
    ulong p2 = u2 + v2;  // [0..4n)
    ulong p3 = u3 + v3;  // [0..4n)
    u0 += F->mod2 - v0;  // [0..4n)
    u1 += F->mod2 - v1;  // [0..4n)
    u2 += F->mod2 - v2;  // [0..4n)
    u3 += F->mod2 - v3;  // [0..4n)

    // left, mod x**2 - 1 | x**2 + 1
    v0 = p0 + p2;             // [0..8n)
    v1 = p1 + p3;             // [0..8n)
    v2 = p0 + F->mod4 - p2;  // [0..8n)
    v3 = p1 + F->mod4 - p3;  // [0..8n)

    // left-left, mod x-1 | x+1
    p0 = v0 + v1;               // [0..16n)
    p1 = v0 + 2*F->mod4 - v1;  // [0..16n)
    if (p0 > 2*F->mod4)
        p0 -= 2*F->mod4;
    if (p0 > F->mod4)
        p0 -= F->mod4;        // [0..4n)
    if (p1 > 2*F->mod4)
        p1 -= 2*F->mod4;
    if (p1 > F->mod4)
        p1 -= F->mod4;        // [0..4n)
    p[0] = p0;
    p[1] = p1;

    // left-right, mod x-I | x+I
    N_MULMOD_PRECOMP_LAZY(v3, F->I, v3, F->Ipre, F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;         // [0..2n)
    p[2] = v2 + v3;             // [0..4n)
    p[3] = v2 + F->mod2 - v3;  // [0..4n)

    // right, mod x**2 - I | x**2 + I
    N_MULMOD_PRECOMP_LAZY(u2, F->I, u2, F->Ipre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(u3, F->I, u3, F->Ipre, F->mod, p_hi, p_lo);
    if (u0 >= F->mod2)
        u0 -= F->mod2;         // [0..2n)
    if (u1 >= F->mod2)
        u1 -= F->mod2;         // [0..2n)
    v0 = u0 + u2;  // [0..4n)
    v1 = u1 + u3;  // [0..4n)
    v2 = u0 + F->mod2 - u2;  // [0..4n)
    v3 = u1 + F->mod2 - u3;  // [0..4n)

    // right-left, mod x - J | x + J
    N_MULMOD_PRECOMP_LAZY(v1, F->J, v1, F->Jpre, F->mod, p_hi, p_lo);
    if (v0 >= F->mod2)
        v0 -= F->mod2;         // [0..2n)
    p[4] = v0 + v1;
    p[5] = v0 + F->mod2 - v1;

    // right-right, mod x - I*J | x + I*J
    N_MULMOD_PRECOMP_LAZY(v3, F->IJ, v3, F->IJpre, F->mod, p_hi, p_lo);
    if (v2 >= F->mod2)
        v2 -= F->mod2;         // [0..2n)
    p[6] = v2 + v3;
    p[7] = v2 + F->mod2 - v3;
}


// input [0..2*n), output [0..4*n)
// order >= 3
void _nmod_fft_dif_rec2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // order == 0: nothing to do
    //if (order == 1)
    //    DFT2_LAZY3_RED(p[0], p[1], F->mod4);
    //else if (order == 2)
    //    DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->Ipre, F->mod, F->mod2);
    //else
    if (order == 3)
        dft8_red_lazy(p, F);
    else
    {
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        const nn_ptr ww = F->tab_w[order-4];
        // in: p0, p1 in [0..2n); out: p0, p1 in [0..2n)
        for (ulong k = 0; k < len/2; k+=4)
        {
            ulong p_hi, p_lo, tmp;
            BUTTERFLY_GS_LAZY(p0[k+0], p1[k+0], F->mod, F->mod2, ww[2*k+0], ww[2*k+1], p_hi, p_lo, tmp);
            BUTTERFLY_GS_LAZY(p0[k+1], p1[k+1], F->mod, F->mod2, ww[2*k+2], ww[2*k+3], p_hi, p_lo, tmp);
            BUTTERFLY_GS_LAZY(p0[k+2], p1[k+2], F->mod, F->mod2, ww[2*k+4], ww[2*k+5], p_hi, p_lo, tmp);
            BUTTERFLY_GS_LAZY(p0[k+3], p1[k+3], F->mod, F->mod2, ww[2*k+6], ww[2*k+7], p_hi, p_lo, tmp);
        }
        _nmod_fft_dif_rec2_lazy(p0, len/2, order-1, F);
        _nmod_fft_dif_rec2_lazy(p1, len/2, order-1, F);
    }
}

// input [0..2*n), output [0..4*n)
// order >= 3
void _nmod_fft_dif_rec2_lazy_new(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // order == 0: nothing to do
    //if (order == 1)
    //    DFT2_LAZY3_RED(p[0], p[1], F->mod4);
    //else if (order == 2)
    //    DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->Ipre, F->mod, F->mod2);
    //else
    if (order == 3)
        dft8_red_lazy(p, F);
    else
    {
        for (ulong k = 0; k < len/2; k+=4)
        {
            // in: p[k], p[len/2+k] in [0..2n)
            // out: p[k], p[len/2+k] in [0..2n)
            ulong p_hi, p_lo, tmp;

            tmp = p[k+0] + F->mod2 - p[len/2+k+0];
            p[k+0] = p[k+0] + p[len/2+k+0];
            if (p[k+0] >= F->mod2)
                p[k+0] -= F->mod2;
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][2*k+1], tmp);
            p[len/2+k+0] = F->tab_w[order-4][2*k+0] * tmp - p_hi * F->mod;

            tmp = p[k+1] + F->mod2 - p[len/2+k+1];
            p[k+1] = p[k+1] + p[len/2+k+1];
            if (p[k+1] >= F->mod2)
                p[k+1] -= F->mod2;
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][2*k+3], tmp);
            p[len/2+k+1] = F->tab_w[order-4][2*k+2] * tmp - p_hi * F->mod;

            tmp = p[k+2] + F->mod2 - p[len/2+k+2];
            p[k+2] = p[k+2] + p[len/2+k+2];
            if (p[k+2] >= F->mod2)
                p[k+2] -= F->mod2;
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][2*k+5], tmp);
            p[len/2+k+2] = F->tab_w[order-4][2*k+4] * tmp - p_hi * F->mod;

            tmp = p[k+3] + F->mod2 - p[len/2+k+3];
            p[k+3] = p[k+3] + p[len/2+k+3];
            if (p[k+3] >= F->mod2)
                p[k+3] -= F->mod2;
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][2*k+7], tmp);
            p[len/2+k+3] = F->tab_w[order-4][2*k+6] * tmp - p_hi * F->mod;
        }
        _nmod_fft_dif_rec2_lazy(p, len/2, order-1, F);
        _nmod_fft_dif_rec2_lazy(p+len/2, len/2, order-1, F);
    }
}

// order >= 3
void _nmod_fft_dif_iter2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // perform FFT layers up to order 3
    ulong llen = len;
    for (ulong ell = order; ell > 3; ell--, llen>>=1)
        for (ulong k = 0; k < len; k+=llen)
        {
            const nn_ptr p0 = p+k;
            const nn_ptr p1 = p+llen/2+k;
            const nn_ptr ww = F->tab_w[ell-4];
            ulong p_hi, p_lo, tmp;
            // in: p0, p1 in [0..2n); out: p0, p1 in [0..2n)
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                BUTTERFLY_GS_LAZY(p0[kk+0], p1[kk+0], F->mod, F->mod2, ww[2*kk+0], ww[2*kk+1], p_hi, p_lo, tmp);
                BUTTERFLY_GS_LAZY(p0[kk+1], p1[kk+1], F->mod, F->mod2, ww[2*kk+2], ww[2*kk+3], p_hi, p_lo, tmp);
                BUTTERFLY_GS_LAZY(p0[kk+2], p1[kk+2], F->mod, F->mod2, ww[2*kk+4], ww[2*kk+5], p_hi, p_lo, tmp);
                BUTTERFLY_GS_LAZY(p0[kk+3], p1[kk+3], F->mod, F->mod2, ww[2*kk+6], ww[2*kk+7], p_hi, p_lo, tmp);
            }
        }

    // perform last FFT layers
    for (ulong k = 0; k < len; k+=8)
        dft8_red_lazy(p+k, F);
}

// input [0..2n),  output [0..4n)
// restricting to order >= 3 is a bit faster for smallish orders
void _nmod_fft_red_rec2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // order == 0: nothing to do
    //if (order == 1)
    //    // in [0..4n), out [0..4n)
    //    DFT2_LAZY3_RED(p[0], p[1], F->mod4);
    //else if (order == 2)
    //    // in [0..2n), out [0..4n)
    //    DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->Ipre, F->mod, F->mod2);
    //else
    if (order == 3)
        // in [0..2n), out [0..4n)
        dft8_red_lazy(p, F);
    else
    {
        // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
        // (general accepts [0..4n) as input for order >= 3)
        ulong tmp;
        for (ulong k = 0; k < len/2; k++)
            DFT2_LAZY2_RED1(p[k+0], p[len/2+k+0], F->mod2, tmp);
        _nmod_fft_red_rec2_lazy(p, len/2, order-1, F);
        _nmod_fft_red_rec2_lazy_general(p+len/2, len/2, order-1, 1, F);
    }
}


// if order < 3, in [0..2n) out [0..4n)
// if order >= 3, in [0..4n) out [0..4n)
// restricting to order >= 3 is a bit faster for smallish orders
void _nmod_fft_red_rec2_lazy_general(nn_ptr p, ulong len, ulong order, ulong node, nmod_fft_ctx_t F)
{
    // order == 0: nothing to do
    //if (order == 1)
    //{
    //    // in [0..2n), out [0..4n)
    //    const ulong u = p[0];
    //    ulong v = p[1];
    //    ulong p_hi, p_lo;
    //    umul_ppmm(p_hi, p_lo, F->tab_w_pre[1][node], v);
    //    v = F->tab_w[1][node] * v - p_hi * F->mod;
    //    p[0] = u + v;
    //    p[1] = u + F->mod2 - v;
    //}
    //else if (order == 2)
    //{
    //    // in [0..2n), out [0..4n)
    //    ulong p_hi, p_lo;
    //    const ulong u0 = p[0];
    //    const ulong u1 = p[1];
    //    ulong v0 = p[2];
    //    ulong v1 = p[3];

    //    ulong w = F->tab_w[1][node];
    //    ulong wpre = F->tab_w_pre[1][node];
    //    umul_ppmm(p_hi, p_lo, wpre, v0);
    //    v0 = w * v0 - p_hi * F->mod;
    //    umul_ppmm(p_hi, p_lo, wpre, v1);
    //    v1 = w * v1 - p_hi * F->mod;
    //    ulong p0 = u0 + v0;  // [0..4n)
    //    if (p0 >= F->mod2)
    //        p0 -= F->mod2;      // [0..2n)
    //    ulong p1 = u1 + v1;  // [0..4n)
    //    v0 = u0 + F->mod2 - v0;  // [0..4n)
    //    if (v0 >= F->mod2)
    //        v0 -= F->mod2;       // [0..2n)
    //    v1 = u1 + F->mod2 - v1;  // [0..4n)

    //    w = F->tab_w[1][2*node];
    //    wpre = F->tab_w_pre[1][2*node];
    //    umul_ppmm(p_hi, p_lo, wpre, p1);
    //    p1 = w * p1 - p_hi * F->mod;
    //    w = F->tab_w[1][2*node+1];
    //    wpre = F->tab_w_pre[1][2*node+1];
    //    umul_ppmm(p_hi, p_lo, wpre, v1);
    //    v1 = w * v1 - p_hi * F->mod;
    //    p[0] = p0 + p1;             // [0..4n)
    //    p[1] = p0 + F->mod2 - p1;  // [0..4n)
    //    p[2] = v0 + v1;             // [0..4n)
    //    p[3] = v0 + F->mod2 - v1;  // [0..4n)
    //}
    //else
    if (order == 3)
    {
        // in [0..4n), out [0..4n)
        ulong p_hi, p_lo;

        ulong u0 = p[0];
        ulong u1 = p[1];
        ulong u2 = p[2];
        ulong u3 = p[3];
        ulong v0 = p[4];
        ulong v1 = p[5];
        ulong v2 = p[6];
        ulong v3 = p[7];

        // mod x**4 - w | x**4 + w
        ulong w = F->tab_w[1][2*node];
        ulong wpre = F->tab_w[1][2*node+1];
        umul_ppmm(p_hi, p_lo, wpre, v0);
        v0 = w * v0 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, v1);
        v1 = w * v1 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, v2);
        v2 = w * v2 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, v3);
        v3 = w * v3 - p_hi * F->mod;
        ulong p0 = u0 + v0;   // [0..6n)
        ulong p1 = u1 + v1;   // [0..6n)
        ulong p2 = u2 + v2;   // [0..6n)
        ulong p3 = u3 + v3;   // [0..6n)
        u0 += F->mod2 - v0;  // [0..6n)
        u1 += F->mod2 - v1;  // [0..6n)
        u2 += F->mod2 - v2;  // [0..6n)
        u3 += F->mod2 - v3;  // [0..6n)

        // left, mod x**2 - sqrt(w) | x**2 + sqrt(w)
        w = F->tab_w[1][4*node];
        wpre = F->tab_w[1][4*node+1];
        umul_ppmm(p_hi, p_lo, wpre, p2);
        p2 = w * p2 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, p3);
        p3 = w * p3 - p_hi * F->mod;
        v0 = p0 + p2;             // [0..8n)
        v1 = p1 + p3;             // [0..8n)
        v2 = p0 + F->mod2 - p2;  // [0..8n)
        v3 = p1 + F->mod2 - p3;  // [0..8n)

        // left-left, mod x - fort(w) | x + fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+1], v1);
        v1 = F->tab_w[1][8*node] * v1 - p_hi * F->mod;
        if (v0 >= F->mod4)
            v0 -= F->mod4;
        if (v0 >= F->mod2)
            v0 -= F->mod2;  // [0..2n)
        p[0] = v0 + v1;             // [0..4n)
        p[1] = v0 + F->mod2 - v1;  // [0..4n)

        // left-right, mod x - I*fort(w) | x+ I*fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+3], v3);
        v3 = F->tab_w[1][8*node+2] * v3 - p_hi * F->mod;
        if (v2 >= F->mod4)
            v2 -= F->mod4;
        if (v2 >= F->mod2)
            v2 -= F->mod2;  // [0..2n)
        p[2] = v2 + v3;              // [0..4n)
        p[3] = v2 + F->mod2 - v3;   // [0..4n)

        // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
        w = F->tab_w[1][4*node+2];
        wpre = F->tab_w[1][4*node+3];
        umul_ppmm(p_hi, p_lo, wpre, u2);
        u2 = w * u2 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, u3);
        u3 = w * u3 - p_hi * F->mod;
        v0 = u0 + u2;             // [0..8n)
        v1 = u1 + u3;             // [0..8n)
        v2 = u0 + F->mod2 - u2;  // [0..8n)
        v3 = u1 + F->mod2 - u3;  // [0..8n)

        // right-left, mod x - J*fort(w) | x + J*fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+5], v1);
        v1 = F->tab_w[1][8*node+4] * v1 - p_hi * F->mod;
        if (v0 >= F->mod4)
            v0 -= F->mod4;
        if (v0 >= F->mod2)
            v0 -= F->mod2;  // [0..2n)
        p[4] = v0 + v1;             // [0..4n)
        p[5] = v0 + F->mod2 - v1;  // [0..4n)

        // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+7], v3);
        v3 = F->tab_w[1][8*node+6] * v3 - p_hi * F->mod;
        if (v2 >= F->mod4)
            v2 -= F->mod4;
        if (v2 >= F->mod2)
            v2 -= F->mod2;  // [0..2n)
        p[6] = v2 + v3;             // [0..4n)
        p[7] = v2 + F->mod2 - v3;  // [0..4n)
    }
    else
    {
        // in: [0..4n), out: [0..4n)
        const ulong w = F->tab_w[1][2*node];
        const ulong wpre = F->tab_w[1][2*node+1];
        for (ulong k = 0; k < len/2; k+=4)
        {
            ulong p_hi, p_lo, u, v;
            u = p[k+0];
            if (u >= F->mod2)
                u -= F->mod2;  // [0..2n)
            v = p[len/2+k+0];
            umul_ppmm(p_hi, p_lo, wpre, v);
            v = w * v - p_hi * F->mod;
            p[k+0] = u + v;                   // [0..4n)
            p[len/2+k+0] = u + F->mod2 - v;  // [0..4n)

            u = p[k+1];
            if (u >= F->mod2)
                u -= F->mod2;  // [0..2n)
            v = p[len/2+k+1];
            umul_ppmm(p_hi, p_lo, wpre, v);
            v = w * v - p_hi * F->mod;
            p[k+1] = u + v;  // [0..4n)
            p[len/2+k+1] = u + F->mod2 - v;  // [0..4n)

            u = p[k+2];
            if (u >= F->mod2)
                u -= F->mod2;  // [0..2n)
            v = p[len/2+k+2];
            umul_ppmm(p_hi, p_lo, wpre, v);
            v = w * v - p_hi * F->mod;
            p[k+2] = u + v;  // [0..4n)
            p[len/2+k+2] = u + F->mod2 - v;  // [0..4n)

            u = p[k+3];
            if (u >= F->mod2)
                u -= F->mod2;  // [0..2n)
            v = p[len/2+k+3];
            umul_ppmm(p_hi, p_lo, wpre, v);
            v = w * v - p_hi * F->mod;
            p[k+3] = u + v;  // [0..4n)
            p[len/2+k+3] = u + F->mod2 - v;  // [0..4n)
        }
        _nmod_fft_red_rec2_lazy_general(p, len/2, order-1, 2*node, F);
        _nmod_fft_red_rec2_lazy_general(p+len/2, len/2, order-1, 2*node+1, F);
    }
}

void _nmod_fft_red_iter2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // perform FFT layers up to order 3
    ulong llen = len;
    for (ulong ell = order; ell > 3; ell--, llen>>=1)
    {
        // k = 0: point is 1, handle separately
        ulong tmp;
        for (ulong kk = 0; kk < llen/2; kk++)
            // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
            DFT2_LAZY2_RED1(p[kk], p[llen/2+kk], F->mod2, tmp);

        ulong node = 1;  // index of current point in tab
        for (ulong k = llen; k < len; k+=llen, node++)
        {
            ulong w = F->tab_w[1][2*node];
            ulong wpre = F->tab_w[1][2*node+1];

            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                // in: p[k+kk], p[llen/2+k+kk] in [0..4n)
                // out: p[k+kk], p[llen/2+k+kk] in [0..4n)
                ulong p_hi, p_lo, u, v;
                u = p[k+kk+0];
                if (u >= F->mod2)
                    u -= F->mod2;  // [0..2n)
                v = p[llen/2+k+kk+0];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod;
                p[k+kk+0] = u + v;                   // [0..4n)
                p[llen/2+k+kk+0] = u + F->mod2 - v;  // [0..4n)

                u = p[k+kk+1];
                if (u >= F->mod2)
                    u -= F->mod2;  // [0..2n)
                v = p[llen/2+k+kk+1];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod;
                p[k+kk+1] = u + v;  // [0..4n)
                p[llen/2+k+kk+1] = u + F->mod2 - v;  // [0..4n)

                u = p[k+kk+2];
                if (u >= F->mod2)
                    u -= F->mod2;  // [0..2n)
                v = p[llen/2+k+kk+2];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod;
                p[k+kk+2] = u + v;  // [0..4n)
                p[llen/2+k+kk+2] = u + F->mod2 - v;  // [0..4n)

                u = p[k+kk+3];
                if (u >= F->mod2)
                    u -= F->mod2;  // [0..2n)
                v = p[llen/2+k+kk+3];
                umul_ppmm(p_hi, p_lo, wpre, v);
                v = w * v - p_hi * F->mod;
                p[k+kk+3] = u + v;  // [0..4n)
                p[llen/2+k+kk+3] = u + F->mod2 - v;  // [0..4n)
            }
        }
    }

    // perform last two FFT layers
    dft8_red_lazy(p, F);  // k == 0..7, in [0..2n) (is ok!), out [0..4n)
    ulong node = 1;
    for (ulong k = 8; k < len; k+=8, node++)
    {
        // in [0..4n), out [0..4n)
        ulong p_hi, p_lo;

        ulong u0 = p[k+0];
        ulong u1 = p[k+1];
        ulong u2 = p[k+2];
        ulong u3 = p[k+3];
        ulong v0 = p[k+4];
        ulong v1 = p[k+5];
        ulong v2 = p[k+6];
        ulong v3 = p[k+7];

        // mod x**4 - w | x**4 + w
        ulong w = F->tab_w[1][2*node];
        ulong wpre = F->tab_w[1][2*node+1];
        umul_ppmm(p_hi, p_lo, wpre, v0);
        v0 = w * v0 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, v1);
        v1 = w * v1 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, v2);
        v2 = w * v2 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, v3);
        v3 = w * v3 - p_hi * F->mod;
        ulong p0 = u0 + v0;   // [0..6n)
        ulong p1 = u1 + v1;   // [0..6n)
        ulong p2 = u2 + v2;   // [0..6n)
        ulong p3 = u3 + v3;   // [0..6n)
        u0 += F->mod2 - v0;  // [0..6n)
        u1 += F->mod2 - v1;  // [0..6n)
        u2 += F->mod2 - v2;  // [0..6n)
        u3 += F->mod2 - v3;  // [0..6n)

        // left, mod x**2 - sqrt(w) | x**2 + sqrt(w)
        w = F->tab_w[1][4*node];
        wpre = F->tab_w[1][4*node+1];
        umul_ppmm(p_hi, p_lo, wpre, p2);
        p2 = w * p2 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, p3);
        p3 = w * p3 - p_hi * F->mod;
        v0 = p0 + p2;             // [0..8n)
        v1 = p1 + p3;             // [0..8n)
        v2 = p0 + F->mod2 - p2;  // [0..8n)
        v3 = p1 + F->mod2 - p3;  // [0..8n)

        // left-left, mod x - fort(w) | x + fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+1], v1);
        v1 = F->tab_w[1][8*node] * v1 - p_hi * F->mod;
        if (v0 >= F->mod4)
            v0 -= F->mod4;
        if (v0 >= F->mod2)
            v0 -= F->mod2;  // [0..2n)
        p[k+0] = v0 + v1;             // [0..4n)
        p[k+1] = v0 + F->mod2 - v1;  // [0..4n)

        // left-right, mod x - I*fort(w) | x+ I*fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+3], v3);
        v3 = F->tab_w[1][8*node+2] * v3 - p_hi * F->mod;
        if (v2 >= F->mod4)
            v2 -= F->mod4;
        if (v2 >= F->mod2)
            v2 -= F->mod2;  // [0..2n)
        p[k+2] = v2 + v3;              // [0..4n)
        p[k+3] = v2 + F->mod2 - v3;   // [0..4n)

        // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
        w = F->tab_w[1][4*node+2];
        wpre = F->tab_w[1][4*node+3];
        umul_ppmm(p_hi, p_lo, wpre, u2);
        u2 = w * u2 - p_hi * F->mod;
        umul_ppmm(p_hi, p_lo, wpre, u3);
        u3 = w * u3 - p_hi * F->mod;
        v0 = u0 + u2;             // [0..8n)
        v1 = u1 + u3;             // [0..8n)
        v2 = u0 + F->mod2 - u2;  // [0..8n)
        v3 = u1 + F->mod2 - u3;  // [0..8n)

        // right-left, mod x - J*fort(w) | x + J*fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+5], v1);
        v1 = F->tab_w[1][8*node+4] * v1 - p_hi * F->mod;
        if (v0 >= F->mod4)
            v0 -= F->mod4;
        if (v0 >= F->mod2)
            v0 -= F->mod2;  // [0..2n)
        p[k+4] = v0 + v1;             // [0..4n)
        p[k+5] = v0 + F->mod2 - v1;  // [0..4n)

        // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
        umul_ppmm(p_hi, p_lo, F->tab_w[1][8*node+7], v3);
        v3 = F->tab_w[1][8*node+6] * v3 - p_hi * F->mod;
        if (v2 >= F->mod4)
            v2 -= F->mod4;
        if (v2 >= F->mod2)
            v2 -= F->mod2;  // [0..2n)
        p[k+6] = v2 + v3;             // [0..4n)
        p[k+7] = v2 + F->mod2 - v3;  // [0..4n)
    }
    //for (ulong k = 0; k < len; k+=4)
    //    DFT4_DIF_SHOUP_LAZY2_RED(p[k+0], p[k+1], p[k+2], p[k+3], F->tab_w[0][1], F->tab_w_pre[0][1], F->mod, F->mod2);
}

// in [0..2n) out [0..4n)
void _nmod_fft_dif_rec4_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        // in [0..2n) out [0..4n)
        DFT2_LAZY3_RED(p[0], p[1], F->mod2);
    else if (order == 2)
        // in [0..2n) out [0..4n)
        DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->Ipre, F->mod, F->mod2);
    else if (order == 3)
        dft8_red_lazy(p, F);
    else
    {
        // in [0..2n) out [0..2n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+(len/4);
        const nn_ptr p2 = p1+(len/4);
        const nn_ptr p3 = p2+(len/4);
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k++)
        {
            BUTTERFLY4_GS_LAZY(p0[k], p1[k], p2[k], p3[k],
                               F->I, F->Ipre, F->tab_w[order-4][2*k], F->tab_w[order-4][2*k+1], F->tab_w[order-4][4*k], F->tab_w[order-4][4*k+1], F->tab_w[order-4][6*k], F->tab_w[order-4][6*k+1],
                               F->mod, F->mod2, F->mod4, p_hi, p_lo);
        }
        _nmod_fft_dif_rec4_lazy(p0, len/4, order-2, F);
        _nmod_fft_dif_rec4_lazy(p1, len/4, order-2, F);
        _nmod_fft_dif_rec4_lazy(p2, len/4, order-2, F);
        _nmod_fft_dif_rec4_lazy(p3, len/4, order-2, F);
    }
}

// in [0..2n) out [0..4n)
void _nmod_fft_dif_rec8_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        // in [0..2n) out [0..4n)
        DFT2_LAZY3_RED(p[0], p[1], F->mod2);
    else if (order == 2)
        // in [0..2n) out [0..4n)
        DFT4_DIF_SHOUP_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->Ipre, F->mod, F->mod2);
    else if (order == 3)
        dft8_red_lazy(p, F);
    else
    {
        // in [0..2n) out [0..2n)
        const nn_ptr pp0 = p;
        const nn_ptr pp1 = pp0+(len/8);
        const nn_ptr pp2 = pp1+(len/8);
        const nn_ptr pp3 = pp2+(len/8);
        const nn_ptr pp4 = pp3+(len/8);
        const nn_ptr pp5 = pp4+(len/8);
        const nn_ptr pp6 = pp5+(len/8);
        const nn_ptr pp7 = pp6+(len/8);
        for (ulong k = 0; k < len/8; k++)
        {
            ulong p_hi, p_lo;
            ulong u0 = pp0[k];
            ulong u1 = pp1[k];
            ulong u2 = pp2[k];
            ulong u3 = pp3[k];
            ulong v0 = pp4[k];
            ulong v1 = pp5[k];
            ulong v2 = pp6[k];
            ulong v3 = pp7[k];

            // mod x**4 - 1 | x**4 + 1
            ulong p0 = u0 + v0;  // [0..4n)
            ulong p1 = u1 + v1;  // [0..4n)
            ulong p2 = u2 + v2;  // [0..4n)
            ulong p3 = u3 + v3;  // [0..4n)
            u0 += F->mod2 - v0;  // [0..4n)
            u1 += F->mod2 - v1;  // [0..4n)
            u2 += F->mod2 - v2;  // [0..4n)
            u3 += F->mod2 - v3;  // [0..4n)

            // left, mod x**2 - 1 | x**2 + 1
            v0 = p0 + p2;             // [0..8n)
            v1 = p1 + p3;             // [0..8n)
            v2 = p0 + F->mod4 - p2;  // [0..8n)
            v3 = p1 + F->mod4 - p3;  // [0..8n)

            // left-left, mod x-1 | x+1
            p0 = v0 + v1;               // [0..16n) // TODO 16n
            p1 = v0 + 2*F->mod4 - v1;  // [0..16n) // TODO 16n
            if (p0 > 2*F->mod4)
                p0 -= 2*F->mod4;
            if (p0 > F->mod4)
                p0 -= F->mod4;        // [0..4n)
            if (p0 > F->mod2)
                p0 -= F->mod2;        // [0..2n)
            pp0[k] = p0;
            // scale p1 <- w**(4*k) * p1
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][8*k+1], p1);
            pp1[k] = F->tab_w[order-4][8*k] * p1 - p_hi * F->mod;  // [0..2n)

            // left-right, mod x-I | x+I
            umul_ppmm(p_hi, p_lo, F->Ipre, v3);
            v3 = F->I * v3 - p_hi * F->mod;  // [0..2n)
            p2 = v2 + v3;             // [0..10n)
            p3 = v2 + F->mod2 - v3;  // [0..10n)
            // scale p2 <-- w**(2*k) * p2, p3 <-- w**(6*k) * p3
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][4*k+1], p2);
            pp2[k] = F->tab_w[order-4][4*k] * p2 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][12*k+1], p3);
            pp3[k] = F->tab_w[order-4][12*k] * p3 - p_hi * F->mod;  // [0..2n)

            // right, mod x**2 - I | x**2 + I
            umul_ppmm(p_hi, p_lo, F->Ipre, u2);
            u2 = F->I * u2 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->Ipre, u3);
            u3 = F->I * u3 - p_hi * F->mod;  // [0..2n)
            v0 = u0 + u2;  // [0..6n)
            v1 = u1 + u3;  // [0..6n)
            v2 = u0 + F->mod2 - u2;  // [0..6n)
            v3 = u1 + F->mod2 - u3;  // [0..6n)

            // right-left, mod x - J | x + J
            umul_ppmm(p_hi, p_lo, F->Jpre, v1);
            v1 = F->J * v1 - p_hi * F->mod;  // [0..2n)
            p0 = v0 + v1;  // [0..8n)
            p1 = v0 + F->mod2 - v1;  // [0..8n)
            // scale p0 <-- w**k * p0, p1 <-- w**(5*k) * p1
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][2*k+1], p0);
            pp4[k] = F->tab_w[order-4][2*k] * p0 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][10*k+1], p1);
            pp5[k] = F->tab_w[order-4][10*k] * p1 - p_hi * F->mod;  // [0..2n)

            // right-right, mod x - I*J | x + I*J
            umul_ppmm(p_hi, p_lo, F->IJpre, v3);
            v3 = F->IJ * v3 - p_hi * F->mod;  // [0..2n)
            p2 = v2 + v3;   // [0..8n)
            p3 = v2 + F->mod2 - v3;   // [0..8n)
            // scale p2 <-- w**(3*k) * p2, p3 <-- w**(7*k) * p3
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][6*k+1], p2);
            pp6[k] = F->tab_w[order-4][6*k] * p2 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->tab_w[order-4][14*k+1], p3);
            pp7[k] = F->tab_w[order-4][14*k] * p3 - p_hi * F->mod;  // [0..2n)
        }
        _nmod_fft_dif_rec8_lazy(pp0, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp1, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp2, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp3, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp4, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp5, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp6, len/8, order-3, F);
        _nmod_fft_dif_rec8_lazy(pp7, len/8, order-3, F);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
