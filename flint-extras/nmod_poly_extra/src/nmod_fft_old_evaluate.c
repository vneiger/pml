#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_vec.h>
#include <immintrin.h>

#include "nmod_poly_fft.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COPIED FROM n_fft_EVALUATE.c */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

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

/*-----------------*/
/* 2-point macros  */
/*-----------------*/

/** 2-point DFT: in-place transform
 *                    [1   1]
 * [a  b]  <-  [a  b] [1  -1]
 * */

// lazy2 red1:
// in [0..2n) x [0..2n), out [0..2n) x [0..4n)
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
// in [0..2n) x [0..2n), out [0..2n) x [0..4n)
// b unchanged, a-b stored in tmp
// n2 is 2*n, tmp is a temporary (ulong)
#define DFT2_LAZY2_RED1_OUT(a, b, n2, tmp) \
    do {                                   \
        tmp = (a) + (n2) - (b);            \
        (a) = (a) + (b);                   \
        if ((a) >= (n2))                   \
            (a) -= (n2);                   \
    } while(0)

// lazy4 red:
// in/out: [0..4n)
// n4 is 4*n
// currently only used in base cases
#define DFT2_LAZY4_RED(a, b, n4) \
    do {                         \
        const ulong tmp = (a);   \
        (a) = (a) + (b);         \
        if ((a) >= (n4))         \
            (a) -= (n4);         \
        (b) = tmp + n4 - (b);    \
        if ((b) >= (n4))         \
            (b) -= (n4);         \
    } while(0)

/** 2-point butterflies
 *
 * Gentleman-Sande: in-place transform
 *                    [1   w]
 * [a  b]  <-  [a  b] [1  -w]
 *
 * Cooley-Tukey: in-place transform
 *                    [1   1]
 * [a  b]  <-  [a  b] [w  -w]
 * */

// Gentleman-Sande butterfly, in-place, lazy
// in/out:  a,b in [0..2n) x [0..2n)
// n2 is 2*n, w_pr is the precomputed data for multiplication by w mod n
// p_hi, p_lo, tmp are temporaries
#define BUTTERFLY2_GS_LAZY(a, b, n, n2, w, w_pr, p_hi, p_lo, tmp) \
    do {                                                         \
        DFT2_LAZY2_RED1_OUT(a, b, n2, tmp);                      \
        N_MULMOD_PRECOMP_LAZY(b, w, tmp, w_pr, n, p_hi, p_lo);   \
    } while(0)

// Cooley-Tukey butterfly, in-place, lazy
// in/out:  a,b in [0..4n) x [0..4n)
// n2 is 2*n, w_pr is the precomputed data for multiplication by w mod n
// p_hi, p_lo, u, v are temporaries
#define BUTTERFLY2_CT_LAZY(a, b, n, n2, w, w_pr, p_hi, p_lo, u, v) \
    do {                                                          \
        u = (a);                                                  \
        if (u >= n2)                                              \
            u -= n2;  /* [0..2n) */                               \
        v = (b);                                                  \
        N_MULMOD_PRECOMP_LAZY(v, w, v, w_pr, n, p_hi, p_lo);      \
        (a) = u + v;                   /* [0..4n) */              \
        (b) = u + F->mod2 - v;         /* [0..4n) */              \
    } while(0)

/*----------------*/
/* 4-point macros */
/*----------------*/

/** 4-point DFT: in-place transform
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 * I is typically a square root of -1 (this property is not exploited)
 * */

// lazy2 red: input in [0..2*n) --> output [0..4*n)
// n2 is 2*n
#define DFT4_LAZY2_RED(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo)                     \
    do {                                                               \
        const ulong v0 = (a);                                          \
        const ulong v1 = (b);                                          \
        const ulong v2 = (c);                                          \
        const ulong v3 = (d);                                          \
        ulong v4 = v0 + v2;                         /* < 4*n */        \
        if (v4 >= (n2))                                                \
            v4 -= (n2);                             /* < 2*n */        \
        ulong v5 = v0 + (n2) - v2;                  /* < 4*n */        \
        if (v5 >= (n2))                                                \
            v5 -= (n2);                             /* < 2*n */        \
        ulong v6 = v1 + v3;                         /* < 4*n */        \
        if (v6 >= (n2))                                                \
            v6 -= (n2);                             /* < 2*n */        \
        ulong v7;                                                      \
        N_MULMOD_PRECOMP_LAZY(v7, (I), v1 + (n2) - v3, (I_pr), (n),    \
                              p_hi, p_lo);                             \
        (a) = v4 + v6;                              /* < 4*n */        \
        (b) = v4 + (n2) - v6;                       /* < 4*n */        \
        (c) = v5 + v7;                              /* < 4*n */        \
        (d) = v5 + (n2) - v7;                       /* < 4*n */        \
    } while(0)


/** 4-point butterflies
 *
 * Gentleman-Sande: in-place transform
 *                              [1  w2    w1    w3]
 *                              [1 -w2  I*w1 -I*w3]
 * [a  b  c  d] <- [a  b  c  d] [1  w2   -w1   -w3]
 *                              [1 -w2 -I*w1  I*w3]
 *
 * Cooley-Tukey: in-place transform   TODO
 *                    [1   1]
 * [a  b]  <-  [a  b] [w  -w]
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 * */

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


// in [0..4n), out [0..4n)
#define BUTTERFLY4_CT_LAZY(a,b,c,d,                                    \
                           w1,w1_pr,w2,w2_pr,w3,w3_pr,                 \
                           n,n2,p_hi,p_lo,tmp)                         \
    do {                                                               \
            ulong u0 = (a);                                            \
            ulong u1 = (b);                                            \
            ulong u2 = (c);                                            \
            ulong u3 = (d);                                            \
            if (u0 >= n2)                                              \
                u0 -= n2;                                              \
            if (u1 >= n2)                                              \
                u1 -= n2;                                              \
                                                                       \
            N_MULMOD_PRECOMP_LAZY(u2, w1, u2, w1_pr, n, p_hi, p_lo);   \
            tmp = u0;                                                  \
            u0 = u0 + u2;                    /* [0..4n) */             \
            u2 = tmp + n2 - u2;              /* [0..4n) */             \
            if (u0 >= n2)                                              \
                u0 -= n2;                    /* [0..2n) */             \
            if (u2 >= n2)                                              \
                u2 -= n2;                    /* [0..2n) */             \
                                                                       \
            N_MULMOD_PRECOMP_LAZY(u3, w1, u3, w1_pr, n, p_hi, p_lo);   \
            tmp = u1;                                                  \
            u1 = u1 + u3;                    /* [0..8n) */             \
            u3 = tmp + n2 - u3;              /* [0..8n) */             \
                                                                       \
            N_MULMOD_PRECOMP_LAZY(u1, w2, u1, w2_pr, n, p_hi, p_lo);   \
            tmp = u0;                                                  \
            (a) = u0 + u1;                    /* [0..4n) */             \
            (b) = tmp + n2 - u1;              /* [0..4n) */             \
                                                                       \
            N_MULMOD_PRECOMP_LAZY(u3, w3, u3, w3_pr, n, p_hi, p_lo);   \
            tmp = u2;                                                  \
            (c) = u2 + u3;                    /* [0..4n) */             \
            (d) = tmp + n2 - u3;              /* [0..4n) */             \
    } while(0)

/*----------------*/
/* 8-point DFT    */
/*----------------*/

// reduction tree 8-point lazy DFT
// lazy red: input in [0..2*n) --> output in [0..4*n)
// max value < 8n
FLINT_FORCE_INLINE void dft8_red_lazy(nn_ptr p, n_fft_old_ctx_t F)
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
    if (v0 >= F->mod4)
        v0 -= F->mod4;
    if (v1 >= F->mod4)
        v1 -= F->mod4;
    p0 = v0 + v1;               // [0..8n)
    p1 = v0 + F->mod4 - v1;     // [0..8n)
    if (p0 >= F->mod4)
        p0 -= F->mod4;
    if (p1 >= F->mod4)
        p1 -= F->mod4;
    p[0] = p0;
    p[1] = p1;

    // left-right, mod x-I | x+I
    N_MULMOD_PRECOMP_LAZY(v3, F->I, v3, F->I_pr, F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;         // [0..2n)
    p[2] = v2 + v3;             // [0..4n)
    p[3] = v2 + F->mod2 - v3;  // [0..4n)

    // right, mod x**2 - I | x**2 + I
    N_MULMOD_PRECOMP_LAZY(u2, F->I, u2, F->I_pr, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(u3, F->I, u3, F->I_pr, F->mod, p_hi, p_lo);
    if (u0 >= F->mod2)
        u0 -= F->mod2;         // [0..2n)
    if (u1 >= F->mod2)
        u1 -= F->mod2;         // [0..2n)
    v0 = u0 + u2;  // [0..4n)
    v1 = u1 + u3;  // [0..4n)
    v2 = u0 + F->mod2 - u2;  // [0..4n)
    v3 = u1 + F->mod2 - u3;  // [0..4n)

    // right-left, mod x - J | x + J
    N_MULMOD_PRECOMP_LAZY(v1, F->J, v1, F->J_pr, F->mod, p_hi, p_lo);
    if (v0 >= F->mod2)
        v0 -= F->mod2;         // [0..2n)
    p[4] = v0 + v1;
    p[5] = v0 + F->mod2 - v1;

    // right-right, mod x - I*J | x + I*J
    N_MULMOD_PRECOMP_LAZY(v3, F->IJ, v3, F->IJ_pr, F->mod, p_hi, p_lo);
    if (v2 >= F->mod2)
        v2 -= F->mod2;         // [0..2n)
    p[6] = v2 + v3;
    p[7] = v2 + F->mod2 - v3;
}

// in [0..4n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft8_red_lazy_general(nn_ptr p, ulong node, n_fft_old_ctx_t F)
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

    // mod x**4 - w | x**4 + w
    ulong w = F->tab_w[1][2*node];
    ulong wpre = F->tab_w[1][2*node+1];
    N_MULMOD_PRECOMP_LAZY(v0, w, v0, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(v1, w, v1, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(v2, w, v2, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(v3, w, v3, wpre, F->mod, p_hi, p_lo);
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
    N_MULMOD_PRECOMP_LAZY(p2, w, p2, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(p3, w, p3, wpre, F->mod, p_hi, p_lo);
    v0 = p0 + p2;             // [0..8n)
    v1 = p1 + p3;             // [0..8n)
    v2 = p0 + F->mod2 - p2;  // [0..8n)
    v3 = p1 + F->mod2 - p3;  // [0..8n)

    // left-left, mod x - fort(w) | x + fort(w)
    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[1][8*node], v1, F->tab_w[1][8*node+1], F->mod, p_hi, p_lo);
    if (v0 >= F->mod4)
        v0 -= F->mod4;
    if (v0 >= F->mod2)
        v0 -= F->mod2;  // [0..2n)
    p[0] = v0 + v1;             // [0..4n)
    p[1] = v0 + F->mod2 - v1;  // [0..4n)

    // left-right, mod x - I*fort(w) | x+ I*fort(w)
    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[1][8*node+2], v3, F->tab_w[1][8*node+3], F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;  // [0..2n)
    p[2] = v2 + v3;              // [0..4n)
    p[3] = v2 + F->mod2 - v3;   // [0..4n)

    // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
    w = F->tab_w[1][4*node+2];
    wpre = F->tab_w[1][4*node+3];
    N_MULMOD_PRECOMP_LAZY(u2, w, u2, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(u3, w, u3, wpre, F->mod, p_hi, p_lo);
    v0 = u0 + u2;             // [0..8n)
    v1 = u1 + u3;             // [0..8n)
    v2 = u0 + F->mod2 - u2;  // [0..8n)
    v3 = u1 + F->mod2 - u3;  // [0..8n)

    // right-left, mod x - J*fort(w) | x + J*fort(w)
    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[1][8*node+4], v1, F->tab_w[1][8*node+5], F->mod, p_hi, p_lo);
    if (v0 >= F->mod4)
        v0 -= F->mod4;
    if (v0 >= F->mod2)
        v0 -= F->mod2;  // [0..2n)
    p[4] = v0 + v1;             // [0..4n)
    p[5] = v0 + F->mod2 - v1;  // [0..4n)

    // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[1][8*node+6], v3, F->tab_w[1][8*node+7], F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;  // [0..2n)
    p[6] = v2 + v3;             // [0..4n)
    p[7] = v2 + F->mod2 - v3;  // [0..4n)
}
/* ALTERNATIVE to dft8 general: */
     //   ulong w2 = F->tab_w[1][2*node];
     //   ulong w2pre = F->tab_w[1][2*node+1];
     //   ulong w = F->tab_w[1][4*node];
     //   ulong wpre = F->tab_w[1][4*node+1];
     //   ulong Iw = F->tab_w[1][4*node+2];
     //   ulong Iwpre = F->tab_w[1][4*node+3];

     //   ulong p_hi, p_lo, u, v;

     //   BUTTERFLY4_CT_LAZY(p[0], p[2], p[4], p[6], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, u);
     //   BUTTERFLY4_CT_LAZY(p[1], p[3], p[5], p[7], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, u);

     //   BUTTERFLY2_CT_LAZY(p[0], p[1], F->mod, F->mod2, F->tab_w[1][8*node], F->tab_w[1][8*node+1], p_hi, p_lo, u, v);
     //   BUTTERFLY2_CT_LAZY(p[2], p[3], F->mod, F->mod2, F->tab_w[1][8*node+2], F->tab_w[1][8*node+3], p_hi, p_lo, u, v);
     //   BUTTERFLY2_CT_LAZY(p[4], p[5], F->mod, F->mod2, F->tab_w[1][8*node+4], F->tab_w[1][8*node+5], p_hi, p_lo, u, v);
     //   BUTTERFLY2_CT_LAZY(p[6], p[7], F->mod, F->mod2, F->tab_w[1][8*node+6], F->tab_w[1][8*node+7], p_hi, p_lo, u, v);


/*------------------*/
/* other base cases */
/*------------------*/

// in [0..2n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft16_red_lazy(nn_ptr p, n_fft_old_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT4_LAZY2_RED(p[0], p[4], p[8], p[12], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_LAZY2_RED(p[1], p[5], p[9 ], p[13], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_LAZY2_RED(p[2], p[6], p[10], p[14], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_LAZY2_RED(p[3], p[7], p[11], p[15], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;

    // next line requires < 2n, hence the four reductions above
    DFT4_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    BUTTERFLY4_CT_LAZY(p[4], p[5], p[6], p[7], F->I, F->I_pr, F->J, F->J_pr, F->IJ, F->IJ_pr, F->mod, F->mod2, p_hi, p_lo, tmp);
    BUTTERFLY4_CT_LAZY(p[8], p[9], p[10], p[11], F->J, F->J_pr, F->tab_w[1][8], F->tab_w[1][9], F->tab_w[1][10], F->tab_w[1][11], F->mod, F->mod2, p_hi, p_lo, tmp);
    BUTTERFLY4_CT_LAZY(p[12], p[13], p[14], p[15], F->IJ, F->IJ_pr, F->tab_w[1][12], F->tab_w[1][13], F->tab_w[1][14], F->tab_w[1][15], F->mod, F->mod2, p_hi, p_lo, tmp);
}

// in [0..2n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft32_red_lazy(nn_ptr p, n_fft_old_ctx_t F)
{
    ulong p_hi, p_lo;

    DFT4_LAZY2_RED(p[0], p[8 ], p[16], p[24], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_LAZY2_RED(p[1], p[9 ], p[17], p[25], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_LAZY2_RED(p[2], p[10], p[18], p[26], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_LAZY2_RED(p[3], p[11], p[19], p[27], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;
    DFT4_LAZY2_RED(p[4], p[12], p[20], p[28], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[4] >= F->mod2)
        p[4] -= F->mod2;
    DFT4_LAZY2_RED(p[5], p[13], p[21], p[29], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[5] >= F->mod2)
        p[5] -= F->mod2;
    DFT4_LAZY2_RED(p[6], p[14], p[22], p[30], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[6] >= F->mod2)
        p[6] -= F->mod2;
    DFT4_LAZY2_RED(p[7], p[15], p[23], p[31], F->I, F->I_pr, F->mod, F->mod2, p_hi, p_lo);
    if (p[7] >= F->mod2)
        p[7] -= F->mod2;

    // next line requires < 2n, hence the 8 reductions above
    dft8_red_lazy(p, F);
    dft8_red_lazy_general(p+8, 1, F);
    dft8_red_lazy_general(p+16, 2, F);
    dft8_red_lazy_general(p+24, 3, F);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* <---- END COPIED FROM n_fft_EVALUATE.c */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/







// input [0..2*n), output [0..4*n)
// depth >= 3
void _n_fft_old_dif_rec2_lazy(nn_ptr p, ulong len, ulong depth, n_fft_old_ctx_t F)
{
    // depth == 0: nothing to do
    //if (depth == 1)
    //    DFT2_LAZY4_RED(p[0], p[1], F->mod4);
    //else if (depth == 2)
    //    DFT4_LAZY2_RED(p[0], p[1], p[2], p[3], F->I, F->I_pr, F->mod, F->mod2);
    //else
    if (depth == 3)
        dft8_red_lazy(p, F);
    else
    {
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        const nn_ptr ww = F->tab_w[depth-4];
        // in: p0, p1 in [0..2n); out: p0, p1 in [0..2n)
        for (ulong k = 0; k < len/2; k+=4)
        {
            ulong p_hi, p_lo, tmp;
            BUTTERFLY2_GS_LAZY(p0[k+0], p1[k+0], F->mod, F->mod2, ww[2*k+0], ww[2*k+1], p_hi, p_lo, tmp);
            BUTTERFLY2_GS_LAZY(p0[k+1], p1[k+1], F->mod, F->mod2, ww[2*k+2], ww[2*k+3], p_hi, p_lo, tmp);
            BUTTERFLY2_GS_LAZY(p0[k+2], p1[k+2], F->mod, F->mod2, ww[2*k+4], ww[2*k+5], p_hi, p_lo, tmp);
            BUTTERFLY2_GS_LAZY(p0[k+3], p1[k+3], F->mod, F->mod2, ww[2*k+6], ww[2*k+7], p_hi, p_lo, tmp);
        }
        _n_fft_old_dif_rec2_lazy(p0, len/2, depth-1, F);
        _n_fft_old_dif_rec2_lazy(p1, len/2, depth-1, F);
    }
}

// depth >= 3
void _n_fft_old_dif_iter2_lazy(nn_ptr p, ulong len, ulong depth, n_fft_old_ctx_t F)
{
    // perform FFT layers up to depth 3
    ulong llen = len;
    for (ulong ell = depth; ell > 3; ell--, llen>>=1)
        for (ulong k = 0; k < len; k+=llen)
        {
            const nn_ptr p0 = p+k;
            const nn_ptr p1 = p+llen/2+k;
            const nn_ptr ww = F->tab_w[ell-4];
            ulong p_hi, p_lo, tmp;
            // in: p0, p1 in [0..2n); out: p0, p1 in [0..2n)
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                BUTTERFLY2_GS_LAZY(p0[kk+0], p1[kk+0], F->mod, F->mod2, ww[2*kk+0], ww[2*kk+1], p_hi, p_lo, tmp);
                BUTTERFLY2_GS_LAZY(p0[kk+1], p1[kk+1], F->mod, F->mod2, ww[2*kk+2], ww[2*kk+3], p_hi, p_lo, tmp);
                BUTTERFLY2_GS_LAZY(p0[kk+2], p1[kk+2], F->mod, F->mod2, ww[2*kk+4], ww[2*kk+5], p_hi, p_lo, tmp);
                BUTTERFLY2_GS_LAZY(p0[kk+3], p1[kk+3], F->mod, F->mod2, ww[2*kk+6], ww[2*kk+7], p_hi, p_lo, tmp);
            }
        }

    // perform last FFT layers
    for (ulong k = 0; k < len; k+=8)
        dft8_red_lazy(p+k, F);
}

// in [0..2n) out [0..4n)
void _n_fft_old_dif_rec4_lazy(nn_ptr p, ulong len, ulong depth, n_fft_old_ctx_t F)
{
    if (depth == 3)
        dft8_red_lazy(p, F);
    else if (depth == 4)
    {
        ulong p_hi, p_lo, tmp;
        DFT2_LAZY2_RED1(p[0], p[0+8], F->mod2, tmp);
        BUTTERFLY2_GS_LAZY(p[2], p[2+8], F->mod, F->mod2, F->J, F->J_pr, p_hi, p_lo, tmp);
        BUTTERFLY2_GS_LAZY(p[4], p[4+8], F->mod, F->mod2, F->I, F->I_pr, p_hi, p_lo, tmp);
        BUTTERFLY2_GS_LAZY(p[6], p[6+8], F->mod, F->mod2, F->IJ, F->IJ_pr, p_hi, p_lo, tmp);
        BUTTERFLY2_GS_LAZY(p[1], p[1+8], F->mod, F->mod2, F->tab_w[0][2*1], F->tab_w[0][2*1+1], p_hi, p_lo, tmp);
        BUTTERFLY2_GS_LAZY(p[3], p[3+8], F->mod, F->mod2, F->tab_w[0][2*3], F->tab_w[0][2*3+1], p_hi, p_lo, tmp);
        BUTTERFLY2_GS_LAZY(p[5], p[5+8], F->mod, F->mod2, F->tab_w[0][2*5], F->tab_w[0][2*5+1], p_hi, p_lo, tmp);
        BUTTERFLY2_GS_LAZY(p[7], p[7+8], F->mod, F->mod2, F->tab_w[0][2*7], F->tab_w[0][2*7+1], p_hi, p_lo, tmp);
        dft8_red_lazy(p, F);
        dft8_red_lazy(p+8, F);
    }
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
                               F->I, F->I_pr, F->tab_w[depth-4][2*k], F->tab_w[depth-4][2*k+1], F->tab_w[depth-4][4*k], F->tab_w[depth-4][4*k+1], F->tab_w[depth-4][6*k], F->tab_w[depth-4][6*k+1],
                               F->mod, F->mod2, F->mod4, p_hi, p_lo);
        }
        _n_fft_old_dif_rec4_lazy(p0, len/4, depth-2, F);
        _n_fft_old_dif_rec4_lazy(p1, len/4, depth-2, F);
        _n_fft_old_dif_rec4_lazy(p2, len/4, depth-2, F);
        _n_fft_old_dif_rec4_lazy(p3, len/4, depth-2, F);
    }
}

// in [0..2n) out [0..4n)
void _n_fft_old_dif_rec8_lazy(nn_ptr p, ulong len, ulong depth, n_fft_old_ctx_t F)
{
    if (depth == 3)
        dft8_red_lazy(p, F);
    else if (depth <= 5)
        _n_fft_old_dif_rec2_lazy(p, len, depth, F);
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
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][8*k+1], p1);
            pp1[k] = F->tab_w[depth-4][8*k] * p1 - p_hi * F->mod;  // [0..2n)

            // left-right, mod x-I | x+I
            umul_ppmm(p_hi, p_lo, F->I_pr, v3);
            v3 = F->I * v3 - p_hi * F->mod;  // [0..2n)
            p2 = v2 + v3;             // [0..10n)
            p3 = v2 + F->mod2 - v3;  // [0..10n)
            // scale p2 <-- w**(2*k) * p2, p3 <-- w**(6*k) * p3
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][4*k+1], p2);
            pp2[k] = F->tab_w[depth-4][4*k] * p2 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][12*k+1], p3);
            pp3[k] = F->tab_w[depth-4][12*k] * p3 - p_hi * F->mod;  // [0..2n)

            // right, mod x**2 - I | x**2 + I
            umul_ppmm(p_hi, p_lo, F->I_pr, u2);
            u2 = F->I * u2 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->I_pr, u3);
            u3 = F->I * u3 - p_hi * F->mod;  // [0..2n)
            v0 = u0 + u2;  // [0..6n)
            v1 = u1 + u3;  // [0..6n)
            v2 = u0 + F->mod2 - u2;  // [0..6n)
            v3 = u1 + F->mod2 - u3;  // [0..6n)

            // right-left, mod x - J | x + J
            umul_ppmm(p_hi, p_lo, F->J_pr, v1);
            v1 = F->J * v1 - p_hi * F->mod;  // [0..2n)
            p0 = v0 + v1;  // [0..8n)
            p1 = v0 + F->mod2 - v1;  // [0..8n)
            // scale p0 <-- w**k * p0, p1 <-- w**(5*k) * p1
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][2*k+1], p0);
            pp4[k] = F->tab_w[depth-4][2*k] * p0 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][10*k+1], p1);
            pp5[k] = F->tab_w[depth-4][10*k] * p1 - p_hi * F->mod;  // [0..2n)

            // right-right, mod x - I*J | x + I*J
            umul_ppmm(p_hi, p_lo, F->IJ_pr, v3);
            v3 = F->IJ * v3 - p_hi * F->mod;  // [0..2n)
            p2 = v2 + v3;   // [0..8n)
            p3 = v2 + F->mod2 - v3;   // [0..8n)
            // scale p2 <-- w**(3*k) * p2, p3 <-- w**(7*k) * p3
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][6*k+1], p2);
            pp6[k] = F->tab_w[depth-4][6*k] * p2 - p_hi * F->mod;  // [0..2n)
            umul_ppmm(p_hi, p_lo, F->tab_w[depth-4][14*k+1], p3);
            pp7[k] = F->tab_w[depth-4][14*k] * p3 - p_hi * F->mod;  // [0..2n)
        }
        _n_fft_old_dif_rec8_lazy(pp0, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp1, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp2, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp3, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp4, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp5, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp6, len/8, depth-3, F);
        _n_fft_old_dif_rec8_lazy(pp7, len/8, depth-3, F);
    }
}





void n_fft_old_ctx_init_set(n_fft_old_ctx_t F, ulong w, ulong depth, ulong mod)
{
    // basic attributes
    F->mod = mod;
    F->mod2 = 2*mod;
    F->mod4 = 4*mod;
    F->depth = depth;
    F->w = w;

    if (depth == 3)
    {
        ulong w_pr_quo, w_pr_rem;
        n_mulmod_precomp_shoup_quo_rem(&w_pr_quo, &w_pr_rem, w, mod);

        // J == w
        F->J = w;
        F->J_pr = w_pr_quo;
        // J**2 == I
        n_mulmod_and_precomp_shoup(&F->I, &F->I_pr, w, w, w_pr_quo, w_pr_rem, w_pr_quo, mod);
        // J**3 == I * J
        n_mulmod_and_precomp_shoup(&F->IJ, &F->IJ_pr, w, F->I, w_pr_quo, w_pr_rem, F->I_pr, mod);
    }
    else
    {
        slong ell = depth-4;  // >= 0
        ulong len = (UWORD(1) << (depth-1));  // len == 2**(ell+3) >= 8, multiple of 4

        // fill array of powers of w
        F->tab_w[ell] = _nmod_vec_init(4*len);
        _n_geometric_sequence_and_opposites_with_precomp(F->tab_w[ell], w, len, mod);

        // copy into other arrays
        // NAIVE VERSION:
        //   -> if depth up to ~10, this is negligible (<~10%) compared to the above
        //   -> beyond, its impact grows and reaches a factor around 2 up to depth 30 at least
        //       (meaning the next for-ell loop takes about as much time as the above for-k loop)
        ell--;
        for (; ell >= 0; ell--)
        {
            len = len >> 1;  // len == 2**(ell+1)
            F->tab_w[ell] = _nmod_vec_init(4*len);
            F->tab_w[ell][0] = F->tab_w[ell+1][0];
            F->tab_w[ell][1] = F->tab_w[ell+1][1];
            for (ulong k = 1; k < 2*len; k++)
            {
                F->tab_w[ell][2*k] = F->tab_w[ell+1][4*k];
                F->tab_w[ell][2*k+1] = F->tab_w[ell+1][4*k+1];
            }
        }

        F->J     = F->tab_w[0][4];
        F->J_pr  = F->tab_w[0][5];
        F->I     = F->tab_w[0][8];
        F->I_pr  = F->tab_w[0][9];
        F->IJ    = F->tab_w[0][12];
        F->IJ_pr = F->tab_w[0][13];
    }
}

void n_fft_old_ctx_clear(n_fft_old_ctx_t F)
{
    for (ulong ell = 0; ell+4 <= F->depth; ell++)
    {
        _nmod_vec_clear(F->tab_w[ell]);
        //_nmod_vec_clear(F->tab_winv[ell]);
    }
}




/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
