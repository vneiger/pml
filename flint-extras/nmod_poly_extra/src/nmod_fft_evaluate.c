#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <immintrin.h>

#include "nmod_poly_fft.h"

#if defined(__GNUC__)
# define FLINT_NO_VECTORIZE __attribute__((optimize("no-tree-vectorize")))
#else
# define FLINT_NO_VECTORIZE
#endif


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

// lazy2:
// in [0..2n), out [0..2n) x [0..4n)
// n2 is 2*n, tmp is a temporary (ulong)
#define DFT2_LAZY2(a, b, n2, tmp) \
    do {                          \
        tmp = (b);                \
        (b) = (a) + (n2) - tmp;   \
        (a) = (a) + tmp;          \
    } while(0)

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
        if (u >= (n2))                                            \
            u -= (n2);  /* [0..2n) */                             \
        v = (b);                                                  \
        N_MULMOD_PRECOMP_LAZY(v, w, v, w_pr, n, p_hi, p_lo);      \
        (a) = u + v;                   /* [0..4n) */              \
        (b) = u + (n2) - v;         /* [0..4n) */                 \
    } while(0)

/*----------------*/
/* 4-point macros */
/*----------------*/


/** 4-point butterflies
 *
 * Gentleman-Sande: in-place transform
 *                              [1  w2    w1    w3]
 *                              [1 -w2  I*w1 -I*w3]
 * [a  b  c  d] <- [a  b  c  d] [1  w2   -w1   -w3]
 *                              [1 -w2 -I*w1  I*w3]
 *
 *                              [1     1   ] [1  w2        ]
 *                              [   1     I] [1 -w2        ]
 *              == [a  b  c  d] [1    -1   ] [       w1  w3]
 *                              [   1    -I] [       w1 -w3]
 **/

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










/** 4-point DFT, node 0
 * * in [0..2n) / out [0..2n) / max 4n
 * * In-place transform
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 * * Corresponds to reducing down the tree with nodes
 *                       x^4 - 1
 *                     /         \
 *             x^2 - 1             x^2 + 1
 *             /     \             /     \
 *         x - 1     x + 1     x - I     x + I
 *  where I is typically a square root of -1
 *  (but this property is not exploited)
 */
#define DFT4_LAZY2_RED(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo)         \
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


/** 4-point DFT, general
 * * in [0..4n) / out [0..4n) / max < 8n
 * * In-place transform
 *                              [ 1          1       1       1]
 *                              [w2        -w2      w3     -w3]
 * [a  b  c  d] <- [a  b  c  d] [w1         w1     -w1     -w1]
 *                              [w1*w2  -w1*w2  -w1*w3   w1*w3]
 * * Corresponds to reducing down the tree with nodes
 *                        x^4 - w1**2
 *                      /             \
 *             x^2 - w1                 x^2 + w1
 *             /      \                 /      \
 *        x - w2      x + w2       x - w3      x + w3
 * typically w2**2 == w1 and w3 == I*w2 (so that w3**2 == -w1) so that this
 * really is the subproduct tree built from the four roots
 *           w2, -w2, I*w2, -I*w2   of x**4 - w1
 */
#define DFT4_LAZY(a, b, c, d,                                 \
                  w1, w1_pr, w2, w2_pr, w3, w3_pr,            \
                  n, n2, p_hi, p_lo, tmp)                     \
do {                                                          \
    ulong u0 = (a);                                           \
    ulong u1 = (b);                                           \
    ulong u2 = (c);                                           \
    ulong u3 = (d);                                           \
    if (u0 >= n2)                                             \
        u0 -= n2;                                             \
    if (u1 >= n2)                                             \
        u1 -= n2;                                             \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u2, w1, u2, w1_pr, n, p_hi, p_lo);  \
    tmp = u0;                                                 \
    u0 = u0 + u2;                    /* [0..4n) */            \
    u2 = tmp + n2 - u2;              /* [0..4n) */            \
    if (u0 >= n2)                                             \
        u0 -= n2;                    /* [0..2n) */            \
    if (u2 >= n2)                                             \
        u2 -= n2;                    /* [0..2n) */            \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u3, w1, u3, w1_pr, n, p_hi, p_lo);  \
    tmp = u1;                                                 \
    u1 = u1 + u3;                    /* [0..8n) */            \
    u3 = tmp + n2 - u3;              /* [0..8n) */            \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u1, w2, u1, w2_pr, n, p_hi, p_lo);  \
    tmp = u0;                                                 \
    (a) = u0 + u1;                    /* [0..4n) */           \
    (b) = tmp + n2 - u1;              /* [0..4n) */           \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u3, w3, u3, w3_pr, n, p_hi, p_lo);  \
    tmp = u2;                                                 \
    (c) = u2 + u3;                    /* [0..4n) */           \
    (d) = tmp + n2 - u3;              /* [0..4n) */           \
} while(0)

/*----------------*/
/* 8-point DFT    */
/*----------------*/

// reduction tree 8-point lazy DFT
// lazy red: input in [0..2*n) --> output in [0..4*n)
// max value < 8n
FLINT_FORCE_INLINE void dft8_red_lazy(ulong * p0,
                                      ulong * p1,
                                      ulong * p2,
                                      ulong * p3,
                                      ulong * p4,
                                      ulong * p5,
                                      ulong * p6,
                                      ulong * p7,
                                      n_fft_ctx_t F)
//{
//    ulong p_hi, p_lo;
//
//    ulong u0 = *p0;
//    ulong u1 = *p1;
//    ulong u2 = *p2;
//    ulong u3 = *p3;
//    ulong v0 = *p4;
//    ulong v1 = *p5;
//    ulong v2 = *p6;
//    ulong v3 = *p7;
//
//    // mod x**4 - 1 | x**4 + 1
//    ulong t0 = u0 + v0;  // [0..4n)
//    ulong t1 = u1 + v1;  // [0..4n)
//    ulong t2 = u2 + v2;  // [0..4n)
//    ulong t3 = u3 + v3;  // [0..4n)
//    u0 += F->mod2 - v0;  // [0..4n)
//    u1 += F->mod2 - v1;  // [0..4n)
//    u2 += F->mod2 - v2;  // [0..4n)
//    u3 += F->mod2 - v3;  // [0..4n)
//
//    // left, mod x**2 - 1 | x**2 + 1
//    v0 = t0 + t2;             // [0..8n)
//    v1 = t1 + t3;             // [0..8n)
//    v2 = t0 + F->mod4 - t2;  // [0..8n)
//    v3 = t1 + F->mod4 - t3;  // [0..8n)
//
//    // left-left, mod x-1 | x+1
//    if (v0 >= F->mod4)
//        v0 -= F->mod4;
//    if (v1 >= F->mod4)
//        v1 -= F->mod4;
//    t0 = v0 + v1;               // [0..8n)
//    t1 = v0 + F->mod4 - v1;     // [0..8n)
//    if (t0 >= F->mod4)
//        t0 -= F->mod4;
//    if (t1 >= F->mod4)
//        t1 -= F->mod4;
//    *p0 = t0;
//    *p1 = t1;
//
//    // left-right, mod x-I | x+I
//    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[2], v3, F->tab_w[3], F->mod, p_hi, p_lo);
//    if (v2 >= F->mod4)
//        v2 -= F->mod4;
//    if (v2 >= F->mod2)
//        v2 -= F->mod2;         // [0..2n)
//    *p2 = v2 + v3;             // [0..4n)
//    *p3 = v2 + F->mod2 - v3;  // [0..4n)
//
//    // right, mod x**2 - I | x**2 + I
//    N_MULMOD_PRECOMP_LAZY(u2, F->tab_w[2], u2, F->tab_w[3], F->mod, p_hi, p_lo);
//    N_MULMOD_PRECOMP_LAZY(u3, F->tab_w[2], u3, F->tab_w[3], F->mod, p_hi, p_lo);
//    if (u0 >= F->mod2)
//        u0 -= F->mod2;         // [0..2n)
//    if (u1 >= F->mod2)
//        u1 -= F->mod2;         // [0..2n)
//    v0 = u0 + u2;  // [0..4n)
//    v1 = u1 + u3;  // [0..4n)
//    v2 = u0 + F->mod2 - u2;  // [0..4n)
//    v3 = u1 + F->mod2 - u3;  // [0..4n)
//
//    // right-left, mod x - J | x + J
//    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[4], v1, F->tab_w[5], F->mod, p_hi, p_lo);
//    if (v0 >= F->mod2)
//        v0 -= F->mod2;         // [0..2n)
//    *p4 = v0 + v1;
//    *p5 = v0 + F->mod2 - v1;
//
//    // right-right, mod x - I*J | x + I*J
//    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[6], v3, F->tab_w[7], F->mod, p_hi, p_lo);
//    if (v2 >= F->mod2)
//        v2 -= F->mod2;         // [0..2n)
//    *p6 = v2 + v3;
//    *p7 = v2 + F->mod2 - v3;
//}

//{
//     ulong p_hi, p_lo, u, v;
//
//     DFT4_LAZY2_RED(*p0, *p2, *p4, *p6,
//                    F->tab_w[2], F->tab_w[3],
//                    F->mod, F->mod2, p_hi, p_lo);
//     DFT4_LAZY2_RED(*p1, *p3, *p5, *p7,
//                    F->tab_w[2], F->tab_w[3],
//                    F->mod, F->mod2, p_hi, p_lo);
//
//     BUTTERFLY2_CT_LAZY(*p0, *p1, F->mod, F->mod2, F->tab_w[0], F->tab_w[1], p_hi, p_lo, u, v);
//     BUTTERFLY2_CT_LAZY(*p2, *p3, F->mod, F->mod2, F->tab_w[2], F->tab_w[3], p_hi, p_lo, u, v);
//     BUTTERFLY2_CT_LAZY(*p4, *p5, F->mod, F->mod2, F->tab_w[4], F->tab_w[5], p_hi, p_lo, u, v);
//     BUTTERFLY2_CT_LAZY(*p6, *p7, F->mod, F->mod2, F->tab_w[6], F->tab_w[7], p_hi, p_lo, u, v);
//}

   //     // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
   //     // (general accepts [0..4n) as input for depth >= 3)
   //     ulong tmp;
   //     for (ulong k = 0; k < len/2; k++)
   //         DFT2_LAZY2_RED1(p[k+0], p[len/2+k+0], F->mod2, tmp);
   //     _n_fft_red_rec2_lazy(p, len/2, depth-1, F);
   //     _n_fft_red_rec2_lazy_general(p+len/2, len/2, depth-1, 1, F);
{
    ulong p_hi, p_lo, tmp;

    DFT2_LAZY2_RED1(*p0, *p4, F->mod2, tmp);
    DFT2_LAZY2_RED1(*p1, *p5, F->mod2, tmp);
    DFT2_LAZY2_RED1(*p2, *p6, F->mod2, tmp);
    DFT2_LAZY2_RED1(*p3, *p7, F->mod2, tmp);
    DFT4_LAZY2_RED(*p0, *p1, *p2, *p3, F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    DFT4_LAZY(*p4, *p5, *p6, *p7,
              F->tab_w[2], F->tab_w[3],
              F->tab_w[4], F->tab_w[5],
              F->tab_w[6], F->tab_w[7],
              F->mod, F->mod2, p_hi, p_lo, tmp);
}

// in [0..4n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft8_red_lazy_general(ulong * p0,
                                              ulong * p1,
                                              ulong * p2,
                                              ulong * p3,
                                              ulong * p4,
                                              ulong * p5,
                                              ulong * p6,
                                              ulong * p7,
                                              ulong node,
                                              n_fft_ctx_t F)
/* FIRST VERSION of dft8_red_lazy_general */
//{
//    ulong p_hi, p_lo;
//
//    ulong u0 = *p0;
//    ulong u1 = *p1;
//    ulong u2 = *p2;
//    ulong u3 = *p3;
//    ulong v0 = *p4;
//    ulong v1 = *p5;
//    ulong v2 = *p6;
//    ulong v3 = *p7;
//
//    // mod x**4 - w | x**4 + w
//    ulong w = F->tab_w[2*node];
//    ulong wpre = F->tab_w[2*node+1];
//    N_MULMOD_PRECOMP_LAZY(v0, w, v0, wpre, F->mod, p_hi, p_lo);
//    N_MULMOD_PRECOMP_LAZY(v1, w, v1, wpre, F->mod, p_hi, p_lo);
//    N_MULMOD_PRECOMP_LAZY(v2, w, v2, wpre, F->mod, p_hi, p_lo);
//    N_MULMOD_PRECOMP_LAZY(v3, w, v3, wpre, F->mod, p_hi, p_lo);
//    ulong pp0 = u0 + v0;   // [0..6n)
//    ulong pp1 = u1 + v1;   // [0..6n)
//    ulong pp2 = u2 + v2;   // [0..6n)
//    ulong pp3 = u3 + v3;   // [0..6n)
//    u0 += F->mod2 - v0;  // [0..6n)
//    u1 += F->mod2 - v1;  // [0..6n)
//    u2 += F->mod2 - v2;  // [0..6n)
//    u3 += F->mod2 - v3;  // [0..6n)
//
//    // left, mod x**2 - sqrt(w) | x**2 + sqrt(w)
//    w = F->tab_w[4*node];
//    wpre = F->tab_w[4*node+1];
//    N_MULMOD_PRECOMP_LAZY(pp2, w, pp2, wpre, F->mod, p_hi, p_lo);
//    N_MULMOD_PRECOMP_LAZY(pp3, w, pp3, wpre, F->mod, p_hi, p_lo);
//    v0 = pp0 + pp2;             // [0..8n)
//    v1 = pp1 + pp3;             // [0..8n)
//    v2 = pp0 + F->mod2 - pp2;  // [0..8n)
//    v3 = pp1 + F->mod2 - pp3;  // [0..8n)
//
//    // left-left, mod x - fort(w) | x + fort(w)
//    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[8*node], v1, F->tab_w[8*node+1], F->mod, p_hi, p_lo);
//    if (v0 >= F->mod4)
//        v0 -= F->mod4;
//    if (v0 >= F->mod2)
//        v0 -= F->mod2;  // [0..2n)
//    *p0 = v0 + v1;             // [0..4n)
//    *p1 = v0 + F->mod2 - v1;  // [0..4n)
//
//    // left-right, mod x - I*fort(w) | x+ I*fort(w)
//    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[8*node+2], v3, F->tab_w[8*node+3], F->mod, p_hi, p_lo);
//    if (v2 >= F->mod4)
//        v2 -= F->mod4;
//    if (v2 >= F->mod2)
//        v2 -= F->mod2;  // [0..2n)
//    *p2 = v2 + v3;              // [0..4n)
//    *p3 = v2 + F->mod2 - v3;   // [0..4n)
//
//    // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
//    w = F->tab_w[4*node+2];
//    wpre = F->tab_w[4*node+3];
//    N_MULMOD_PRECOMP_LAZY(u2, w, u2, wpre, F->mod, p_hi, p_lo);
//    N_MULMOD_PRECOMP_LAZY(u3, w, u3, wpre, F->mod, p_hi, p_lo);
//    v0 = u0 + u2;             // [0..8n)
//    v1 = u1 + u3;             // [0..8n)
//    v2 = u0 + F->mod2 - u2;  // [0..8n)
//    v3 = u1 + F->mod2 - u3;  // [0..8n)
//
//    // right-left, mod x - J*fort(w) | x + J*fort(w)
//    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[8*node+4], v1, F->tab_w[8*node+5], F->mod, p_hi, p_lo);
//    if (v0 >= F->mod4)
//        v0 -= F->mod4;
//    if (v0 >= F->mod2)
//        v0 -= F->mod2;  // [0..2n)
//    *p4 = v0 + v1;             // [0..4n)
//    *p5 = v0 + F->mod2 - v1;  // [0..4n)
//
//    // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
//    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[8*node+6], v3, F->tab_w[8*node+7], F->mod, p_hi, p_lo);
//    if (v2 >= F->mod4)
//        v2 -= F->mod4;
//    if (v2 >= F->mod2)
//        v2 -= F->mod2;  // [0..2n)
//    *p6 = v2 + v3;             // [0..4n)
//    *p7 = v2 + F->mod2 - v3;  // [0..4n)
//}

/* ALTERNATIVE to dft8 general: */
//{
//     ulong p_hi, p_lo, u, v;
//
//     DFT4_LAZY(*p0, *p2, *p4, *p6,
//                        F->tab_w[2*node], F->tab_w[2*node+1],
//                        F->tab_w[4*node], F->tab_w[4*node+1],
//                        F->tab_w[4*node+2], F->tab_w[4*node+3],
//                        F->mod, F->mod2, p_hi, p_lo, u);
//     DFT4_LAZY(*p1, *p3, *p5, *p7,
//                        F->tab_w[2*node], F->tab_w[2*node+1],
//                        F->tab_w[4*node], F->tab_w[4*node+1],
//                        F->tab_w[4*node+2], F->tab_w[4*node+3],
//                        F->mod, F->mod2, p_hi, p_lo, u);
//
//     BUTTERFLY2_CT_LAZY(*p0, *p1, F->mod, F->mod2, F->tab_w[8*node], F->tab_w[8*node+1], p_hi, p_lo, u, v);
//     BUTTERFLY2_CT_LAZY(*p2, *p3, F->mod, F->mod2, F->tab_w[8*node+2], F->tab_w[8*node+3], p_hi, p_lo, u, v);
//     BUTTERFLY2_CT_LAZY(*p4, *p5, F->mod, F->mod2, F->tab_w[8*node+4], F->tab_w[8*node+5], p_hi, p_lo, u, v);
//     BUTTERFLY2_CT_LAZY(*p6, *p7, F->mod, F->mod2, F->tab_w[8*node+6], F->tab_w[8*node+7], p_hi, p_lo, u, v);
//}

// SECOND ALTERNATIVE to dft8 general:
{
    ulong p_hi, p_lo, u, v;

    const ulong w = F->tab_w[2*node];
    const ulong wpre = F->tab_w[2*node+1];

    BUTTERFLY2_CT_LAZY(*p0, *p4, F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
    BUTTERFLY2_CT_LAZY(*p1, *p5, F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
    BUTTERFLY2_CT_LAZY(*p2, *p6, F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
    BUTTERFLY2_CT_LAZY(*p3, *p7, F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);

    DFT4_LAZY(*p0, *p1, *p2, *p3,
              F->tab_w[4*node], F->tab_w[4*node+1],
              F->tab_w[8*node], F->tab_w[8*node+1],
              F->tab_w[8*node+2], F->tab_w[8*node+3],
              F->mod, F->mod2, p_hi, p_lo, u);
    DFT4_LAZY(*p4, *p5, *p6, *p7,
              F->tab_w[4*node+2], F->tab_w[4*node+3],
              F->tab_w[8*node+4], F->tab_w[8*node+5],
              F->tab_w[8*node+6], F->tab_w[8*node+7],
              F->mod, F->mod2, p_hi, p_lo, u);
}

/*------------------*/
/* other base cases */
/*------------------*/

// in [0..2n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft16_red_lazy(nn_ptr p, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT4_LAZY2_RED(p[0], p[4], p[ 8], p[12], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_LAZY2_RED(p[1], p[5], p[ 9], p[13], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_LAZY2_RED(p[2], p[6], p[10], p[14], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_LAZY2_RED(p[3], p[7], p[11], p[15], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;

    // next line requires < 2n, hence the four reductions above
    DFT4_LAZY2_RED(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    DFT4_LAZY(p[4], p[5], p[6], p[7], F->tab_w[2], F->tab_w[3], F->tab_w[4], F->tab_w[5], F->tab_w[6], F->tab_w[7], F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[8], p[9], p[10], p[11], F->tab_w[4], F->tab_w[5], F->tab_w[8], F->tab_w[9], F->tab_w[10], F->tab_w[11], F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[12], p[13], p[14], p[15], F->tab_w[6], F->tab_w[7], F->tab_w[12], F->tab_w[13], F->tab_w[14], F->tab_w[15], F->mod, F->mod2, p_hi, p_lo, tmp);
}

FLINT_FORCE_INLINE void dft16_red_lazy_general(nn_ptr p, ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    ulong w2 = F->tab_w[2*node];
    ulong w2pre = F->tab_w[2*node+1];
    ulong w = F->tab_w[4*node];
    ulong wpre = F->tab_w[4*node+1];
    ulong Iw = F->tab_w[4*node+2];
    ulong Iwpre = F->tab_w[4*node+3];
    DFT4_LAZY(p[0], p[4], p[8], p[12], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[1], p[5], p[9], p[13], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[2], p[6], p[10], p[14], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[3], p[7], p[11], p[15], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node];
    w2pre = F->tab_w[8*node+1];
    w = F->tab_w[16*node];
    wpre = F->tab_w[16*node+1];
    Iw = F->tab_w[16*node+2];
    Iwpre = F->tab_w[16*node+3];
    DFT4_LAZY(p[0], p[1], p[2], p[3], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+2];
    w2pre = F->tab_w[8*node+3];
    w = F->tab_w[16*node+4];
    wpre = F->tab_w[16*node+5];
    Iw = F->tab_w[16*node+6];
    Iwpre = F->tab_w[16*node+7];
    DFT4_LAZY(p[4], p[5], p[6], p[7], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+4];
    w2pre = F->tab_w[8*node+5];
    w = F->tab_w[16*node+8];
    wpre = F->tab_w[16*node+9];
    Iw = F->tab_w[16*node+10];
    Iwpre = F->tab_w[16*node+11];
    DFT4_LAZY(p[8], p[9], p[10], p[11], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+6];
    w2pre = F->tab_w[8*node+7];
    w = F->tab_w[16*node+12];
    wpre = F->tab_w[16*node+13];
    Iw = F->tab_w[16*node+14];
    Iwpre = F->tab_w[16*node+15];
    DFT4_LAZY(p[12], p[13], p[14], p[15], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
}

/* ALTERNATIVE to dft16 general: */
//{
//    const ulong w = F->tab_w[2*node];
//    const ulong wpre = F->tab_w[2*node+1];
//    ulong p_hi, p_lo, u, v;
//    BUTTERFLY2_CT_LAZY(p[0], p[8], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[1], p[9], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[2], p[10], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[3], p[11], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[4], p[12], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[5], p[13], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[6], p[14], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    BUTTERFLY2_CT_LAZY(p[7], p[15], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
//    dft8_red_lazy_general(p, 2*node, F);
//    dft8_red_lazy_general(p+8, 2*node+1, F);
//}

// in [0..2n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft32_red_lazy(nn_ptr p, n_fft_ctx_t F)
{
    ulong p_hi, p_lo;

    DFT4_LAZY2_RED(p[0], p[8 ], p[16], p[24], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_LAZY2_RED(p[1], p[9 ], p[17], p[25], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_LAZY2_RED(p[2], p[10], p[18], p[26], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_LAZY2_RED(p[3], p[11], p[19], p[27], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;
    DFT4_LAZY2_RED(p[4], p[12], p[20], p[28], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[4] >= F->mod2)
        p[4] -= F->mod2;
    DFT4_LAZY2_RED(p[5], p[13], p[21], p[29], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[5] >= F->mod2)
        p[5] -= F->mod2;
    DFT4_LAZY2_RED(p[6], p[14], p[22], p[30], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[6] >= F->mod2)
        p[6] -= F->mod2;
    DFT4_LAZY2_RED(p[7], p[15], p[23], p[31], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[7] >= F->mod2)
        p[7] -= F->mod2;

    // next line requires < 2n, hence the 8 reductions above
    dft8_red_lazy(        p+ 0, p+ 1, p+ 2, p+ 3, p+ 4, p+ 5, p+ 6, p+ 7,    F);
    dft8_red_lazy_general(p+ 8, p+ 9, p+10, p+11, p+12, p+13, p+14, p+15, 1, F);
    dft8_red_lazy_general(p+16, p+17, p+18, p+19, p+20, p+21, p+22, p+23, 2, F);
    dft8_red_lazy_general(p+24, p+25, p+26, p+27, p+28, p+29, p+30, p+31, 3, F);
}


// in [0..2n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft32_red_lazy_general(nn_ptr p, ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    ulong w2 = F->tab_w[2*node];
    ulong w2pre = F->tab_w[2*node+1];
    ulong w = F->tab_w[4*node];
    ulong wpre = F->tab_w[4*node+1];
    ulong Iw = F->tab_w[4*node+2];
    ulong Iwpre = F->tab_w[4*node+3];
    DFT4_LAZY(p[0], p[ 8], p[16], p[24], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[1], p[ 9], p[17], p[25], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[2], p[10], p[18], p[26], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[3], p[11], p[19], p[27], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[4], p[12], p[20], p[28], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[5], p[13], p[21], p[29], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[6], p[14], p[22], p[30], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[7], p[15], p[23], p[31], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    // next line requires < 2n, hence the four reductions above
    dft8_red_lazy_general(p+ 0, p+ 1, p+ 2, p+ 3, p+ 4, p+ 5, p+ 6, p+ 7, 4*node, F);
    dft8_red_lazy_general(p+ 8, p+ 9, p+10, p+11, p+12, p+13, p+14, p+15, 4*node+1, F);
    dft8_red_lazy_general(p+16, p+17, p+18, p+19, p+20, p+21, p+22, p+23, 4*node+2, F);
    dft8_red_lazy_general(p+24, p+25, p+26, p+27, p+28, p+29, p+30, p+31, 4*node+3, F);
}






// if depth < 3, in [0..2n) out [0..4n)
// if depth >= 3, in [0..4n) out [0..4n)
void _n_fft_red_rec2_lazy_general(nn_ptr p, ulong len, ulong depth, ulong node, n_fft_ctx_t F)
{
    if (depth == 3)
        dft8_red_lazy_general(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, node, F);
    else if (depth == 4)
        dft16_red_lazy_general(p, node, F);
    else if (depth == 5)
        dft32_red_lazy_general(p, node, F);
    else
    {
        // in: [0..4n), out: [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+len/2;
        const ulong w = F->tab_w[2*node];
        const ulong wpre = F->tab_w[2*node+1];
        ulong p_hi, p_lo, u, v;
        for (ulong k = 0; k < len/2; k+=4)
        {
            BUTTERFLY2_CT_LAZY(p0[k+0], p1[k+0], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            BUTTERFLY2_CT_LAZY(p0[k+1], p1[k+1], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            BUTTERFLY2_CT_LAZY(p0[k+2], p1[k+2], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            BUTTERFLY2_CT_LAZY(p0[k+3], p1[k+3], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
        }
        _n_fft_red_rec2_lazy_general(p0, len/2, depth-1, 2*node, F);
        _n_fft_red_rec2_lazy_general(p1, len/2, depth-1, 2*node+1, F);
    }
}

// input [0..2n),  output [0..4n)
void _n_fft_red_rec2_lazy(nn_ptr p, ulong len, ulong depth, n_fft_ctx_t F)
{
    // depth == 0: nothing to do
    if (depth == 1) // in [0..2n), out [0..4n)
    {
        ulong tmp;
        DFT2_LAZY2(p[0], p[1], F->mod2, tmp);
    }
    else if (depth == 2) // in [0..2n), out [0..4n)
    {
        ulong p_hi, p_lo;
        DFT4_LAZY2_RED(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 3)
        dft8_red_lazy(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, F);  // in [0..2n), out [0..4n)
    else if (depth == 4)
        dft16_red_lazy(p, F);  // in [0..2n), out [0..4n)
    else if (depth == 5)
        dft32_red_lazy(p, F);  // in [0..2n), out [0..4n)
    else
    {
        // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
        // (general accepts [0..4n) as input for depth >= 3)
        ulong tmp;
        for (ulong k = 0; k < len/2; k++)
            DFT2_LAZY2_RED1(p[k+0], p[len/2+k+0], F->mod2, tmp);
        _n_fft_red_rec2_lazy(p, len/2, depth-1, F);
        _n_fft_red_rec2_lazy_general(p+len/2, len/2, depth-1, 1, F);
    }
}



void _n_fft_red_iter2_lazy(nn_ptr p, ulong len, ulong depth, n_fft_ctx_t F)
{
    // perform FFT layers up to depth 3
    ulong llen = len;
    for (ulong ell = depth; ell > 3; ell--, llen>>=1)
    {
        // k = 0: point is 1, handle separately
        ulong tmp;
        // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
        for (ulong kk = 0; kk+3 < llen/2; kk+=4)
        {
            DFT2_LAZY2_RED1(p[kk+0], p[llen/2+kk+0], F->mod2, tmp);
            DFT2_LAZY2_RED1(p[kk+1], p[llen/2+kk+1], F->mod2, tmp);
            DFT2_LAZY2_RED1(p[kk+2], p[llen/2+kk+2], F->mod2, tmp);
            DFT2_LAZY2_RED1(p[kk+3], p[llen/2+kk+3], F->mod2, tmp);
        }

        ulong node = 1;  // index of current point in tab
        for (ulong k = llen; k < len; k+=llen, node++)
        {
            ulong w = F->tab_w[2*node];
            ulong wpre = F->tab_w[2*node+1];

            const nn_ptr p0 = p+k;
            const nn_ptr p1 = p+llen/2+k;
            ulong p_hi, p_lo, u, v;
            for (ulong kk = 0; kk < llen/2; kk+=4)
            {
                BUTTERFLY2_CT_LAZY(p0[kk+0], p1[kk+0], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
                BUTTERFLY2_CT_LAZY(p0[kk+1], p1[kk+1], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
                BUTTERFLY2_CT_LAZY(p0[kk+2], p1[kk+2], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
                BUTTERFLY2_CT_LAZY(p0[kk+3], p1[kk+3], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            }
        }
    }

    // perform last two FFT layers
    dft8_red_lazy(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, F);  // k == 0..7, in [0..2n) (is ok!), out [0..4n)
    ulong node = 1;
    for (ulong k = 8; k < len; k+=8, node++)
        dft8_red_lazy_general(p+k+0, p+k+1, p+k+2, p+k+3, p+k+4, p+k+5, p+k+6, p+k+7, node, F);  // in [0..4n), out [0..4n)
}

void _n_fft_red_rec4_lazy_general(nn_ptr p, ulong len, ulong depth, ulong node, n_fft_ctx_t F)
{
    if (depth == 3)
        dft8_red_lazy_general(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, node, F);
    else if (depth == 4)
        dft16_red_lazy_general(p, node, F);
    else if (depth == 5)  // TODO unclear this helps (no acceleration on argiope)
        dft32_red_lazy_general(p, node, F);
    else
    {
        // in: [0..4n), out: [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+len/4;
        const nn_ptr p2 = p+2*len/4;
        const nn_ptr p3 = p+3*len/4;
        const ulong w2 = F->tab_w[2*node];
        const ulong w2pre = F->tab_w[2*node+1];
        const ulong w = F->tab_w[4*node];
        const ulong wpre = F->tab_w[4*node+1];
        const ulong Iw = F->tab_w[4*node+2];
        const ulong Iwpre = F->tab_w[4*node+3];
        ulong p_hi, p_lo, tmp;

        for (ulong k = 0; k < len/4; k+=4)
        {
            DFT4_LAZY(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_LAZY(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_LAZY(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_LAZY(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
        }

        _n_fft_red_rec4_lazy_general(p0, len/4, depth-2, 4*node, F);
        _n_fft_red_rec4_lazy_general(p1, len/4, depth-2, 4*node+1, F);
        _n_fft_red_rec4_lazy_general(p2, len/4, depth-2, 4*node+2, F);
        _n_fft_red_rec4_lazy_general(p3, len/4, depth-2, 4*node+3, F);
    }
}

// input [0..2*n), output [0..4*n)
// depth >= 3
void _n_fft_red_rec4_lazy(nn_ptr p, ulong len, ulong depth, n_fft_ctx_t F)
{
    // depth == 0: nothing to do
    //if (depth == 1)
    //    // in [0..4n), out [0..4n)
    //    DFT2_LAZY4_RED(p[0], p[1], F->mod4);
    //else if (depth == 2)
    //    // in [0..2n), out [0..4n)
    //    DFT4_LAZY2_RED(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
    //else
    if (depth == 3)
        dft8_red_lazy(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, F);  // in [0..2n), out [0..4n)
    else if (depth == 4)
        dft16_red_lazy(p, F);  // in [0..2n), out [0..4n)
    else if (depth == 5)   // TODO unclear this helps (no acceleration on argiope)
        dft32_red_lazy(p, F);  // in [0..2n), out [0..4n)
    else
    {
        // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
        // (general accepts [0..4n) as input for depth >= 3)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/4;
        const nn_ptr p2 = p + 2*len/4;
        const nn_ptr p3 = p + 3*len/4;
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k++)
        {
            DFT4_LAZY2_RED(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            if (p[k] >= F->mod2)
                p[k] -= F->mod2;
        }
        _n_fft_red_rec4_lazy(p0, len/4, depth-2, F);
        _n_fft_red_rec4_lazy_general(p1, len/4, depth-2, 1, F);
        _n_fft_red_rec4_lazy_general(p2, len/4, depth-2, 2, F);
        _n_fft_red_rec4_lazy_general(p3, len/4, depth-2, 3, F);
    }
}









/*-----------------------*/
/*-----------------------*/
/* VERSIONS WITH STRIDES */
/*-----------------------*/
/*-----------------------*/

// in [0..4n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft8_red_lazy_general_stride(nn_ptr p, ulong stride, ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo;

    ulong u0 = p[0*stride];
    ulong u1 = p[1*stride];
    ulong u2 = p[2*stride];
    ulong u3 = p[3*stride];
    ulong v0 = p[4*stride];
    ulong v1 = p[5*stride];
    ulong v2 = p[6*stride];
    ulong v3 = p[7*stride];

    // mod x**4 - w | x**4 + w
    ulong w = F->tab_w[2*node];
    ulong wpre = F->tab_w[2*node+1];
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
    w = F->tab_w[4*node];
    wpre = F->tab_w[4*node+1];
    N_MULMOD_PRECOMP_LAZY(p2, w, p2, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(p3, w, p3, wpre, F->mod, p_hi, p_lo);
    v0 = p0 + p2;             // [0..8n)
    v1 = p1 + p3;             // [0..8n)
    v2 = p0 + F->mod2 - p2;  // [0..8n)
    v3 = p1 + F->mod2 - p3;  // [0..8n)

    // left-left, mod x - fort(w) | x + fort(w)
    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[8*node], v1, F->tab_w[8*node+1], F->mod, p_hi, p_lo);
    if (v0 >= F->mod4)
        v0 -= F->mod4;
    if (v0 >= F->mod2)
        v0 -= F->mod2;  // [0..2n)
    p[0*stride] = v0 + v1;             // [0..4n)
    p[1*stride] = v0 + F->mod2 - v1;  // [0..4n)

    // left-right, mod x - I*fort(w) | x+ I*fort(w)
    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[8*node+2], v3, F->tab_w[8*node+3], F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;  // [0..2n)
    p[2*stride] = v2 + v3;              // [0..4n)
    p[3*stride] = v2 + F->mod2 - v3;   // [0..4n)

    // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
    w = F->tab_w[4*node+2];
    wpre = F->tab_w[4*node+3];
    N_MULMOD_PRECOMP_LAZY(u2, w, u2, wpre, F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(u3, w, u3, wpre, F->mod, p_hi, p_lo);
    v0 = u0 + u2;             // [0..8n)
    v1 = u1 + u3;             // [0..8n)
    v2 = u0 + F->mod2 - u2;  // [0..8n)
    v3 = u1 + F->mod2 - u3;  // [0..8n)

    // right-left, mod x - J*fort(w) | x + J*fort(w)
    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[8*node+4], v1, F->tab_w[8*node+5], F->mod, p_hi, p_lo);
    if (v0 >= F->mod4)
        v0 -= F->mod4;
    if (v0 >= F->mod2)
        v0 -= F->mod2;  // [0..2n)
    p[4*stride] = v0 + v1;             // [0..4n)
    p[5*stride] = v0 + F->mod2 - v1;  // [0..4n)

    // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[8*node+6], v3, F->tab_w[8*node+7], F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;  // [0..2n)
    p[6*stride] = v2 + v3;             // [0..4n)
    p[7*stride] = v2 + F->mod2 - v3;  // [0..4n)
}

FLINT_FORCE_INLINE void dft8_red_lazy_stride(nn_ptr p, ulong stride, n_fft_ctx_t F)
{
    ulong p_hi, p_lo;
    ulong u0 = p[0*stride];
    ulong u1 = p[1*stride];
    ulong u2 = p[2*stride];
    ulong u3 = p[3*stride];
    ulong v0 = p[4*stride];
    ulong v1 = p[5*stride];
    ulong v2 = p[6*stride];
    ulong v3 = p[7*stride];

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
    p[0*stride] = p0;
    p[1*stride] = p1;

    // left-right, mod x-I | x+I
    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[2], v3, F->tab_w[3], F->mod, p_hi, p_lo);
    if (v2 >= F->mod4)
        v2 -= F->mod4;
    if (v2 >= F->mod2)
        v2 -= F->mod2;         // [0..2n)
    p[2*stride] = v2 + v3;             // [0..4n)
    p[3*stride] = v2 + F->mod2 - v3;  // [0..4n)

    // right, mod x**2 - I | x**2 + I
    N_MULMOD_PRECOMP_LAZY(u2, F->tab_w[2], u2, F->tab_w[3], F->mod, p_hi, p_lo);
    N_MULMOD_PRECOMP_LAZY(u3, F->tab_w[2], u3, F->tab_w[3], F->mod, p_hi, p_lo);
    if (u0 >= F->mod2)
        u0 -= F->mod2;         // [0..2n)
    if (u1 >= F->mod2)
        u1 -= F->mod2;         // [0..2n)
    v0 = u0 + u2;  // [0..4n)
    v1 = u1 + u3;  // [0..4n)
    v2 = u0 + F->mod2 - u2;  // [0..4n)
    v3 = u1 + F->mod2 - u3;  // [0..4n)

    // right-left, mod x - J | x + J
    N_MULMOD_PRECOMP_LAZY(v1, F->tab_w[4], v1, F->tab_w[5], F->mod, p_hi, p_lo);
    if (v0 >= F->mod2)
        v0 -= F->mod2;         // [0..2n)
    p[4*stride] = v0 + v1;
    p[5*stride] = v0 + F->mod2 - v1;

    // right-right, mod x - I*J | x + I*J
    N_MULMOD_PRECOMP_LAZY(v3, F->tab_w[6], v3, F->tab_w[7], F->mod, p_hi, p_lo);
    if (v2 >= F->mod2)
        v2 -= F->mod2;         // [0..2n)
    p[6*stride] = v2 + v3;
    p[7*stride] = v2 + F->mod2 - v3;
}

// in [0..2n), out [0..4n), max value < 8n
FLINT_FORCE_INLINE void dft16_red_lazy_stride(nn_ptr p, ulong stride, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT4_LAZY2_RED(p[0*stride], p[4*stride], p[8*stride], p[12*stride], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[0*stride] >= F->mod2)
        p[0*stride] -= F->mod2;
    DFT4_LAZY2_RED(p[1*stride], p[5*stride], p[9 ], p[13*stride], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[1*stride] >= F->mod2)
        p[1*stride] -= F->mod2;
    DFT4_LAZY2_RED(p[2*stride], p[6*stride], p[10*stride], p[14*stride], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[2*stride] >= F->mod2)
        p[2*stride] -= F->mod2;
    DFT4_LAZY2_RED(p[3*stride], p[7*stride], p[11*stride], p[15*stride], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[3*stride] >= F->mod2)
        p[3*stride] -= F->mod2;

    // next line requires < 2n, hence the four reductions above
    DFT4_LAZY2_RED(p[0*stride], p[1*stride], p[2*stride], p[3*stride], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    DFT4_LAZY(p[4*stride], p[5*stride], p[6*stride], p[7*stride], F->tab_w[2], F->tab_w[3], F->tab_w[4], F->tab_w[5], F->tab_w[6], F->tab_w[7], F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[8*stride], p[9*stride], p[10*stride], p[11*stride], F->tab_w[4], F->tab_w[5], F->tab_w[8], F->tab_w[9], F->tab_w[10], F->tab_w[11], F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[12*stride], p[13*stride], p[14*stride], p[15*stride], F->tab_w[6], F->tab_w[7], F->tab_w[12], F->tab_w[13], F->tab_w[14], F->tab_w[15], F->mod, F->mod2, p_hi, p_lo, tmp);
}


FLINT_FORCE_INLINE void dft16_red_lazy_general_stride(nn_ptr p, ulong stride, ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    ulong w2 = F->tab_w[2*node];
    ulong w2pre = F->tab_w[2*node+1];
    ulong w = F->tab_w[4*node];
    ulong wpre = F->tab_w[4*node+1];
    ulong Iw = F->tab_w[4*node+2];
    ulong Iwpre = F->tab_w[4*node+3];
    DFT4_LAZY(p[0*stride], p[4*stride], p[8*stride], p[12*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[1*stride], p[5*stride], p[9*stride], p[13*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[2*stride], p[6*stride], p[10*stride], p[14*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY(p[3*stride], p[7*stride], p[11*stride], p[15*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node];
    w2pre = F->tab_w[8*node+1];
    w = F->tab_w[16*node];
    wpre = F->tab_w[16*node+1];
    Iw = F->tab_w[16*node+2];
    Iwpre = F->tab_w[16*node+3];
    DFT4_LAZY(p[0*stride], p[1*stride], p[2*stride], p[3*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+2];
    w2pre = F->tab_w[8*node+3];
    w = F->tab_w[16*node+4];
    wpre = F->tab_w[16*node+5];
    Iw = F->tab_w[16*node+6];
    Iwpre = F->tab_w[16*node+7];
    DFT4_LAZY(p[4*stride], p[5*stride], p[6*stride], p[7*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+4];
    w2pre = F->tab_w[8*node+5];
    w = F->tab_w[16*node+8];
    wpre = F->tab_w[16*node+9];
    Iw = F->tab_w[16*node+10];
    Iwpre = F->tab_w[16*node+11];
    DFT4_LAZY(p[8*stride], p[9*stride], p[10*stride], p[11*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+6];
    w2pre = F->tab_w[8*node+7];
    w = F->tab_w[16*node+12];
    wpre = F->tab_w[16*node+13];
    Iw = F->tab_w[16*node+14];
    Iwpre = F->tab_w[16*node+15];
    DFT4_LAZY(p[12*stride], p[13*stride], p[14*stride], p[15*stride], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
}

// if depth < 3, in [0..2n) out [0..4n)
// if depth >= 3, in [0..4n) out [0..4n)
void _n_fft_red_rec2_lazy_general_stride(nn_ptr p, ulong len, ulong depth, ulong node, n_fft_ctx_t F)
{
    if (depth == 3)
    {
        for (ulong k = 0; k < len/8; k++)
            dft8_red_lazy_general_stride(p+k, len/8, node, F);
    }
    else if (depth == 4)
    {
        for (ulong k = 0; k < len/16; k++)
            dft16_red_lazy_general_stride(p+k, len/16, node, F);
    }
    //else if (depth == 5)
    //    dft32_red_lazy_general(p, node, F);
    else
    {
        // in: [0..4n), out: [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+len/2;
        const ulong w = F->tab_w[2*node];
        const ulong wpre = F->tab_w[2*node+1];
        ulong p_hi, p_lo, u, v;
        for (ulong k = 0; k < len/2; k+=4)
        {
            BUTTERFLY2_CT_LAZY(p0[k+0], p1[k+0], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            BUTTERFLY2_CT_LAZY(p0[k+1], p1[k+1], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            BUTTERFLY2_CT_LAZY(p0[k+2], p1[k+2], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
            BUTTERFLY2_CT_LAZY(p0[k+3], p1[k+3], F->mod, F->mod2, w, wpre, p_hi, p_lo, u, v);
        }
        _n_fft_red_rec2_lazy_general_stride(p0, len/2, depth-1, 2*node, F);
        _n_fft_red_rec2_lazy_general_stride(p1, len/2, depth-1, 2*node+1, F);
    }
}

// input [0..2n),  output [0..4n)
void _n_fft_red_rec2_lazy_stride(nn_ptr p, ulong len, ulong depth, n_fft_ctx_t F)
{
    // depth == 0: nothing to do
    if (depth == 1) // in [0..2n), out [0..4n)
    {
        ulong tmp;
        for (ulong k = 0; k < len/2; k++)
            DFT2_LAZY2(p[k], p[k+len/2], F->mod2, tmp);
    }
    else if (depth == 2) // in [0..2n), out [0..4n)
    {
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k++)
            DFT4_LAZY2_RED(p[k], p[len/4 + k], p[2*len/4 + k], p[3*len/4 + k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 3)
    {
        for (ulong k = 0; k < len/8; k++)
            dft8_red_lazy_stride(p+k, len/8, F);
    }
    else if (depth == 4)
    {
        for (ulong k = 0; k < len/16; k++)
            dft16_red_lazy_stride(p+k, len/16, F);  // in [0..2n), out [0..4n)
    }
    //else if (depth == 5)
    //    dft32_red_lazy(p, F);  // in [0..2n), out [0..4n)
    else
    {
        // input [0..2n) x [0..2n), output [0..2n) x [0..4n)
        // (general accepts [0..4n) as input for depth >= 3)
        ulong tmp;
        for (ulong k = 0; k < len/2; k++)
            DFT2_LAZY2_RED1(p[k+0], p[len/2+k+0], F->mod2, tmp);
        _n_fft_red_rec2_lazy_stride(p, len/2, depth-1, F);
        _n_fft_red_rec2_lazy_general_stride(p+len/2, len/2, depth-1, 1, F);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
