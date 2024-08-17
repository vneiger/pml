#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include "nmod_poly_integer_fft.h"
#include "nmod_integer_fft_base_cases.c"

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

/** 4-points DFT (DIF style):
 * (Decimation In Frequency returns bit reversed order evaluations)
 * returns (p(1), p(-1), p(I), p(-I)) where p(x) = a + b*x + c*x**2 + d*x**3
 * and I is typically a square root of -1 (this property is not exploited)
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 * */

// using NMOD_MUL_PRENORM (Inorm is I << mod.norm):
#define DFT4_DIF_PRENORM(a,b,c,d,Inorm,mod)             \
    do {                                                 \
        const ulong p0 = (a);                        \
        const ulong p1 = (b);                        \
        const ulong p2 = (c);                        \
        const ulong p3 = (d);                        \
        const ulong p4 = nmod_add(p0, p2, (mod));    \
        const ulong p5 = nmod_sub(p0, p2, (mod));    \
        const ulong p6 = nmod_add(p1, p3, (mod));    \
        ulong p7;                                    \
        NMOD_MUL_PRENORM(p7,                             \
                         nmod_sub(p1, p3, (mod)),        \
                         (Inorm), (mod));                \
        (a) = nmod_add(p4, p6, (mod));                   \
        (b) = nmod_sub(p4, p6, (mod));                   \
        (c) = nmod_add(p5, p7, (mod));                   \
        (d) = nmod_sub(p5, p7, (mod));                   \
    } while(0)

// using n_mulmod_shoup  (Ipre is I with Shoup's precomputation)
#define DFT4_DIF_SHOUP(a,b,c,d,I,Ipre,mod)              \
    do {                                                 \
        const ulong p0 = (a);                        \
        const ulong p1 = (b);                        \
        const ulong p2 = (c);                        \
        const ulong p3 = (d);                        \
        const ulong p4 = nmod_add(p0, p2, (mod));    \
        const ulong p5 = nmod_sub(p0, p2, (mod));    \
        const ulong p6 = nmod_add(p1, p3, (mod));    \
        const ulong p7 =                             \
                n_mulmod_shoup((I),                      \
                               nmod_sub(p1, p3, (mod)),  \
                               (Ipre), (mod).n);         \
        (a) = nmod_add(p4, p6, (mod));                   \
        (b) = nmod_sub(p4, p6, (mod));                   \
        (c) = nmod_add(p5, p7, (mod));                   \
        (d) = nmod_sub(p5, p7, (mod));                   \
    } while(0)

// lazy1: input in [0..n) --> output [0..4*n)
#define DFT4_DIF_SHOUP_LAZY1(a,b,c,d,I,Ipre,n)                     \
    do {                                                           \
        const ulong p0 = (a);                                  \
        const ulong p1 = (b);                                  \
        const ulong p2 = (c);                                  \
        const ulong p3 = (d);                                  \
        const ulong p4 = p0 + p2;              /* < 2*n */     \
        const ulong p5 = p0 + (n) - p2;        /* < 2*n */     \
        const ulong p6 = p1 + p3;              /* < 2*n */     \
        const ulong p7 =                       /* < 2*n */     \
             n_mulmod_shoup_lazy((I), p1 + (n) - p3, (Ipre), (n)); \
        (a) = p4 + p6;                             /* < 4*n */     \
        (b) = p4 + 2*(n) - p6;                     /* < 4*n */     \
        (c) = p5 + p7;                             /* < 4*n */     \
        (d) = p5 + 2*(n) - p7;                     /* < 4*n */     \
    } while(0)

// using n_mulmod_shoup  (Ipre is I with Shoup's precomputation)
// lazy red: input in [0..2*n) --> output [0..4*n)
// n2 is 2*n
#define DFT4_DIF_SHOUP_LAZY2_RED(a,b,c,d,I,Ipre,n,n2)                  \
    do {                                                               \
        const ulong p0 = (a);                                      \
        const ulong p1 = (b);                                      \
        const ulong p2 = (c);                                      \
        const ulong p3 = (d);                                      \
        ulong p4 = p0 + p2;                     /* < 4*n */        \
        if (p4 >= (n2))                                                \
            p4 -= (n2);                             /* < 2*n */        \
        ulong p5 = p0 + (n2) - p2;              /* < 4*n */        \
        if (p5 >= (n2))                                                \
            p5 -= (n2);                             /* < 2*n */        \
        ulong p6 = p1 + p3;                     /* < 4*n */        \
        if (p6 >= (n2))                                                \
            p6 -= (n2);                             /* < 2*n */        \
        const ulong p7 =                        /* < 2*n */        \
             n_mulmod_shoup_lazy((I), p1 + (n2) - p3, (Ipre), (n));    \
        (a) = p4 + p6;                              /* < 4*n */        \
        (b) = p4 + (n2) - p6;                       /* < 4*n */        \
        (c) = p5 + p7;                              /* < 4*n */        \
        (d) = p5 + (n2) - p7;                       /* < 4*n */        \
    } while(0)

// using n_mulmod_shoup  (Ipre is I with Shoup's precomputation)
// lazy red: input in [0..4*n) --> output [0..4*n)
// n4 is 4*n
#define DFT4_DIF_SHOUP_LAZY3_RED(a,b,c,d,I,Ipre,n)                     \
    do {                                                               \
        const ulong p0 = (a);                                      \
        const ulong p1 = (b);                                      \
        const ulong p2 = (c);                                      \
        const ulong p3 = (d);                                      \
        ulong p4 = p0 + p2;                    /* < 8*n */         \
        if (p4 >= 4*(n))                                               \
            p4 -= 4*(n);                             /* < 4*n */       \
        if (p4 >= 2*(n))                                               \
            p4 -= 2*(n);                             /* < 2*n */       \
        ulong p5 = p0 + 4*(n) - p2;            /* < 8*n */         \
        if (p5 >= 4*(n))                                               \
            p5 -= 4*(n);                             /* < 4*n */       \
        if (p5 >= 2*(n))                                               \
            p5 -= 2*(n);                             /* < 2*n */       \
        ulong p6 = p1 + p3;                    /* < 8*n */         \
        if (p6 >= 4*(n))                                               \
            p6 -= 4*(n);                             /* < 4*n */       \
        if (p6 >= 2*(n))                                               \
            p6 -= 2*(n);                             /* < 2*n */       \
        const ulong p7 =                       /* < 2*n */         \
             n_mulmod_shoup_lazy((I), p1 + 4*(n) - p3, (Ipre), (n));   \
        (a) = p4 + p6;                             /* < 4*n */         \
        (b) = p4 + 2*(n) - p6;                     /* < 4*n */         \
        (c) = p5 + p7;                             /* < 4*n */         \
        (d) = p5 + 2*(n) - p7;                     /* < 4*n */         \
    } while(0)

/*------------------------------------------------------------*/
/* 8-point DFT                                                */
/*------------------------------------------------------------*/

// remainder tree version
// lazy red: input in [0..2*n) --> output in [0..4*n)

/** one level of DFT, DIF style, unrolling 4 by 4:
 * does DFT2 and multiply by powers of w
 * */
#define DIF2_PRENORM_UNROLL4(p0,p1,ww,mod)             \
    do {                                               \
        const ulong u0 = (p0)[0];                  \
        const ulong u1 = (p0)[1];                  \
        const ulong u2 = (p0)[2];                  \
        const ulong u3 = (p0)[3];                  \
        ulong v0 = (p1)[0];                        \
        ulong v1 = (p1)[1];                        \
        ulong v2 = (p1)[2];                        \
        ulong v3 = (p1)[3];                        \
        (p0)[0] = nmod_add(u0, v0, (mod));             \
        (p0)[1] = nmod_add(u1, v1, (mod));             \
        (p0)[2] = nmod_add(u2, v2, (mod));             \
        (p0)[3] = nmod_add(u3, v3, (mod));             \
        v0 = nmod_sub(u0, v0, (mod)) << (mod).norm;    \
        v1 = nmod_sub(u1, v1, (mod)) << (mod).norm;    \
        v2 = nmod_sub(u2, v2, (mod)) << (mod).norm;    \
        v3 = nmod_sub(u3, v3, (mod)) << (mod).norm;    \
        NMOD_MUL_PRENORM((p1)[0], v0, (ww)[0], (mod)); \
        NMOD_MUL_PRENORM((p1)[1], v1, (ww)[1], (mod)); \
        NMOD_MUL_PRENORM((p1)[2], v2, (ww)[2], (mod)); \
        NMOD_MUL_PRENORM((p1)[3], v3, (ww)[3], (mod)); \
    } while (0)

#define DIF2_SHOUP_UNROLL4(p0,p1,ww,wwpre,mod)                      \
    do {                                                            \
        const ulong u0 = (p0)[0];                               \
        const ulong u1 = (p0)[1];                               \
        const ulong u2 = (p0)[2];                               \
        const ulong u3 = (p0)[3];                               \
        ulong v0 = (p1)[0];                                     \
        ulong v1 = (p1)[1];                                     \
        ulong v2 = (p1)[2];                                     \
        ulong v3 = (p1)[3];                                     \
        (p0)[0] = nmod_add(u0, v0, (mod));                          \
        (p0)[1] = nmod_add(u1, v1, (mod));                          \
        (p0)[2] = nmod_add(u2, v2, (mod));                          \
        (p0)[3] = nmod_add(u3, v3, (mod));                          \
        v0 = nmod_sub(u0, v0, (mod));                               \
        v1 = nmod_sub(u1, v1, (mod));                               \
        v2 = nmod_sub(u2, v2, (mod));                               \
        v3 = nmod_sub(u3, v3, (mod));                               \
        (p1)[0] = n_mulmod_shoup((ww)[0], v0, (wwpre)[0], (mod).n); \
        (p1)[1] = n_mulmod_shoup((ww)[1], v1, (wwpre)[1], (mod).n); \
        (p1)[2] = n_mulmod_shoup((ww)[2], v2, (wwpre)[2], (mod).n); \
        (p1)[3] = n_mulmod_shoup((ww)[3], v3, (wwpre)[3], (mod).n); \
    } while (0)



/*****************************************
*  TEMPORARIES: BENCH WITHOUT REDUCTION  *
*****************************************/

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

/** fft dif evaluation, in place
 * - len == 2**order, order >= 0
 * - input p has length at least len, seen as polynomial
 * - output p contains, in 0...len-1, the evaluations of trunc(p, len)
 *       at w**k, for 0 <= k < len, in bit reverse order
 * TODO : ASSUMES order suitable for w,
 *          and tables for 'order' already precomputed in F
 */
void _nmod_poly_dif_inplace_radix2_rec_prenorm(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_NMOD(p[0], p[1], F->mod);
    else if (order == 2)
        DFT4_DIF_PRENORM(p[0], p[1], p[2], p[3], F->tab_w[0][1] << F->mod.norm, F->mod);
    else
    {
        for (ulong k = 0; k < len/2; k++)
        {
            DFT2_NMOD(p[k], p[len/2+k], F->mod);
            NMOD_MUL_PRENORM(p[len/2+k], p[len/2+k] << F->mod.norm, F->tab_w[order-2][k], F->mod);
        }
        _nmod_poly_dif_inplace_radix2_rec_prenorm(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_prenorm(p+len/2, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_rec_prenorm_unroll4(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_NMOD(p[0], p[1], F->mod);
    else if (order == 2)
        DFT4_DIF_PRENORM(p[0], p[1], p[2], p[3], F->tab_w[0][1] << F->mod.norm, F->mod);
    else
    {
        // here order >= 3, len >= 8, unroll
        for (ulong k = 0; k < len/2; k+=4)
            DIF2_PRENORM_UNROLL4(p+k, p+len/2+k, F->tab_w[order-2]+k, F->mod);
        _nmod_poly_dif_inplace_radix2_rec_prenorm_unroll4(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_prenorm_unroll4(p+len/2, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_rec_shoup(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_NMOD(p[0], p[1], F->mod);
    else if (order == 2)
        DFT4_DIF_SHOUP(p[0], p[1], p[2], p[3], F->tab_w[0][1], F->tab_w_pre[0][1], F->mod);
    else
    {
        for (ulong k = 0; k < len/2; k++)
        {
            DFT2_NMOD(p[k], p[len/2+k], F->mod);
            p[len/2+k] = n_mulmod_shoup(F->tab_w[order-2][k], p[len/2+k], F->tab_w_pre[order-2][k], F->mod.n);
        }
        _nmod_poly_dif_inplace_radix2_rec_shoup(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_shoup(p+len/2, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_rec_shoup_unroll4(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_NMOD(p[0], p[1], F->mod);
    else if (order == 2)
        DFT4_DIF_SHOUP(p[0], p[1], p[2], p[3], F->tab_w[0][1], F->tab_w_pre[0][1], F->mod);
    else
    {
        // here order >= 3, len >= 8, try some unrolling
        for (ulong k = 0; k < len/2; k+=4)
            DIF2_SHOUP_UNROLL4(p+k, p+len/2+k, F->tab_w[order-2]+k, F->tab_w_pre[order-2]+k, F->mod);
        _nmod_poly_dif_inplace_radix2_rec_shoup_unroll4(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_shoup_unroll4(p+len/2, len/2, order-1, F);
    }

}

void _nmod_poly_dif_inplace_radix2_rec_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        DFT2_BENCH(p[0], p[1]);
    else if (order == 2)
        DFT4_DIF_BENCH(p[0], p[1], p[2], p[3], F->tab_w[0][1]);
    else
    {
        const nn_ptr ww = F->tab_w[order-2];
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
//            p1[k+0] = ww[k] * (u0 - v0) + u0 * v0 + v1 * F->mod.n;
//
//            p1[k+1] = ww[k+1] * (u1 - v1) + v1 * u1 + v2 * F->mod.n;
//
//            p1[k+2] = ww[k+2] * (u2 - v2) + v2 * u2 + v3 * (u2 - v2);
//
//            p1[k+3] = ww[k+3] * (u3 - v3) + u3 * v3 + u1 * F->mod.n;
//        }
//        _nmod_poly_dif_inplace_radix2_rec_bench(p, len/2, order-1, F);
//        _nmod_poly_dif_inplace_radix2_rec_bench(p+len/2, len/2, order-1, F);
//    }
//}

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

// iterative version, using NMOD_MUL_PRENORM
// TODO does not support order==1
void _nmod_poly_dif_inplace_radix2_iter_prenorm(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
        for (ulong k = 0; k < len; k+=llen)
            for (ulong kk = 0; kk < llen/2; kk+=4)
                DIF2_PRENORM_UNROLL4(p+k+kk, p+llen/2+k+kk, F->tab_w[ell-2]+kk, F->mod);
    // perform last two FFT layers
    for (ulong k = 0; k < len; k+=4)
        DFT4_DIF_PRENORM(p[k+0], p[k+1], p[k+2], p[k+3], F->tab_w[0][1] << F->mod.norm, F->mod);
}

// iterative version, using n_mulmod_shoup
// TODO does not support order==1
void _nmod_poly_dif_inplace_radix2_iter_shoup(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
        for (ulong k = 0; k < len; k+=llen)
            for (ulong kk = 0; kk < llen/2; kk+=4)
                DIF2_SHOUP_UNROLL4(p+k+kk, p+llen/2+k+kk, F->tab_w[ell-2]+kk, F->tab_w_pre[ell-2]+kk, F->mod);
    // perform last two FFT layers
    for (ulong k = 0; k < len; k+=4)
        DFT4_DIF_SHOUP(p[k+0], p[k+1], p[k+2], p[k+3], F->tab_w[0][1], F->tab_w_pre[0][1], F->mod);
}

void _nmod_poly_dif_inplace_radix2_iter_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell--, llen>>=1)
    {
        const nn_ptr ww = F->tab_w[ell-2];
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
        DFT4_DIF_BENCH(p[k+0], p[k+1], p[k+2], p[k+3], F->tab_w[0][1]);
}

/***************************
*  DIF radix 4 recursive  *
***************************/

// radix 4, rec, using prenorm
void _nmod_poly_dif_inplace_radix4_rec(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        ulong tmp = p[0];
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
        const ulong p0 = p[0];
        const ulong p1 = p[1];
        const ulong p2 = p[2];
        const ulong p3 = p[3];
        const ulong p4 = nmod_add(p0, p2, F->mod);
        const ulong p5 = nmod_sub(p0, p2, F->mod);
        const ulong p6 = nmod_add(p1, p3, F->mod);
        ulong p7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        p[0] = nmod_add(p4, p6, F->mod);
        p[1] = nmod_sub(p4, p6, F->mod);
        p[2] = nmod_add(p5, p7, F->mod);
        p[3] = nmod_sub(p5, p7, F->mod);
    }
    else if (order == 3)
    {
        const ulong p0 = p[0];
        const ulong q0 = p[1];
        const ulong p1 = p[2];
        const ulong q1 = p[3];
        const ulong p2 = p[4];
        const ulong q2 = p[5];
        const ulong p3 = p[6];
        const ulong q3 = p[7];
        const ulong p4 = nmod_add(p0, p2, F->mod);
        const ulong q4 = nmod_add(q0, q2, F->mod);
        const ulong p5 = nmod_sub(p0, p2, F->mod);
        const ulong q5 = nmod_sub(q0, q2, F->mod);
        const ulong p6 = nmod_add(p1, p3, F->mod);
        const ulong q6 = nmod_add(q1, q3, F->mod);
        ulong p7;
        ulong q7;
        NMOD_MUL_PRENORM(p7, nmod_sub(p1, p3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        NMOD_MUL_PRENORM(q7, nmod_sub(q1, q3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
        const ulong pp0 = nmod_add(p4, p6, F->mod);
        const ulong pp1 = nmod_add(q4, q6, F->mod);
        const ulong pp2 = nmod_sub(p4, p6, F->mod); // multiplied by F->tab_w[1][2*0] = 1
        ulong pp3;
        NMOD_MUL_PRENORM(pp3, nmod_sub(q4, q6, F->mod) << F->mod.norm, F->tab_w[1][2], F->mod);
        const ulong pp4 = nmod_add(p5, p7, F->mod); // multiplied by F->tab_w[1][0] = 1
        ulong pp5;
        NMOD_MUL_PRENORM(pp5, nmod_add(q5, q7, F->mod) << F->mod.norm, F->tab_w[1][1], F->mod);
        const ulong pp6 = nmod_sub(p5, p7, F->mod); // multiplied by F->tab_w[1][3*0] = 1
        ulong pp7;
        NMOD_MUL_PRENORM(pp7, nmod_sub(q5, q7, F->mod) << F->mod.norm, F->tab_w[1][3], F->mod);
        p[0] = nmod_add(pp0, pp1, F->mod);
        p[1] = nmod_sub(pp0, pp1, F->mod);
        p[2] = nmod_add(pp2, pp3, F->mod);
        p[3] = nmod_sub(pp2, pp3, F->mod);
        p[4] = nmod_add(pp4, pp5, F->mod);
        p[5] = nmod_sub(pp4, pp5, F->mod);
        p[6] = nmod_add(pp6, pp7, F->mod);
        p[7] = nmod_sub(pp6, pp7, F->mod);
    }
    else
    {
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+(len/4);
        const nn_ptr p2 = p1+(len/4);
        const nn_ptr p3 = p2+(len/4);
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
            const ulong u0 = p0[k];
            const ulong u1 = p1[k];
            const ulong u2 = p2[k];
            const ulong u3 = p3[k];
            const ulong u4 = nmod_add(u0, u2, F->mod);
            const ulong u5 = nmod_sub(u0, u2, F->mod);
            const ulong u6 = nmod_add(u1, u3, F->mod);
            ulong u7;
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


// radix 4, rec, bench
void _nmod_poly_dif_inplace_radix4_rec_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        //                  [1   1]
        // p <- [p[0] p[1]] [1  -1]
        ulong tmp = p[0];
        p[0] = p[0] + p[1];
        p[1] = tmp + p[1];
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
        const ulong p0 = p[0];
        const ulong p1 = p[1];
        const ulong p2 = p[2];
        const ulong p3 = p[3];
        const ulong p4 = p0 + p2;
        const ulong p5 = p0 - p2;
        const ulong p6 = p1 + p3;
        const ulong p7 = F->tab_w[0][1] * (p1 - p3);
        p[0] = p4 + p6;
        p[1] = p4 - p6;
        p[2] = p5 + p7;
        p[3] = p5 - p7;
    }
    else if (order == 3)
    {
        const ulong p0 = p[0];
        const ulong q0 = p[1];
        const ulong p1 = p[2];
        const ulong q1 = p[3];
        const ulong p2 = p[4];
        const ulong q2 = p[5];
        const ulong p3 = p[6];
        const ulong q3 = p[7];
        const ulong p7 = F->tab_w[0][1] * (p1 - p3);
        const ulong q7 = F->tab_w[0][1] * (q1 - q3);
        const ulong pp0 = p0 + p1 + p2 + p3;
        const ulong pp1 = q0 + q1 + q2 + q3;
        const ulong pp2 = p0 + p2 - (p1 + p3);
        const ulong pp3 = F->tab_w[1][2] * (q0 + q2 - (q1 + q3));
        const ulong pp4 = p0 - p2 + p7;
        const ulong pp5 = F->tab_w[1][1] * (q0 - q2 + q7);
        const ulong pp6 = p0 - p2 - p7;
        const ulong pp7 = F->tab_w[1][3] * (q0 - q2 - q7);
        p[0] = pp0 + pp1;
        p[1] = pp0 - pp1;
        p[2] = pp2 + pp3;
        p[3] = pp2 - pp3;
        p[4] = pp4 + pp5;
        p[5] = pp4 - pp5;
        p[6] = pp6 + pp7;
        p[7] = pp6 - pp7;
    }
    //{
    //    ulong p0 = p[0];
    //    ulong q0 = p[1];
    //    ulong p1 = p[2];
    //    ulong q1 = p[3];
    //    ulong p2 = p[4];
    //    ulong q2 = p[5];
    //    ulong p3 = p[6];
    //    ulong q3 = p[7];
    //    const ulong p4 = p0 + p2;
    //    const ulong q4 = q0 + q2;
    //    const ulong p5 = p0 - p2;
    //    const ulong q5 = q0 - q2;
    //    const ulong p6 = p1 + p3;
    //    const ulong q6 = q1 + q3;
    //    const ulong p7 = F->tab_w[0][1] * (p1 - p3);
    //    const ulong q7 = F->tab_w[0][1] * (q1 - q3);
    //    p0 = p4 + p6;
    //    p1 = q4 + q6;
    //    p2 = p4 - p6; // multiplied by F->tab_w[1][2*0] = 1
    //    p3 = F->tab_w[1][2] * (q4 - q6);
    //    q0 = p5 + p7; // multiplied by F->tab_w[1][0] = 1
    //    q1 = F->tab_w[1][1] * (q5 + q7);
    //    q2 = p5 - p7; // multiplied by F->tab_w[1][3*0] = 1
    //    q3 = F->tab_w[1][3] * (q5 - q7);
    //    p[0] = p0 + p1;
    //    p[1] = p0 - p1;
    //    p[2] = p2 + p3;
    //    p[3] = p2 - p3;
    //    p[4] = q4 + q5;
    //    p[5] = q4 - q5;
    //    p[6] = q6 + q7;
    //    p[7] = q6 - q7;
    //}
    else
    {
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+(len/4);
        const nn_ptr p2 = p1+(len/4);
        const nn_ptr p3 = p2+(len/4);
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
            const ulong u0 = p0[k];
            const ulong u1 = p1[k];
            const ulong u2 = p2[k];
            const ulong u3 = p3[k];
            const ulong u4 = u0 + u2;
            const ulong u5 = u0 - u2;
            const ulong u6 = u1 + u3;
            const ulong u7 = F->tab_w[0][1] * (u1 - u3);
            p0[k] = u4 + u6;
            p1[k] = F->tab_w[order-2][2*k] * (u4 - u6);
            p2[k] = F->tab_w[order-2][k] * (u5 + u7);
            if (3*k < len/2)
                p3[k] = F->tab_w[order-2][3*k] * (u5 - u7);
            else
                p3[k] = (F->mod.n - F->tab_w[order-2][3*k-len/2]) * (u5 - u7);
        }
        _nmod_poly_dif_inplace_radix4_rec_bench(p0, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec_bench(p1, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec_bench(p2, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec_bench(p3, len/4, order-2, F);
    }
}

/***************************
*  DIF radix 4 iterative  *
***************************/

// radix 4, iter, using prenorm
void _nmod_poly_dif_inplace_radix4_iter(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // perform FFT layers up to order 2
    ulong llen = len;
    for (ulong ell = order; ell > 2; ell-=2, llen>>=2)
    {
        const nn_ptr ww = F->tab_w[ell-2];
        for (ulong k = 0; k < len; k+=llen)
        {
            const nn_ptr p0 = p+k;
            const nn_ptr p1 = p0+(llen/4);
            const nn_ptr p2 = p1+(llen/4);
            const nn_ptr p3 = p2+(llen/4);
            for (ulong kk = 0; kk < llen/4; kk++)
            {
                const ulong u0 = p0[kk];
                const ulong u1 = p1[kk];
                const ulong u2 = p2[kk];
                const ulong u3 = p3[kk];
                const ulong u4 = nmod_add(u0, u2, F->mod);
                const ulong u5 = nmod_sub(u0, u2, F->mod);
                const ulong u6 = nmod_add(u1, u3, F->mod);
                ulong u7;
                NMOD_MUL_PRENORM(u7, nmod_sub(u1, u3, F->mod) << F->mod.norm, F->tab_w[0][1], F->mod);
                p0[kk] = nmod_add(u4, u6, F->mod);
                NMOD_MUL_PRENORM(p1[kk], nmod_sub(u4, u6, F->mod) << F->mod.norm, ww[2*kk], F->mod);
                NMOD_MUL_PRENORM(p2[kk], nmod_add(u5, u7, F->mod) << F->mod.norm, ww[kk], F->mod);
                if (3*kk < llen/2)
                    NMOD_MUL_PRENORM(p3[kk], nmod_sub(u5, u7, F->mod) << F->mod.norm, ww[3*kk], F->mod);
                else
                    NMOD_MUL_PRENORM(p3[kk], nmod_sub(u5, u7, F->mod) << F->mod.norm, F->mod.n - ww[3*kk-llen/2], F->mod);
            }
        }
    }
    // perform last two FFT layers
    if (llen == 4)
        for (ulong k = 0; k < len; k+=4)
        {
            const ulong p0 = p[k+0];
            const ulong p1 = p[k+1];
            const ulong p2 = p[k+2];
            const ulong p3 = p[k+3];
            const ulong p4 = nmod_add(p0, p2, F->mod);
            const ulong p5 = nmod_sub(p0, p2, F->mod);
            const ulong p6 = nmod_add(p1, p3, F->mod);
            const ulong tmp7 = nmod_sub(p1, p3, F->mod);
            ulong p7;
            NMOD_MUL_PRENORM(p7, tmp7 << F->mod.norm, F->tab_w[0][1], F->mod);
            p[k+0] = nmod_add(p4, p6, F->mod);
            p[k+1] = nmod_sub(p4, p6, F->mod);
            p[k+2] = nmod_add(p5, p7, F->mod);
            p[k+3] = nmod_sub(p5, p7, F->mod);
        }
    else // llen == 2
        for (ulong k = 0; k < len; k+=2)
        {
            const ulong tmp = p[k+0];
            p[k+0] = nmod_add(p[k+0], p[k+1], F->mod);
            p[k+1] = nmod_sub(tmp, p[k+1], F->mod);
        }
}


/*******************************
*  Reduction tree viewpoint   *
*******************************/

/** fft evaluation, in place
 * - len == 2**order, order >= 0
 * - input p has length at least len, seen as polynomial
 * - output p contains, in 0...len-1, the evaluations of trunc(p, len)
 *       at w**k, for 0 <= k < len, in bit reverse order
 * TODO : ASSUMES order suitable for w,
 *          and RED tables for 'order' already precomputed in F
 */
// note: unroll4 does not seem to gain anything
void _nmod_poly_red_inplace_radix2_rec_prenorm(nn_ptr p, ulong len, ulong order, ulong node, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        if (node==0)  // w == 1
            DFT2_NMOD(p[0], p[1], F->mod);
        else
        {
            const ulong w = F->tab_w[1][node];
            const ulong u = p[0];
            ulong v = p[1];
            NMOD_MUL_PRENORM(v, v << F->mod.norm, w, F->mod);
            p[0] = nmod_add(u, v, F->mod);
            p[1] = nmod_sub(u, v, F->mod);
        }
    }
    else if (order == 2)
    {
        if (node==0)  // w == 1
        {
            const ulong u0 = p[0];
            const ulong u1 = p[1];
            ulong v0 = p[2];
            ulong v1 = p[3];

            const ulong p0 = nmod_add(u0, v0, F->mod);
            const ulong p1 = nmod_add(u1, v1, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);
            v1 = nmod_sub(u1, v1, F->mod);

            NMOD_MUL_PRENORM(v1, v1 << F->mod.norm, F->tab_w[1][1], F->mod);
            p[0] = nmod_add(p0, p1, F->mod);
            p[1] = nmod_sub(p0, p1, F->mod);
            p[2] = nmod_add(v0, v1, F->mod);
            p[3] = nmod_sub(v0, v1, F->mod);
        }
        else
        {
            const ulong u0 = p[0];
            const ulong u1 = p[1];
            ulong v0 = p[2];
            ulong v1 = p[3];

            NMOD_MUL_PRENORM(v0, v0 << F->mod.norm, F->tab_w[1][node], F->mod);
            NMOD_MUL_PRENORM(v1, v1 << F->mod.norm, F->tab_w[1][node], F->mod);
            const ulong p0 = nmod_add(u0, v0, F->mod);
            ulong p1 = nmod_add(u1, v1, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);  // p2
            v1 = nmod_sub(u1, v1, F->mod);  // p3

            NMOD_MUL_PRENORM(p1, p1 << F->mod.norm, F->tab_w[1][2*node], F->mod);
            NMOD_MUL_PRENORM(v1, v1 << F->mod.norm, F->tab_w[1][2*node+1], F->mod);
            p[0] = nmod_add(p0, p1, F->mod);
            p[1] = nmod_sub(p0, p1, F->mod);
            p[2] = nmod_add(v0, v1, F->mod);
            p[3] = nmod_sub(v0, v1, F->mod);
        }
    }
    else if (order == 3)
    {
        if (node==0)  // w == 1
        {
            ulong u0 = p[0];
            ulong u1 = p[1];
            ulong u2 = p[2];
            ulong u3 = p[3];
            ulong v0 = p[4];
            ulong v1 = p[5];
            ulong v2 = p[6];
            ulong v3 = p[7];

            // mod x**4 - 1 | x**4 + 1
            ulong p0 = nmod_add(u0, v0, F->mod);
            ulong p1 = nmod_add(u1, v1, F->mod);
            ulong p2 = nmod_add(u2, v2, F->mod);
            ulong p3 = nmod_add(u3, v3, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);
            v1 = nmod_sub(u1, v1, F->mod);
            v2 = nmod_sub(u2, v2, F->mod);
            v3 = nmod_sub(u3, v3, F->mod);

            // left, mod x**2 - 1 | x**2 + 1
            u0 = nmod_add(p0, p2, F->mod);
            u1 = nmod_add(p1, p3, F->mod);
            u2 = nmod_sub(p0, p2, F->mod);
            u3 = nmod_sub(p1, p3, F->mod);

            // left-left, mod x-1 | x+1
            p0 = nmod_add(u0, u1, F->mod);
            p1 = nmod_sub(u0, u1, F->mod);

            // left-right, mod x-I | x+I
            NMOD_MUL_PRENORM(u3, u3 << F->mod.norm, F->tab_w[1][1], F->mod);
            p2 = nmod_add(u2, u3, F->mod);
            p3 = nmod_sub(u2, u3, F->mod);

            // right, mod x**2 - I | x**2 + I
            NMOD_MUL_PRENORM(v2, v2 << F->mod.norm, F->tab_w[1][1], F->mod);
            NMOD_MUL_PRENORM(v3, v3 << F->mod.norm, F->tab_w[1][1], F->mod);
            u0 = nmod_add(v0, v2, F->mod);
            u1 = nmod_add(v1, v3, F->mod);
            u2 = nmod_sub(v0, v2, F->mod);
            u3 = nmod_sub(v1, v3, F->mod);

            // right-left, mod x - J | x + J
            NMOD_MUL_PRENORM(u1, u1 << F->mod.norm, F->tab_w[1][2], F->mod);
            v0 = nmod_add(u0, u1, F->mod);
            v1 = nmod_sub(u0, u1, F->mod);

            // right-right, mod x - I*J | x + I*J
            NMOD_MUL_PRENORM(u3, u3 << F->mod.norm, F->tab_w[1][3], F->mod);
            v2 = nmod_add(u2, u3, F->mod);
            v3 = nmod_sub(u2, u3, F->mod);

            p[0] = p0;
            p[1] = p1;
            p[2] = p2;
            p[3] = p3;
            p[4] = v0;
            p[5] = v1;
            p[6] = v2;
            p[7] = v3;
        }
        else
        {
            ulong u0 = p[0];
            ulong u1 = p[1];
            ulong u2 = p[2];
            ulong u3 = p[3];
            ulong v0 = p[4];
            ulong v1 = p[5];
            ulong v2 = p[6];
            ulong v3 = p[7];

            // w = F->tab_w[1][node]
            ulong w = F->tab_w[1][node] << F->mod.norm;

            // mod x**4 - w | x**4 + w
            NMOD_MUL_PRENORM(v0, v0, w, F->mod);
            NMOD_MUL_PRENORM(v1, v1, w, F->mod);
            NMOD_MUL_PRENORM(v2, v2, w, F->mod);
            NMOD_MUL_PRENORM(v3, v3, w, F->mod);
            ulong p0 = nmod_add(u0, v0, F->mod);
            ulong p1 = nmod_add(u1, v1, F->mod);
            ulong p2 = nmod_add(u2, v2, F->mod);
            ulong p3 = nmod_add(u3, v3, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);
            v1 = nmod_sub(u1, v1, F->mod);
            v2 = nmod_sub(u2, v2, F->mod);
            v3 = nmod_sub(u3, v3, F->mod);

            // left, mod x**2 - sqrt(w) | x**2 + sqrt(w)
            w = F->tab_w[1][2*node] << F->mod.norm;
            NMOD_MUL_PRENORM(p2, p2, w, F->mod);
            NMOD_MUL_PRENORM(p3, p3, w, F->mod);
            u0 = nmod_add(p0, p2, F->mod);
            u1 = nmod_add(p1, p3, F->mod);
            u2 = nmod_sub(p0, p2, F->mod);
            u3 = nmod_sub(p1, p3, F->mod);

            // left-left, mod x - fort(w) | x + fort(w)
            NMOD_MUL_PRENORM(u1, u1, F->tab_w[1][4*node] << F->mod.norm, F->mod);
            p0 = nmod_add(u0, u1, F->mod);
            p1 = nmod_sub(u0, u1, F->mod);

            // left-right, mod x - I*fort(w) | x+ I*fort(w)
            NMOD_MUL_PRENORM(u3, u3 << F->mod.norm, F->tab_w[1][4*node+1], F->mod);
            p2 = nmod_add(u2, u3, F->mod);
            p3 = nmod_sub(u2, u3, F->mod);

            // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
            w = F->tab_w[1][2*node+1] << F->mod.norm;
            NMOD_MUL_PRENORM(v2, v2, w, F->mod);
            NMOD_MUL_PRENORM(v3, v3, w, F->mod);
            u0 = nmod_add(v0, v2, F->mod);
            u1 = nmod_add(v1, v3, F->mod);
            u2 = nmod_sub(v0, v2, F->mod);
            u3 = nmod_sub(v1, v3, F->mod);

            // right-left, mod x - J*fort(w) | x + J*fort(w)
            NMOD_MUL_PRENORM(u1, u1 << F->mod.norm, F->tab_w[1][4*node+2], F->mod);
            v0 = nmod_add(u0, u1, F->mod);
            v1 = nmod_sub(u0, u1, F->mod);

            // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
            NMOD_MUL_PRENORM(u3, u3 << F->mod.norm, F->tab_w[1][4*node+3], F->mod);
            v2 = nmod_add(u2, u3, F->mod);
            v3 = nmod_sub(u2, u3, F->mod);

            p[0] = p0;
            p[1] = p1;
            p[2] = p2;
            p[3] = p3;
            p[4] = v0;
            p[5] = v1;
            p[6] = v2;
            p[7] = v3;
        }
    }
    else
    {
        if (node==0)  // w == 1
            for (ulong k = 0; k < len/2; k++)
                DFT2_NMOD(p[k], p[len/2+k], F->mod);
        else
        {
            const ulong w = F->tab_w[1][node];
            for (ulong k = 0; k < len/2; k++)
            {
                const ulong u = p[k];
                ulong v = p[len/2+k];
                NMOD_MUL_PRENORM(v, v << F->mod.norm, w, F->mod);
                p[k] = nmod_add(u, v, F->mod);
                p[len/2+k] = nmod_sub(u, v, F->mod);
            }
        }
        _nmod_poly_red_inplace_radix2_rec_prenorm(p, len/2, order-1, 2*node, F);
        _nmod_poly_red_inplace_radix2_rec_prenorm(p+len/2, len/2, order-1, 2*node+1, F);
    }
}

// note: unroll4 does not seem to gain anything
void _nmod_poly_red_inplace_radix2_rec_shoup(nn_ptr p, ulong len, ulong order, ulong node, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        if (node==0)  // w == 1
            DFT2_NMOD(p[0], p[1], F->mod);
        else
        {
            const ulong w = F->tab_w[1][node];
            const ulong wpre = F->tab_w_pre[1][node];
            const ulong u = p[0];
            ulong v = p[1];
            v = n_mulmod_shoup(w, v, wpre, F->mod.n);
            p[0] = nmod_add(u, v, F->mod);
            p[1] = nmod_sub(u, v, F->mod);
        }
    }
    else if (order == 2)
    {
        if (node==0)  // w == 1
        {
            const ulong u0 = p[0];
            const ulong u1 = p[1];
            ulong v0 = p[2];
            ulong v1 = p[3];

            const ulong p0 = nmod_add(u0, v0, F->mod);
            const ulong p1 = nmod_add(u1, v1, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);
            v1 = nmod_sub(u1, v1, F->mod);

            const ulong w = F->tab_w[1][1];
            const ulong wpre = F->tab_w_pre[1][1];
            v1 = n_mulmod_shoup(w, v1, wpre, F->mod.n);
            p[0] = nmod_add(p0, p1, F->mod);
            p[1] = nmod_sub(p0, p1, F->mod);
            p[2] = nmod_add(v0, v1, F->mod);
            p[3] = nmod_sub(v0, v1, F->mod);
        }
        else
        {
            ulong w = F->tab_w[1][node];
            ulong wpre = F->tab_w_pre[1][node];
            const ulong u0 = p[0];
            const ulong u1 = p[1];
            ulong v0 = p[2];
            ulong v1 = p[3];

            v0 = n_mulmod_shoup(w, v0, wpre, F->mod.n);
            v1 = n_mulmod_shoup(w, v1, wpre, F->mod.n);
            const ulong p0 = nmod_add(u0, v0, F->mod);
            ulong p1 = nmod_add(u1, v1, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);  // p2
            v1 = nmod_sub(u1, v1, F->mod);  // p3

            w = F->tab_w[1][2*node];
            wpre = F->tab_w_pre[1][2*node];
            p1 = n_mulmod_shoup(w, p1, wpre, F->mod.n);
            w = F->tab_w[1][2*node+1];
            wpre = F->tab_w_pre[1][2*node+1];
            v1 = n_mulmod_shoup(w, v1, wpre, F->mod.n);
            p[0] = nmod_add(p0, p1, F->mod);
            p[1] = nmod_sub(p0, p1, F->mod);
            p[2] = nmod_add(v0, v1, F->mod);
            p[3] = nmod_sub(v0, v1, F->mod);
        }
    }
    else if (order == 3)
    {
        if (node==0)  // w == 1
        {
            ulong u0 = p[0];
            ulong u1 = p[1];
            ulong u2 = p[2];
            ulong u3 = p[3];
            ulong v0 = p[4];
            ulong v1 = p[5];
            ulong v2 = p[6];
            ulong v3 = p[7];

            // mod x**4 - 1 | x**4 + 1
            ulong p0 = nmod_add(u0, v0, F->mod);
            ulong p1 = nmod_add(u1, v1, F->mod);
            ulong p2 = nmod_add(u2, v2, F->mod);
            ulong p3 = nmod_add(u3, v3, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);
            v1 = nmod_sub(u1, v1, F->mod);
            v2 = nmod_sub(u2, v2, F->mod);
            v3 = nmod_sub(u3, v3, F->mod);

            // left, mod x**2 - 1 | x**2 + 1
            u0 = nmod_add(p0, p2, F->mod);
            u1 = nmod_add(p1, p3, F->mod);
            u2 = nmod_sub(p0, p2, F->mod);
            u3 = nmod_sub(p1, p3, F->mod);

            // left-left, mod x-1 | x+1
            p0 = nmod_add(u0, u1, F->mod);
            p1 = nmod_sub(u0, u1, F->mod);

            // left-right, mod x-I | x+I
            u3 = n_mulmod_shoup(F->tab_w[1][1], u3, F->tab_w_pre[1][1], F->mod.n);
            p2 = nmod_add(u2, u3, F->mod);
            p3 = nmod_sub(u2, u3, F->mod);

            // right, mod x**2 - I | x**2 + I
            v2 = n_mulmod_shoup(F->tab_w[1][1], v2, F->tab_w_pre[1][1], F->mod.n);
            v3 = n_mulmod_shoup(F->tab_w[1][1], v3, F->tab_w_pre[1][1], F->mod.n);
            u0 = nmod_add(v0, v2, F->mod);
            u1 = nmod_add(v1, v3, F->mod);
            u2 = nmod_sub(v0, v2, F->mod);
            u3 = nmod_sub(v1, v3, F->mod);

            // right-left, mod x - J | x + J
            u1 = n_mulmod_shoup(F->tab_w[1][2], u1, F->tab_w_pre[1][2], F->mod.n);
            v0 = nmod_add(u0, u1, F->mod);
            v1 = nmod_sub(u0, u1, F->mod);

            // right-right, mod x - I*J | x + I*J
            u3 = n_mulmod_shoup(F->tab_w[1][3], u3, F->tab_w_pre[1][3], F->mod.n);
            v2 = nmod_add(u2, u3, F->mod);
            v3 = nmod_sub(u2, u3, F->mod);

            p[0] = p0;
            p[1] = p1;
            p[2] = p2;
            p[3] = p3;
            p[4] = v0;
            p[5] = v1;
            p[6] = v2;
            p[7] = v3;
        }
        else
        {
            ulong u0 = p[0];
            ulong u1 = p[1];
            ulong u2 = p[2];
            ulong u3 = p[3];
            ulong v0 = p[4];
            ulong v1 = p[5];
            ulong v2 = p[6];
            ulong v3 = p[7];

            // mod x**4 - w | x**4 + w
            ulong w = F->tab_w[1][node];
            ulong wpre = F->tab_w_pre[1][node];
            v0 = n_mulmod_shoup(w, v0, wpre, F->mod.n);
            v1 = n_mulmod_shoup(w, v1, wpre, F->mod.n);
            v2 = n_mulmod_shoup(w, v2, wpre, F->mod.n);
            v3 = n_mulmod_shoup(w, v3, wpre, F->mod.n);
            ulong p0 = nmod_add(u0, v0, F->mod);
            ulong p1 = nmod_add(u1, v1, F->mod);
            ulong p2 = nmod_add(u2, v2, F->mod);
            ulong p3 = nmod_add(u3, v3, F->mod);
            v0 = nmod_sub(u0, v0, F->mod);
            v1 = nmod_sub(u1, v1, F->mod);
            v2 = nmod_sub(u2, v2, F->mod);
            v3 = nmod_sub(u3, v3, F->mod);

            // left, mod x**2 - sqrt(w) | x**2 + sqrt(w)
            w = F->tab_w[1][2*node];
            wpre = F->tab_w_pre[1][2*node];
            p2 = n_mulmod_shoup(w, p2, wpre, F->mod.n);
            p3 = n_mulmod_shoup(w, p3, wpre, F->mod.n);
            u0 = nmod_add(p0, p2, F->mod);
            u1 = nmod_add(p1, p3, F->mod);
            u2 = nmod_sub(p0, p2, F->mod);
            u3 = nmod_sub(p1, p3, F->mod);

            // left-left, mod x - fort(w) | x + fort(w)
            u1 = n_mulmod_shoup(F->tab_w[1][4*node], u1, F->tab_w_pre[1][4*node], F->mod.n);
            p0 = nmod_add(u0, u1, F->mod);
            p1 = nmod_sub(u0, u1, F->mod);

            // left-right, mod x - I*fort(w) | x+ I*fort(w)
            u3 = n_mulmod_shoup(F->tab_w[1][4*node+1], u3, F->tab_w_pre[1][4*node+1], F->mod.n);
            p2 = nmod_add(u2, u3, F->mod);
            p3 = nmod_sub(u2, u3, F->mod);

            // right, mod x**2 - I*sqrt(w) | x**2 + I*sqrt(w)
            w = F->tab_w[1][2*node+1];
            wpre = F->tab_w_pre[1][2*node+1];
            v2 = n_mulmod_shoup(w, v2, wpre, F->mod.n);
            v3 = n_mulmod_shoup(w, v3, wpre, F->mod.n);
            u0 = nmod_add(v0, v2, F->mod);
            u1 = nmod_add(v1, v3, F->mod);
            u2 = nmod_sub(v0, v2, F->mod);
            u3 = nmod_sub(v1, v3, F->mod);

            // right-left, mod x - J*fort(w) | x + J*fort(w)
            u1 = n_mulmod_shoup(F->tab_w[1][4*node+2], u1, F->tab_w_pre[1][4*node+2], F->mod.n);
            v0 = nmod_add(u0, u1, F->mod);
            v1 = nmod_sub(u0, u1, F->mod);

            // right-right, mod x - I*J*fort(w) | x + I*J*fort(w)
            u3 = n_mulmod_shoup(F->tab_w[1][4*node+3], u3, F->tab_w_pre[1][4*node+3], F->mod.n);
            v2 = nmod_add(u2, u3, F->mod);
            v3 = nmod_sub(u2, u3, F->mod);

            p[0] = p0;
            p[1] = p1;
            p[2] = p2;
            p[3] = p3;
            p[4] = v0;
            p[5] = v1;
            p[6] = v2;
            p[7] = v3;
        }
    }
    else
    {
        if (node==0)  // w == 1
            for (ulong k = 0; k < len/2; k++)
                DFT2_NMOD(p[k], p[len/2+k], F->mod);
        else
        {
            const ulong w = F->tab_w[1][node];
            const ulong wpre = F->tab_w_pre[1][node];
            for (ulong k = 0; k < len/2; k++)
            {
                const ulong u = p[k];
                ulong v = p[len/2+k];
                v = n_mulmod_shoup(w, v, wpre, F->mod.n);
                p[k] = nmod_add(u, v, F->mod);
                p[len/2+k] = nmod_sub(u, v, F->mod);
            }
        }
        _nmod_poly_red_inplace_radix2_rec_shoup(p, len/2, order-1, 2*node, F);
        _nmod_poly_red_inplace_radix2_rec_shoup(p+len/2, len/2, order-1, 2*node+1, F);
    }
}


/**********
*  Lazy  *
**********/

// input [0..2*n), output [0..4*n)
void _nmod_poly_dif_inplace_radix2_rec_shoup_lazy(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
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
        _nmod_poly_dif_inplace_radix2_rec_shoup_lazy(p, len/2, order-1, F);
        _nmod_poly_dif_inplace_radix2_rec_shoup_lazy(p+len/2, len/2, order-1, F);
    }
}

void _nmod_poly_dif_inplace_radix2_iter_shoup_lazy(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
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

void _nmod_poly_red_inplace_radix2_rec_shoup_lazy(nn_ptr p, ulong len, ulong order, ulong node, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
    {
        if (node==0)  // w == 1
            DFT2_LAZY3_RED(p[0], p[1], F->modn4);
        else
        {
            // in [0..2n), out [0..4n)
            const ulong u = p[0];
            const ulong v = n_mulmod_shoup_lazy(F->tab_w[1][node], p[1], F->tab_w_pre[1][node], F->mod.n);
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
            const ulong u0 = p[0];
            const ulong u1 = p[1];
            ulong v0 = p[2];
            ulong v1 = p[3];

            ulong w = F->tab_w[1][node];
            ulong wpre = F->tab_w_pre[1][node];
            v0 = n_mulmod_shoup_lazy(w, v0, wpre, F->mod.n);
            v1 = n_mulmod_shoup_lazy(w, v1, wpre, F->mod.n);
            ulong p0 = u0 + v0;  // [0..4n)
            if (p0 >= F->modn2)
                p0 -= F->modn2;      // [0..2n)
            ulong p1 = u1 + v1;  // [0..4n)
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
            ulong u0 = p[0];
            ulong u1 = p[1];
            ulong u2 = p[2];
            ulong u3 = p[3];
            ulong v0 = p[4];
            ulong v1 = p[5];
            ulong v2 = p[6];
            ulong v3 = p[7];

            // mod x**4 - w | x**4 + w
            ulong w = F->tab_w[1][node];
            ulong wpre = F->tab_w_pre[1][node];
            v0 = n_mulmod_shoup_lazy(w, v0, wpre, F->mod.n);
            v1 = n_mulmod_shoup_lazy(w, v1, wpre, F->mod.n);
            v2 = n_mulmod_shoup_lazy(w, v2, wpre, F->mod.n);
            v3 = n_mulmod_shoup_lazy(w, v3, wpre, F->mod.n);
            ulong p0 = u0 + v0;   // [0..4n)
            ulong p1 = u1 + v1;   // [0..4n)
            ulong p2 = u2 + v2;   // [0..4n)
            ulong p3 = u3 + v3;   // [0..4n)
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
            const ulong w = F->tab_w[1][node];
            const ulong wpre = F->tab_w_pre[1][node];
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
        _nmod_poly_red_inplace_radix2_rec_shoup_lazy(p, len/2, order-1, 2*node, F);
        _nmod_poly_red_inplace_radix2_rec_shoup_lazy(p+len/2, len/2, order-1, 2*node+1, F);
    }
}

// in [0..2n) out [0..4n)
void _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F)
{
    // order == 0: nothing to do
    if (order == 1)
        // in [0..2n) out [0..4n)
        DFT2_LAZY2_RED(p[0], p[1], F->modn2);
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
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+(len/4);
        const nn_ptr p2 = p1+(len/4);
        const nn_ptr p3 = p2+(len/4);
        for (ulong k = 0; k < len/4; k++)
        {
            const ulong u0 = p0[k];
            const ulong u1 = p1[k];
            const ulong u2 = p2[k];
            const ulong u3 = p3[k];
            ulong u4 = u0 + u2;  // [0..4n)
            ulong u5 = u0 + F->modn2 - u2;  // [0..4n)
            ulong u6 = u1 + u3;  // [0..4n)
            ulong u7 = u1 + F->modn2 - u3;

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
        _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(p0, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(p1, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(p2, len/4, order-2, F);
        _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(p3, len/4, order-2, F);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
