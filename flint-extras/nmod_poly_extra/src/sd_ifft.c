/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <flint/fft_small.h>

/*
    N is supposed to be a good fit for the number of points to process per loop
    in the radix 4 inverse butterflies.
*/

#define N 8
#define VECND vec8d
#define VECNOP(op) vec8d_##op


/************************** inverse butterfly *********************************
    2*a0 =      (b0 + b1)
    2*a1 = w^-1*(b0 - b1)
    W := -w^-1
*/
#define RADIX_2_REVERSE_PARAM_J_IS_Z(V, Q) \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_REVERSE_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    y0 = V##_add(x0, x1); \
    y1 = V##_sub(x0, x1); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_reduce_to_pm1n(y1, n, ninv); \
    V##_store(X0, y0); \
    V##_store(X1, y1); \
}

#define RADIX_2_REVERSE_PARAM_J_IS_NZ(V, Q, j_mr, j_bits) \
    V W = V##_set_d(Q->w2tab[j_bits][j_mr]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_REVERSE_MOTH_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1, y0, y1; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    y0 = V##_add(x0, x1); \
    y1 = V##_sub(x1, x0); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_mulmod(y1, W, n, ninv); \
    V##_store(X0, y0); \
    V##_store(X1, y1); \
}

#define _RADIX_2_REVERSE_PARAM_J_IS_Z(...)  RADIX_2_REVERSE_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_REVERSE_MOTH_J_IS_Z(...)   RADIX_2_REVERSE_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_REVERSE_PARAM_J_IS_NZ(...) RADIX_2_REVERSE_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_2_REVERSE_MOTH_J_IS_NZ(...)  RADIX_2_REVERSE_MOTH_J_IS_NZ(__VA_ARGS__)

/************************* inverse butterfly **********************************
    4*a0 =            (b0 + b1) +        (b2 + b3)
    4*a1 =       w^-1*(b0 - b1) - i*w^-1*(b2 - b3)
    4*a2 = w^-2*(     (b0 + b1) -        (b2 + b3))
    4*a3 = w^-2*(w^-1*(b0 - b1) + i*w^-1*(b2 - b3))
    W  := -w^-1
    W2 := -w^-2
    IW := i*w^-1
*/
#define RADIX_4_REVERSE_PARAM_J_IS_Z(V, Q) \
    V IW = V##_set_d(Q->w2tab[0][1]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_MOTH_J_IS_Z(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    y0 = V##_add(x0, x1); \
    y1 = V##_add(x2, x3); \
    y2 = V##_sub(x0, x1); \
    y3 = V##_sub(x2, x3); \
    y0 = V##_reduce_to_pm1n(y0, n, ninv); \
    y1 = V##_reduce_to_pm1n(y1, n, ninv); \
    y2 = V##_reduce_to_pm1n(y2, n, ninv); \
    y3 = V##_mulmod(y3, IW, n, ninv); \
    V##_store(X0, V##_add(y0, y1)); \
    V##_store(X2, V##_sub(y0, y1)); \
    V##_store(X1, V##_sub(y2, y3)); \
    V##_store(X3, V##_add(y2, y3)); \
}

#define RADIX_4_REVERSE_PARAM_J_IS_NZ(V, Q, j_mr, j_bits) \
    V W  = V##_set_d(Q->w2tab[1+j_bits][2*j_mr+1]); \
    V W2 = V##_set_d(Q->w2tab[0+j_bits][j_mr]); \
    V IW = V##_set_d(Q->w2tab[1+j_bits][2*j_mr+0]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv); \

#define RADIX_4_REVERSE_MOTH_J_IS_NZ(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    y0 = V##_add(x0, x1); \
    y1 = V##_add(x2, x3); \
    y2 = V##_sub(x0, x1); \
    y3 = V##_sub(x3, x2); \
    y2 = V##_mulmod(y2, W, n, ninv); \
    y3 = V##_mulmod(y3, IW, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y3, y2); \
    V##_store(X1, x1); \
    x2 = V##_sub(y1, y0); \
    x3 = V##_add(y3, y2); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x2 = V##_mulmod(x2, W2, n, ninv); \
    x3 = V##_mulmod(x3, W2, n, ninv); \
    V##_store(X0, x0); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

#define _RADIX_4_REVERSE_PARAM_J_IS_Z(...)  RADIX_4_REVERSE_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_REVERSE_MOTH_J_IS_Z(...)   RADIX_4_REVERSE_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_REVERSE_PARAM_J_IS_NZ(...) RADIX_4_REVERSE_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_4_REVERSE_MOTH_J_IS_NZ(...)  RADIX_4_REVERSE_MOTH_J_IS_NZ(__VA_ARGS__)

/************ basecase inverse transform of size BLK_SZ **********************/
/*
    The basecases below 4 are disabled because the ifft expects input in
    slightly-worse-than-bit-reversed order as in basecase_4.
*/
#define DEFINE_IT(j_is_0) \
FLINT_FORCE_INLINE void sd_ifft_basecase_4_##j_is_0(\
    const sd_fft_lctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    vec4d n    = vec4d_set_d(Q->p); \
    vec4d ninv = vec4d_set_d(Q->pinv); \
    vec4d W, W2, IW, u, v; \
    vec4d x0, x1, x2, x3, y0, y1, y2, y3; \
 \
    x0 = vec4d_load(X + 0); \
    x1 = vec4d_load(X + 4); \
    x2 = vec4d_load(X + 8); \
    x3 = vec4d_load(X + 12); \
 \
    if (j_is_0) \
    { \
        W  = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[0][3], Q->w2tab[0][7], Q->w2tab[0][5]); \
        IW = vec4d_set_d4( Q->w2tab[0][1], Q->w2tab[0][2], Q->w2tab[0][6], Q->w2tab[0][4]); \
        W2 = vec4d_set_d4(-Q->w2tab[0][0], Q->w2tab[0][1], Q->w2tab[0][3], Q->w2tab[0][2]); \
    } \
    else \
    { \
        W2 = vec4d_load_aligned(Q->w2tab[2+j_bits] + 4*j_mr);   /* a b c d */ \
        W2 = vec4d_permute_3_2_1_0(W2);                         /* d c b a */ \
        u = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_mr+0);  /* 0 1 2 3 */ \
        v = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_mr+4);  /* 4 5 6 7 */ \
        W  = vec4d_unpackhi_permute_3_1_2_0(u, v);              /* 7 5 3 1 */ \
        IW = vec4d_unpacklo_permute_3_1_2_0(u, v);              /* 6 4 2 0 */ \
    } \
 \
    y0 = vec4d_add(x0, x1); \
    y1 = vec4d_add(x2, x3); \
    y2 = vec4d_sub(x0, x1); \
    y3 = vec4d_sub(x3, x2); \
    y2 = vec4d_mulmod(y2, W, n, ninv); \
    y3 = vec4d_mulmod(y3, IW, n, ninv); \
    x0 = vec4d_add(y0, y1); \
    x1 = vec4d_sub(y3, y2); \
    x2 = vec4d_sub(y1, y0); \
    x3 = vec4d_add(y3, y2); \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x2 = vec4d_mulmod(x2, W2, n, ninv); \
    x3 = vec4d_mulmod(x3, W2, n, ninv); \
 \
    VEC4D_TRANSPOSE(x0, x1, x2, x3, x0, x1, x2, x3); \
 \
    if (j_is_0) \
    { \
        IW = vec4d_set_d(Q->w2tab[0][1]); \
        y0 = vec4d_add(x0, x1); \
        y1 = vec4d_add(x2, x3); \
        y2 = vec4d_sub(x0, x1); \
        y3 = vec4d_sub(x2, x3); \
        y0 = vec4d_reduce_to_pm1n(y0, n, ninv); \
        y1 = vec4d_reduce_to_pm1n(y1, n, ninv); \
        y2 = vec4d_reduce_to_pm1n(y2, n, ninv); \
        y3 = vec4d_mulmod(y3, IW, n, ninv); \
        vec4d_store(X+0, vec4d_add(y0, y1)); \
        vec4d_store(X+8, vec4d_sub(y0, y1)); \
        vec4d_store(X+4, vec4d_sub(y2, y3)); \
        vec4d_store(X+12, vec4d_add(y2, y3)); \
    } \
    else \
    { \
        W  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+1]); \
        IW = vec4d_set_d(Q->w2tab[1+j_bits][2*j_mr+0]); \
        W2 = vec4d_set_d(Q->w2tab[0+j_bits][j_mr]); \
        y0 = vec4d_add(x0, x1); \
        y1 = vec4d_add(x2, x3); \
        y2 = vec4d_sub(x0, x1); \
        y3 = vec4d_sub(x3, x2); \
        y2 = vec4d_mulmod(y2, W, n, ninv); \
        y3 = vec4d_mulmod(y3, IW, n, ninv); \
        x0 = vec4d_add(y0, y1); \
        x1 = vec4d_sub(y3, y2); \
        x2 = vec4d_sub(y1, y0); \
        x3 = vec4d_add(y3, y2); \
        x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
        x2 = vec4d_mulmod(x2, W2, n, ninv); \
        x3 = vec4d_mulmod(x3, W2, n, ninv); \
        vec4d_store(X+0, x0); \
        vec4d_store(X+4, x1); \
        vec4d_store(X+8, x2); \
        vec4d_store(X+12, x3); \
    } \
}

DEFINE_IT(0)
DEFINE_IT(1)
#undef DEFINE_IT


/* use with n = m-2 and m >= 6 */
#define EXTEND_BASECASE(n, m) \
void CAT3(sd_ifft_basecase, m, 1)(const sd_fft_lctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    FLINT_ASSERT(j_bits == 0); \
    CAT3(sd_ifft_basecase, n, 1)(Q, X+0*l, 0, 0); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+1*l, 0, 1); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+2*l, 1, 2); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+3*l, 0, 2); \
    { \
        _RADIX_4_REVERSE_PARAM_J_IS_Z(VECND, Q) \
        ulong i = 0; do { \
            _RADIX_4_REVERSE_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += N, i < l); \
        FLINT_ASSERT(i == l); \
    } \
} \
void CAT3(sd_ifft_basecase, m, 0)(const sd_fft_lctx_t Q, double* X, ulong j_mr, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    FLINT_ASSERT(j_bits != 0); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+0*l, 4*j_mr+3, 2+j_bits); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+1*l, 4*j_mr+2, 2+j_bits); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+2*l, 4*j_mr+1, 2+j_bits); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+3*l, 4*j_mr+0, 2+j_bits); \
    { \
        _RADIX_4_REVERSE_PARAM_J_IS_NZ(VECND, Q, j_mr, j_bits) \
        ulong i = 0; do { \
            _RADIX_4_REVERSE_MOTH_J_IS_NZ(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
        } while (i += N, i < l); \
        FLINT_ASSERT(i == l); \
    } \
}

EXTEND_BASECASE(4, 6)
EXTEND_BASECASE(6, 8)
#undef EXTEND_BASECASE


/* use with n = m-1 and m >= 5 */
#define EXTEND_BASECASE_ODD(n, m)                                          \
void CAT3(sd_ifft_basecase, m, 1)(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 1); \
    CAT3(sd_ifft_basecase, n, 1)(Q, X+0*l, 0, 0); \
    CAT3(sd_ifft_basecase, n, 0)(Q, X+1*l, 0, 1); \
    _RADIX_2_REVERSE_PARAM_J_IS_Z(VECND, Q) \
    ulong i = 0; do { \
        _RADIX_2_REVERSE_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i) \
    } while (i += N, i < l); \
}
EXTEND_BASECASE_ODD(4, 5)
EXTEND_BASECASE_ODD(6, 7)
#undef EXTEND_BASECASE_ODD


#undef RADIX_2_REVERSE_PARAM_J_IS_Z
#undef RADIX_2_REVERSE_MOTH_J_IS_Z
#undef RADIX_2_REVERSE_PARAM_J_IS_NZ
#undef RADIX_2_REVERSE_MOTH_J_IS_NZ
#undef RADIX_4_REVERSE_PARAM_J_IS_Z
#undef RADIX_4_REVERSE_MOTH_J_IS_Z
#undef RADIX_4_REVERSE_PARAM_J_IS_NZ
#undef RADIX_4_REVERSE_MOTH_J_IS_NZ


#undef _RADIX_2_REVERSE_PARAM_J_IS_Z
#undef _RADIX_2_REVERSE_MOTH_J_IS_Z
#undef _RADIX_2_REVERSE_PARAM_J_IS_NZ
#undef _RADIX_2_REVERSE_MOTH_J_IS_NZ
#undef _RADIX_4_REVERSE_PARAM_J_IS_Z
#undef _RADIX_4_REVERSE_MOTH_J_IS_Z
#undef _RADIX_4_REVERSE_PARAM_J_IS_NZ
#undef _RADIX_4_REVERSE_MOTH_J_IS_NZ

#undef N
#undef VECND
#undef VECNOP
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* note: sd_ifft.c in the flint source tree does not undef    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                     new functions                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void sd_ifft_main(
    const sd_fft_lctx_t Q,
    ulong I, /* starting index */
    ulong S, /* stride */
    ulong k, /* transform length BLK_SZ*2^k */
    ulong j);

void do_sd_ifft_basecase_4_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits)
{
    sd_ifft_basecase_4_1(Q, X, j_r, j_bits);
}


void sd_ifft(sd_fft_lctx_t Q, double *X, ulong k)
{
    Q->data = X;
    sd_ifft_main(Q, 0, 1, k - LG_BLK_SZ, 0);
}
