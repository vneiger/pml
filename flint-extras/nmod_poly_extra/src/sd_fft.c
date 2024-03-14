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
    in the radix 4 butterflies.

        16x 4-wide AVX registers  => N = 8
        32x 2-wide NEON registers => N = 8
        32x 8-wide AVX512 registers  => N = ?
*/

#define N 8
#define VECND vec8d
#define VECNOP(op) vec8d_##op

/********************* forward butterfly **************************************
    b0 = a0 + w*a1
    b1 = a0 - w*a1
*/

#define RADIX_2_FORWARD_PARAM_J_IS_Z(V, Q) \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_FORWARD_MOTH_J_IS_Z(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x1 = V##_reduce_to_pm1n(x1, n, ninv); \
    V##_store(X0, V##_add(x0, x1)); \
    V##_store(X1, V##_sub(x0, x1)); \
}

#define RADIX_2_FORWARD_PARAM_J_IS_NZ(V, Q, j_r, j_bits) \
    V w = V##_set_d(Q->w2tab[j_bits][j_r]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_2_FORWARD_MOTH_J_IS_NZ(V, X0, X1) \
{ \
    V x0, x1; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x1 = V##_mulmod(x1, w, n, ninv); \
    V##_store(X0, V##_add(x0, x1)); \
    V##_store(X1, V##_sub(x0, x1)); \
}

/* for when the V arguments above needs "evaluation" */
#define _RADIX_2_FORWARD_PARAM_J_IS_Z(...)  RADIX_2_FORWARD_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_FORWARD_MOTH_J_IS_Z(...)   RADIX_2_FORWARD_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_2_FORWARD_PARAM_J_IS_NZ(...) RADIX_2_FORWARD_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_2_FORWARD_MOTH_J_IS_NZ(...)  RADIX_2_FORWARD_MOTH_J_IS_NZ(__VA_ARGS__)

/********************* forward butterfly **************************************
    b0 = a0 + w^2*a2 +   w*(a1 + w^2*a3)
    b1 = a0 + w^2*a2 -   w*(a1 + w^2*a3)
    b2 = a0 - w^2*a2 + i*w*(a1 - w^2*a3)
    b3 = a0 - w^2*a2 - i*w*(a1 - w^2*a3)
*/

#define RADIX_4_FORWARD_PARAM_J_IS_Z(V, Q) \
    V iw = V##_set_d(Q->w2tab[1][0]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_4_FORWARD_MOTH_J_IS_Z(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    x2 = V##_reduce_to_pm1n(x2, n, ninv); \
    x3 = V##_reduce_to_pm1n(x3, n, ninv); \
    y0 = V##_add(x0, x2); \
    y1 = V##_add(x1, x3); \
    y2 = V##_sub(x0, x2); \
    y3 = V##_sub(x1, x3); \
    y1 = V##_reduce_to_pm1n(y1, n, ninv); \
    y3 = V##_mulmod(y3, iw, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y0, y1); \
    x2 = V##_add(y2, y3); \
    x3 = V##_sub(y2, y3); \
    V##_store(X0, x0); \
    V##_store(X1, x1); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

#define RADIX_4_FORWARD_PARAM_J_IS_NZ(V, Q, j_r, j_bits) \
    FLINT_ASSERT(j_bits > 0); \
    V w  = V##_set_d(Q->w2tab[1+j_bits][2*j_r]); \
    V w2 = V##_set_d(Q->w2tab[0+j_bits][j_r]); \
    V iw = V##_set_d(Q->w2tab[1+j_bits][2*j_r+1]); \
    V n    = V##_set_d(Q->p); \
    V ninv = V##_set_d(Q->pinv);

#define RADIX_4_FORWARD_MOTH_J_IS_NZ(V, X0, X1, X2, X3) \
{ \
    V x0, x1, x2, x3, y0, y1, y2, y3; \
    x0 = V##_load(X0); \
    x0 = V##_reduce_to_pm1n(x0, n, ninv); \
    x1 = V##_load(X1); \
    x2 = V##_load(X2); \
    x3 = V##_load(X3); \
    x2 = V##_mulmod(x2, w2, n, ninv); \
    x3 = V##_mulmod(x3, w2, n, ninv); \
    y0 = V##_add(x0, x2); \
    y1 = V##_add(x1, x3); \
    y2 = V##_sub(x0, x2); \
    y3 = V##_sub(x1, x3); \
    y1 = V##_mulmod(y1, w, n, ninv); \
    y3 = V##_mulmod(y3, iw, n, ninv); \
    x0 = V##_add(y0, y1); \
    x1 = V##_sub(y0, y1); \
    x2 = V##_add(y2, y3); \
    x3 = V##_sub(y2, y3); \
    V##_store(X0, x0); \
    V##_store(X1, x1); \
    V##_store(X2, x2); \
    V##_store(X3, x3); \
}

#define _RADIX_4_FORWARD_PARAM_J_IS_Z(...)  RADIX_4_FORWARD_PARAM_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_FORWARD_MOTH_J_IS_Z(...)   RADIX_4_FORWARD_MOTH_J_IS_Z(__VA_ARGS__)
#define _RADIX_4_FORWARD_PARAM_J_IS_NZ(...) RADIX_4_FORWARD_PARAM_J_IS_NZ(__VA_ARGS__)
#define _RADIX_4_FORWARD_MOTH_J_IS_NZ(...)  RADIX_4_FORWARD_MOTH_J_IS_NZ(__VA_ARGS__)

/**************** basecase transform of size BLK_SZ **************************/
/*
    The basecases below 4 are disabled because the fft is expected to be
    produced in the slightly-worse-than-bit-reversed order of basecase_4.
*/
#define DEFINE_IT(j_is_0) \
FLINT_FORCE_INLINE void CAT(sd_fft_basecase_4, j_is_0)( \
    const sd_fft_lctx_t Q, \
    double* X, \
    ulong j_r, \
    ulong j_bits) \
{ \
    vec4d n    = vec4d_set_d(Q->p); \
    vec4d ninv = vec4d_set_d(Q->pinv); \
    vec4d w, w2, iw; \
    vec4d x0, x1, x2, x3, y0, y1, y2, y3, u, v; \
 \
    /* will abuse the fact that Q->w2tab[0] points to consecutive entries */ \
    FLINT_ASSERT(SD_FFT_CTX_INIT_DEPTH >= 4); \
 \
    x0 = vec4d_load(X+0); \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x1 = vec4d_load(X+4); \
    x2 = vec4d_load(X+8); \
    x3 = vec4d_load(X+12); \
 \
    if (j_is_0) \
    { \
        iw = vec4d_set_d(Q->w2tab[0][1]); \
 \
        x2 = vec4d_reduce_to_pm1n(x2, n, ninv); \
        x3 = vec4d_reduce_to_pm1n(x3, n, ninv); \
        y0 = vec4d_add(x0, x2); \
        y1 = vec4d_add(x1, x3); \
        y2 = vec4d_sub(x0, x2); \
        y3 = vec4d_sub(x1, x3); \
        y1 = vec4d_reduce_to_pm1n(y1, n, ninv); \
        y3 = vec4d_mulmod(y3, iw, n, ninv); \
        x0 = vec4d_add(y0, y1); \
        x1 = vec4d_sub(y0, y1); \
        x2 = vec4d_add(y2, y3); \
        x3 = vec4d_sub(y2, y3); \
    } \
    else \
    { \
        w  = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r]); \
        w2 = vec4d_set_d(Q->w2tab[0+j_bits][j_r]); \
        iw = vec4d_set_d(Q->w2tab[1+j_bits][2*j_r+1]); \
 \
        x2 = vec4d_mulmod(x2, w2, n, ninv); \
        x3 = vec4d_mulmod(x3, w2, n, ninv); \
        y0 = vec4d_add(x0, x2); \
        y1 = vec4d_add(x1, x3); \
        y2 = vec4d_sub(x0, x2); \
        y3 = vec4d_sub(x1, x3); \
        y1 = vec4d_mulmod(y1, w, n, ninv); \
        y3 = vec4d_mulmod(y3, iw, n, ninv); \
        x0 = vec4d_add(y0, y1); \
        x1 = vec4d_sub(y0, y1); \
        x2 = vec4d_add(y2, y3); \
        x3 = vec4d_sub(y2, y3); \
    } \
 \
    if (j_is_0) \
    { \
        u = vec4d_load_aligned(Q->w2tab[0] + 0); \
        v = vec4d_load_aligned(Q->w2tab[0] + 4); \
        w2 = u; \
    } \
    else \
    { \
        u  = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_r + 0); \
        v  = vec4d_load_aligned(Q->w2tab[3+j_bits] + 8*j_r + 4); \
        w2 = vec4d_load_aligned(Q->w2tab[2+j_bits] + 4*j_r + 0); \
    } \
    w  = vec4d_unpack_lo_permute_0_2_1_3(u, v); \
    iw = vec4d_unpack_hi_permute_0_2_1_3(u, v); \
 \
    VEC4D_TRANSPOSE(x0, x1, x2, x3, x0, x1, x2, x3); \
 \
    x0 = vec4d_reduce_to_pm1n(x0, n, ninv); \
    x2 = vec4d_mulmod(x2, w2, n, ninv); \
    x3 = vec4d_mulmod(x3, w2, n, ninv); \
    y0 = vec4d_add(x0, x2); \
    y1 = vec4d_add(x1, x3); \
    y2 = vec4d_sub(x0, x2); \
    y3 = vec4d_sub(x1, x3); \
    y1 = vec4d_mulmod(y1, w, n, ninv); \
    y3 = vec4d_mulmod(y3, iw, n, ninv); \
    x0 = vec4d_add(y0, y1); \
    x1 = vec4d_sub(y0, y1); \
    x2 = vec4d_add(y2, y3); \
    x3 = vec4d_sub(y2, y3); \
 \
    /* another VEC4D_TRANSPOSE here would put the output in bit-reversed */ \
    /* but this slow down is not necessary */ \
 \
    vec4d_store(X+0, x0); \
    vec4d_store(X+4, x1); \
    vec4d_store(X+8, x2); \
    vec4d_store(X+12, x3); \
}

DEFINE_IT(0)
DEFINE_IT(1)
#undef DEFINE_IT

/* use with n = m-2 and m >= 6 */
#define EXTEND_BASECASE(n, m) \
void CAT3(sd_fft_basecase, m, 1)(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    _RADIX_4_FORWARD_PARAM_J_IS_Z(VECND, Q) \
    ulong i = 0; do { \
        _RADIX_4_FORWARD_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 1)(Q, X+0*l, 0, 0); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 0, 1); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+2*l, 0, 2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+3*l, 1, 2); \
} \
void CAT3(sd_fft_basecase, m, 0)(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 2); \
    _RADIX_4_FORWARD_PARAM_J_IS_NZ(VECND, Q, j_r, j_bits) \
    ulong i = 0; do { \
        _RADIX_4_FORWARD_MOTH_J_IS_NZ(VECND, X+0*l+i, X+1*l+i, X+2*l+i, X+3*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+0*l, 4*j_r+0, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 4*j_r+1, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+2*l, 4*j_r+2, j_bits+2); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+3*l, 4*j_r+3, j_bits+2); \
}
EXTEND_BASECASE(4, 6)
EXTEND_BASECASE(6, 8)
#undef EXTEND_BASECASE


/* use with n = m-1 and m >= 5 */
#define EXTEND_BASECASE_ODD(n, m)                                          \
void CAT3(sd_fft_basecase, m, 1)(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits) \
{ \
    ulong l = n_pow2(m - 1); \
    _RADIX_2_FORWARD_PARAM_J_IS_Z(VECND, Q) \
    ulong i = 0; do { \
        _RADIX_2_FORWARD_MOTH_J_IS_Z(VECND, X+0*l+i, X+1*l+i) \
    } while (i += N, i < l); \
    CAT3(sd_fft_basecase, n, 1)(Q, X+0*l, 0, 0); \
    CAT3(sd_fft_basecase, n, 0)(Q, X+1*l, 0, 1); \
}
EXTEND_BASECASE_ODD(4, 5)
EXTEND_BASECASE_ODD(6, 7)
#undef EXTEND_BASECASE_ODD

#undef RADIX_2_FORWARD_PARAM_J_IS_Z
#undef RADIX_2_FORWARD_MOTH_J_IS_Z
#undef RADIX_2_FORWARD_PARAM_J_IS_NZ
#undef RADIX_2_FORWARD_MOTH_J_IS_NZ
#undef RADIX_4_FORWARD_PARAM_J_IS_Z
#undef RADIX_4_FORWARD_MOTH_J_IS_Z
#undef RADIX_4_FORWARD_PARAM_J_IS_NZ
#undef RADIX_4_FORWARD_MOTH_J_IS_NZ


#undef _RADIX_2_FORWARD_PARAM_J_IS_Z
#undef _RADIX_2_FORWARD_MOTH_J_IS_Z
#undef _RADIX_2_FORWARD_PARAM_J_IS_NZ
#undef _RADIX_2_FORWARD_MOTH_J_IS_NZ
#undef _RADIX_4_FORWARD_PARAM_J_IS_Z
#undef _RADIX_4_FORWARD_MOTH_J_IS_Z
#undef _RADIX_4_FORWARD_PARAM_J_IS_NZ
#undef _RADIX_4_FORWARD_MOTH_J_IS_NZ

#undef N
#undef VECND
#undef VECNOP

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                     new functions                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void sd_fft_main(
    const sd_fft_lctx_t Q,
    ulong I,    /* starting index */
    ulong S,    /* stride */
    ulong k,    /* 1 transform of length BLK_SZ*2^k */
    ulong j);    /* twist param */

void do_sd_fft_basecase_4_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits)
{
    sd_fft_basecase_4_1(Q, X, j_r, j_bits);
}


void sd_fft(sd_fft_lctx_t Q, double *X, ulong k)
{
    Q->data = X;
    sd_fft_main(Q, 0, 1, k - LG_BLK_SZ, 0);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
