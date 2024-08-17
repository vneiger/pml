#include <flint/ulong_extras.h>  // for umul_ppmm
#include "nmod_poly_integer_fft.h"

// returns a*b % n  in [0..2*n)
FLINT_FORCE_INLINE ulong n_mulmod_shoup_lazy(ulong a, ulong b, ulong apre, ulong n)
{
    ulong p_hi, p_lo;
    umul_ppmm(p_hi, p_lo, apre, b);
    return a * b - p_hi * n;
}


// reduction tree version
// lazy red: input in [0..2*n) --> output in [0..4*n)
FLINT_FORCE_INLINE void dft8_red_lazy(nn_ptr p, nmod_integer_fft_t F)
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

// DIF version
// lazy red: input in [0..2*n) --> output in [0..8*n)
FLINT_FORCE_INLINE void dft8_dif_lazy(nn_ptr p, nmod_integer_fft_t F)
{
    // in [0..2n) out [0..8n)
    ulong p0 = p[0];
    ulong q0 = p[1];
    ulong p1 = p[2];
    ulong q1 = p[3];
    ulong p2 = p[4];
    ulong q2 = p[5];
    ulong p3 = p[6];
    ulong q3 = p[7];

    ulong p4 = p0 + p2;  // [0..4n)
    ulong q4 = q0 + q2;  // [0..4n)
    ulong p6 = p1 + p3;  // [0..4n)
    ulong q6 = q1 + q3;  // [0..4n)
    ulong p5 = p0 + F->modn2 - p2;  // [0..4n)
    ulong q5 = q0 + F->modn2 - q2;  // [0..4n)
    if (p4 >= F->modn2)
        p4 -= F->modn2;  // [0..2n)
    if (p6 >= F->modn2)
        p6 -= F->modn2;  // [0..2n)
    if (q4 >= F->modn2)
        q4 -= F->modn2;  // [0..2n)
    if (q6 >= F->modn2)
        q6 -= F->modn2;  // [0..2n)

    ulong p_hi, p_lo;
    ulong p7 = (p1 + F->modn2 - p3);
    // {p,q}7 = ({p,q}1 + F->modn2 - {p,q}3) * F->tab_w[0][1]
    umul_ppmm(p_hi, p_lo, F->Ipre, p7);
    p7 = F->I * p7 - p_hi * F->mod.n;  // [0..2n)
    ulong q7 = q1 + F->modn2 - q3;
    umul_ppmm(p_hi, p_lo, F->Ipre, q7);
    q7 = F->I * q7 - p_hi * F->mod.n;  // [0..2n)

    p0 = p4 + p6;  // [0..4n)
    p1 = q4 + q6;  // [0..4n)
    p2 = p4 + F->modn2 - p6;  // [0..4n)

    p3 = q4 + F->modn2 - q6;
    umul_ppmm(p_hi, p_lo, F->Ipre, p3);
    p3 = F->I * p3 - p_hi * F->mod.n;  // [0..2n)

    q0 = p5 + p7;  // [0..6n)

    q1 = q5 + q7;  // [0..6n)
    umul_ppmm(p_hi, p_lo, F->Jpre, q1);
    q1 = F->J * q1 - p_hi * F->mod.n;  // [0..2n)

    q2 = p5 + F->modn2 - p7;  // [0..6n)
    q3 = q5 + F->modn2 - q7;  // [0..6n)
    umul_ppmm(p_hi, p_lo, F->Jpre, q3);
    q3 = F->IJ * q3 - p_hi * F->mod.n;  // [0..2n)

    p[0] = p0 + p1;  // [0..8n)
    p[1] = p0 + F->modn4 - p1;  // [0..8n)
    p[2] = p2 + p3;  // [0..6n)
    p[3] = p2 + F->modn2 - p3;  // [0..6n)
    p[4] = q0 + q1;  // [0..8n)
    p[5] = q0 + F->modn2 - q1;  // [0..8n)
    p[6] = q2 + q3;  // [0..8n)
    p[7] = q2 + F->modn2 - q3;  // [0..8n)
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
