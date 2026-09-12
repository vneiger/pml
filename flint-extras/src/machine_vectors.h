/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef __MACHINE_VECTORS__H
#define __MACHINE_VECTORS__H

#include <flint/flint.h>

#include "pml.h"

#if PML_HAVE_MACHINE_VECTORS
# include <flint/machine_vectors.h>
#endif

#if PML_HAVE_AVX2
# if defined(__GNUC__)
#  include <immintrin.h>
# elif defined(_MSC_VER)
#  include <intrin.h>
# endif
#endif

/*-------------------------------------------------*/
/* some basic additions to FLINT's machine vectors */
/*-------------------------------------------------*/

#if PML_HAVE_MACHINE_VECTORS


#if PML_HAVE_AVX2  /* already in flint's machine_vectors for NEON/ARM64 */
/* returns a + b mod n, assuming a,b reduced mod n            */
FLINT_FORCE_INLINE vec1n vec1n_addmod(vec1n a, vec1n b, vec1n n)
{
    return n - b > a ? a + b : a + b - n;
}
#endif  /* PML_HAVE_AVX2 */

/* returns a + b mod n, assuming a,b reduced mod n            */
FLINT_FORCE_INLINE vec1d vec1d_addmod(vec1d a, vec1d b, vec1d n)
{
    return a + b - n >= 0 ? a + b - n : a + b;
}

#if PML_HAVE_AVX2
/* returns a + b mod n, assuming a,b reduced mod n            */
FLINT_FORCE_INLINE vec4d vec4d_addmod(vec4d a, vec4d b, vec4d n)
{
    return vec4d_reduce_2n_to_n(vec4d_add(a, b), n);
}
#endif  /* PML_HAVE_AVX2 */

/* loads a vec4n from a and converts it to double             */
#if PML_HAVE_AVX512
FLINT_FORCE_INLINE vec4d vec4d_load_unaligned_nn_ptr(nn_ptr a)
{
    return  _mm256_setr_m128d( _mm_cvtepi64_pd(_mm_loadu_si128((vec2n *) a)),
                               _mm_cvtepi64_pd(_mm_loadu_si128((vec2n *) (a + 2))) );
}
#else
FLINT_FORCE_INLINE vec4d vec4d_load_unaligned_nn_ptr(nn_ptr a)
{
    return vec4n_convert_limited_vec4d(vec4n_load_unaligned(a));
}
#endif

/* converts a vec4d to vec4n and stores it                    */
FLINT_FORCE_INLINE void vec4d_store_unaligned_nn_ptr(nn_ptr dest, vec4d a)
{
    vec4n_store_unaligned(dest, vec4d_convert_limited_vec4n(a));
}

#endif



/*----------------------------------------------------*/
/* some AVX-only additions to FLINT's machine vectors */
/*----------------------------------------------------*/

#if PML_HAVE_AVX2

FLINT_FORCE_INLINE void vec4n_store_aligned(ulong* z, vec4n a)
{
    _mm256_store_si256((__m256i*) z, a);
}

/* reduce_to_pm1n(a, n, ninv): return a mod n in (-n,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_to_pm1no(vec2d a, vec2d n, vec2d ninv)
{
    return _mm_fnmadd_pd(_mm_round_pd(_mm_mul_pd(a, ninv), 4), n, a);
}

/* reduce_pm1no_to_0n(a, n): return a mod n in [0,n) assuming a in (-n,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_pm1no_to_0n(vec2d a, vec2d n)
{
    return _mm_blendv_pd(a, _mm_add_pd(a, n), a);
}

/* reduce_to_0n(a, n, ninv): return a mod n in [0,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_to_0n(vec2d a, vec2d n, vec2d ninv)
{
    return vec2d_reduce_pm1no_to_0n(vec2d_reduce_to_pm1no(a, n, ninv), n);
}

FLINT_FORCE_INLINE vec2d vec2d_set_d2(double a1, double a0)
{
    return _mm_set_pd(a0, a1);
}

#define vec4n_bit_shift_right_45(a) vec4n_bit_shift_right((a), 45)


// avx2 horizontal sum  (already in FLINT)
FLINT_FORCE_INLINE ulong _mm256_hsum(__m256i a) {
    __m256i a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    __m256i sum_lo = _mm256_add_epi64(a, a_hi);
    __m128i sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (ulong) _mm_cvtsi128_si64(sum);
}

#endif /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
// avx512 horizontal sum
FLINT_FORCE_INLINE ulong _mm512_hsum(__m512i a) {
    return _mm512_reduce_add_epi64(a);
}
#endif /* PML_HAVE_AVX512 */

/*----------------------------------------------------*/
/* small transposes of ulong blocks                   */
/*----------------------------------------------------*/

/* TODO these ulong lane shuffles, and the 4x4 transpose built from them, are
 * the exact counterparts of vec4d_unpacklo, vec4d_unpackhi,
 * vec4d_permute2_0_2, vec4d_permute2_1_3 and VEC4D_TRANSPOSE, which
 * <flint/machine_vectors.h> already provides for both its AVX2 and its NEON
 * backend; the ulong lane type simply has no lane-combining primitive there
 * yet.  They belong next to their double counterparts in FLINT and should
 * move there; they are kept here in the meantime. */
#if PML_HAVE_MACHINE_VECTORS

#if PML_HAVE_AVX2

/* return {a[0], b[0], a[2], b[2]} */
FLINT_FORCE_INLINE vec4n vec4n_unpacklo(vec4n a, vec4n b)
{
    return _mm256_unpacklo_epi64(a, b);
}

/* return {a[1], b[1], a[3], b[3]} */
FLINT_FORCE_INLINE vec4n vec4n_unpackhi(vec4n a, vec4n b)
{
    return _mm256_unpackhi_epi64(a, b);
}

/* permute2_i0_i1(a, b): return {v[i0], v[i1]}
                     |   v[0]     |    v[1]    |    v[2]    |   v[3]     |
                       a[0], a[1]   a[2], a[3]   b[0], b[1]   b[2], b[3]  */
FLINT_FORCE_INLINE vec4n vec4n_permute2_0_2(vec4n a, vec4n b)
{
    return _mm256_permute2x128_si256(a, b, 0 + 16 * 2);
}

FLINT_FORCE_INLINE vec4n vec4n_permute2_1_3(vec4n a, vec4n b)
{
    return _mm256_permute2x128_si256(a, b, 1 + 16 * 3);
}

#else  /* NEON: vec4n is a pair of uint64x2_t, as vec4d is a pair of
          float64x2_t, and the same lane shuffles apply */

/* return {a[0], b[0]} */
FLINT_FORCE_INLINE vec2n vec2n_unpacklo(vec2n a, vec2n b)
{
    return vtrn1q_u64(a, b);
}

/* return {a[1], b[1]} */
FLINT_FORCE_INLINE vec2n vec2n_unpackhi(vec2n a, vec2n b)
{
    return vtrn2q_u64(a, b);
}

/* return {a[0], b[0], a[2], b[2]} */
FLINT_FORCE_INLINE vec4n vec4n_unpacklo(vec4n a, vec4n b)
{
    vec4n z = {vec2n_unpacklo(a.e1, b.e1), vec2n_unpacklo(a.e2, b.e2)};
    return z;
}

/* return {a[1], b[1], a[3], b[3]} */
FLINT_FORCE_INLINE vec4n vec4n_unpackhi(vec4n a, vec4n b)
{
    vec4n z = {vec2n_unpackhi(a.e1, b.e1), vec2n_unpackhi(a.e2, b.e2)};
    return z;
}

/* return {a[0], a[1], b[0], b[1]} */
FLINT_FORCE_INLINE vec4n vec4n_permute2_0_2(vec4n a, vec4n b)
{
    vec4n z = {a.e1, b.e1};
    return z;
}

/* return {a[2], a[3], b[2], b[3]} */
FLINT_FORCE_INLINE vec4n vec4n_permute2_1_3(vec4n a, vec4n b)
{
    vec4n z = {a.e2, b.e2};
    return z;
}

#endif  /* PML_HAVE_AVX2 */

/* view the 4 vectors as the rows of a 4x4 matrix */
#define VEC4N_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3)                  \
do {                                                                     \
    vec4n _s0 = vec4n_unpacklo(a0, a1);                                  \
    vec4n _s1 = vec4n_unpackhi(a0, a1);                                  \
    vec4n _s2 = vec4n_unpacklo(a2, a3);                                  \
    vec4n _s3 = vec4n_unpackhi(a2, a3);                                  \
    z0 = vec4n_permute2_0_2(_s0, _s2);                                   \
    z1 = vec4n_permute2_0_2(_s1, _s3);                                   \
    z2 = vec4n_permute2_1_3(_s0, _s2);                                   \
    z3 = vec4n_permute2_1_3(_s1, _s3);                                   \
} while (0)

#endif  /* PML_HAVE_MACHINE_VECTORS */

/* 8x8 transpose of ulong blocks, AVX-512: 24 shuffle uops for 64 words. */
#if PML_HAVE_AVX512

/* view the 8 vectors as the rows of an 8x8 matrix; named arguments rather
 * than arrays, as VEC4D_TRANSPOSE, so that nothing has to live in memory */
#define VEC8N_TRANSPOSE(z0, z1, z2, z3, z4, z5, z6, z7,                  \
                        a0, a1, a2, a3, a4, a5, a6, a7)                  \
do {                                                                     \
    __m512i _p0 = _mm512_unpacklo_epi64(a0, a1);                         \
    __m512i _p1 = _mm512_unpackhi_epi64(a0, a1);                         \
    __m512i _p2 = _mm512_unpacklo_epi64(a2, a3);                         \
    __m512i _p3 = _mm512_unpackhi_epi64(a2, a3);                         \
    __m512i _p4 = _mm512_unpacklo_epi64(a4, a5);                         \
    __m512i _p5 = _mm512_unpackhi_epi64(a4, a5);                         \
    __m512i _p6 = _mm512_unpacklo_epi64(a6, a7);                         \
    __m512i _p7 = _mm512_unpackhi_epi64(a6, a7);                         \
    __m512i _q0 = _mm512_shuffle_i64x2(_p0, _p2, 0x88);                  \
    __m512i _q1 = _mm512_shuffle_i64x2(_p1, _p3, 0x88);                  \
    __m512i _q2 = _mm512_shuffle_i64x2(_p0, _p2, 0xdd);                  \
    __m512i _q3 = _mm512_shuffle_i64x2(_p1, _p3, 0xdd);                  \
    __m512i _q4 = _mm512_shuffle_i64x2(_p4, _p6, 0x88);                  \
    __m512i _q5 = _mm512_shuffle_i64x2(_p5, _p7, 0x88);                  \
    __m512i _q6 = _mm512_shuffle_i64x2(_p4, _p6, 0xdd);                  \
    __m512i _q7 = _mm512_shuffle_i64x2(_p5, _p7, 0xdd);                  \
    z0 = _mm512_shuffle_i64x2(_q0, _q4, 0x88);                           \
    z1 = _mm512_shuffle_i64x2(_q1, _q5, 0x88);                           \
    z2 = _mm512_shuffle_i64x2(_q2, _q6, 0x88);                           \
    z3 = _mm512_shuffle_i64x2(_q3, _q7, 0x88);                           \
    z4 = _mm512_shuffle_i64x2(_q0, _q4, 0xdd);                           \
    z5 = _mm512_shuffle_i64x2(_q1, _q5, 0xdd);                           \
    z6 = _mm512_shuffle_i64x2(_q2, _q6, 0xdd);                           \
    z7 = _mm512_shuffle_i64x2(_q3, _q7, 0xdd);                           \
} while (0)

#endif  /* PML_HAVE_AVX512 */


#endif /* ifndef __MACHINE_VECTORS__H */
