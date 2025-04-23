/*
   Copyright 2024 (C) Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
   */

#include "flint/flint.h"
#include "flint/longlong.h"
#include "flint/profiler.h"
#include "flint/ulong_extras.h"
#include <immintrin.h>

#define __ll_30 (1L << 30)
#define __ll_lowpart30(t) ((ulong) (t) & (__ll_30 - 1))
#define __ll_highpart30(t) ((ulong) (t) >> 30)

/** a[i] >= n[i] ? a[i] - n[i] : a[i] */
FLINT_FORCE_INLINE __m256i _mm256_subtract_if_cmpge(__m256i a, __m256i n)
{
    return _mm256_min_epu64(a, _mm256_sub_epi64(a, n));
}

/** high word of widening 64x64 multiplication, lost carry */
FLINT_FORCE_INLINE __m256i _mm256_mulhi_lazy_epu64(__m256i a, __m256i b)
{
    __m256i ahi = _mm256_shuffle_epi32(a, 0xB1);
    __m256i bhi = _mm256_shuffle_epi32(b, 0xB1);

    __m256i alo_bhi = _mm256_mul_epu32(a, bhi);
    __m256i ahi_blo = _mm256_mul_epu32(ahi, b);

    __m256i mid = _mm256_add_epi64(_mm256_srli_epi64(alo_bhi, 32),
                                   _mm256_srli_epi64(ahi_blo, 32));

    __m256i ahi_bhi = _mm256_mul_epu32(ahi, bhi);
    return _mm256_add_epi64(ahi_bhi, mid);
}

/** Shoup's modular multiplication with precomputation, lazy */
FLINT_FORCE_INLINE __m256i
_mm256_mulmod_precomp_lazy(__m256i a, __m256i b, __m256i a_pr, __m256i n, __m256i n2)
{
    __m256i mulhi = _mm256_mulhi_lazy_epu64(b, a_pr);
    __m256i mullo1 = _mm256_mullo_epi64(b, a);
    __m256i mullo2 = _mm256_mullo_epi64(mulhi, n);
    __m256i mul = _mm256_sub_epi64(mullo1, mullo2);
    return _mm256_subtract_if_cmpge(mul, n2);
}

#define _MM256_SUBTRACT_IF_CMPGE(res, a, n) \
{ \
    res = _mm256_min_epu64(a, _mm256_sub_epi64(a, n)); \
}

#define _MM256_MULHI_LAZY_EPU64(res, a, b) \
{ \
    __m256i ahi = _mm256_shuffle_epi32(a, 0xB1); \
    __m256i bhi = _mm256_shuffle_epi32(b, 0xB1); \
 \
    __m256i alo_bhi = _mm256_mul_epu32(a, bhi); \
    __m256i ahi_blo = _mm256_mul_epu32(ahi, b); \
 \
    __m256i mid = _mm256_add_epi64(_mm256_srli_epi64(alo_bhi, 32), \
                                   _mm256_srli_epi64(ahi_blo, 32)); \
 \
    __m256i ahi_bhi = _mm256_mul_epu32(ahi, bhi); \
    res = _mm256_add_epi64(ahi_bhi, mid); \
}

#define _MM256_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n, n2) \
{ \
    __m256i mulhi; \
    _MM256_MULHI_LAZY_EPU64(mulhi, b, a_pr); \
    __m256i mullo1 = _mm256_mullo_epi64(b, a); \
    __m256i mullo2 = _mm256_mullo_epi64(mulhi, n); \
    __m256i mul = _mm256_sub_epi64(mullo1, mullo2); \
    _MM256_SUBTRACT_IF_CMPGE(res, mul, n2); \
}

#define _MM256_MULMOD_PRECOMP_LAZY_V2(res, a, b, a_pr, veca_pr_hi, n, n2) \
{ \
    __m256i v0, v1, v2; \
    v0 = _mm256_srli_epi64(b, 32);           \
    /* v3 = _mm256_shuffle_epi32(a_pr, 0xB1); */    \
    /* v0 = _mm256_shuffle_epi32(b,    0xB1); */    \
 \
    v1 = _mm256_mul_epu32(a_pr, v0); \
    v2 = _mm256_mul_epu32(veca_pr_hi, b); \
 \
    v1 = _mm256_add_epi64(_mm256_srli_epi64(v1, 32), \
                          _mm256_srli_epi64(v2, 32)); \
 \
    v2 = _mm256_mul_epu32(veca_pr_hi, v0); \
    v2 = _mm256_add_epi64(v2, v1); \
    \
    v0 = _mm256_mullo_epi64(b, a); \
    v1 = _mm256_mullo_epi64(v2, n); \
    v2 = _mm256_sub_epi64(v0, v1); \
    res = _mm256_min_epu64(v2, _mm256_sub_epi64(a, n2)); \
}

#define _MM512_MULMOD_PRECOMP_LAZY_V2(res, a, b, a_pr, veca_pr_hi, n, n2) \
{ \
    __m512i v0, v1, v2; \
    v0 = _mm512_srli_epi64(b, 32);           \
    /* v3 = _mm512_shuffle_epi32(a_pr, 0xB1); */    \
    /* v0 = _mm512_shuffle_epi32(b,    0xB1); */    \
 \
    v1 = _mm512_mul_epu32(a_pr, v0); \
    v2 = _mm512_mul_epu32(veca_pr_hi, b); \
 \
    v1 = _mm512_add_epi64(_mm512_srli_epi64(v1, 32), \
                          _mm512_srli_epi64(v2, 32)); \
 \
    v2 = _mm512_mul_epu32(veca_pr_hi, v0); \
    v2 = _mm512_add_epi64(v2, v1); \
    \
    v0 = _mm512_mullo_epi64(b, a); \
    v1 = _mm512_mullo_epi64(v2, n); \
    v2 = _mm512_sub_epi64(v0, v1); \
    res = _mm512_min_epu64(v2, _mm512_sub_epi64(a, n2)); \
}

#define NB_ITER 16384

// pretending modulus is special, like 2**60 - 2**30 - 1
void sample(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randlimb(state);  // a is arbitrary

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, d);  // must be < d

        prof_start();
        const ulong ahi = __ll_highpart(a);
        const ulong alo = __ll_lowpart(a);
        for (ulong j = 0; j < NB_ITER; j++)
        {
            //const ulong aj_pr = n_mulmod_precomp_shoup(array[j], d);
            //array[j] = n_mulmod_shoup(array[j], a, aj_pr, d);
            ulong ajhi = __ll_highpart30(array[j+0]);
            ulong ajlo = __ll_lowpart30(array[j+0]);
            ulong prod0 = ajlo * alo;
            ulong prod1 = ajlo * ahi;
            ulong prod2 = ajhi * alo;
            ulong prod3 = ajhi * ahi;
            prod1 = prod1 + prod2 + prod3;
            prod2 = prod1 >> 30;
            prod0 = prod0 + prod3 + prod2 + (__ll_lowpart30(prod1 + prod2) << 30);
            array[j] = prod0;
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

// pretending modulus is special, like 2**60 - 2**30 - 1
void sample_avx2(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randlimb(state);  // a is arbitrary

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, d);  // must be < d

        const __m256i ahi = _mm256_set1_epi64x(__ll_highpart(a));
        const __m256i alo = _mm256_set1_epi64x(__ll_lowpart(a));
        const __m256i mask = _mm256_set1_epi64x((1L << 30) - 1L);
        prof_start();
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            __m256i aj = _mm256_loadu_si256((__m256i*)(array+j));
            __m256i ajhi = _mm256_srli_epi64(aj, 30);
            __m256i ajlo = _mm256_and_si256(aj, mask);
            __m256i prod0 = _mm256_mul_epu32(ajlo, alo);
            __m256i prod1 = _mm256_mul_epu32(ajhi, alo);
            __m256i prod2 = _mm256_mul_epu32(ajlo, ahi);
            __m256i prod3 = _mm256_mul_epu32(ajhi, ahi);
            prod1 = _mm256_add_epi64(_mm256_add_epi64(prod1, prod2), prod3);
            prod2 = _mm256_srli_epi64(prod1, 30);
            prod0 = _mm256_add_epi64(_mm256_add_epi64(prod0, prod3), prod2);
            prod1 = _mm256_and_si256(_mm256_add_epi64(prod1, prod2), mask);
            prod0 = _mm256_add_epi64(prod0, _mm256_slli_epi64(prod1, 30));
            _mm256_storeu_si256((__m256i*)(array+j), prod0);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

// pretending modulus is special, like 2**60 - 2**30 - 1
void sample_avx2_2(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr brray = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup

        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_randint(state, d);  // must be < d
            brray[j] = n_randint(state, d);  // must be < d
        }

        const __m256i mask = _mm256_set1_epi64x((1L << 30) - 1L);
        prof_start();
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            __m256i aj = _mm256_loadu_si256((__m256i*)(array+j));
            __m256i bj = _mm256_loadu_si256((__m256i*)(brray+j));
            __m256i ajhi = _mm256_srli_epi64(aj, 30);
            __m256i ajlo = _mm256_and_si256(aj, mask);
            __m256i ahi = _mm256_srli_epi64(bj, 30);
            __m256i alo = _mm256_and_si256(bj, mask);
            __m256i prod0 = _mm256_mul_epu32(ajlo, alo);
            __m256i prod1 = _mm256_mul_epu32(ajhi, alo);
            __m256i prod2 = _mm256_mul_epu32(ajlo, ahi);
            __m256i prod3 = _mm256_mul_epu32(ajhi, ahi);
            prod1 = _mm256_add_epi64(_mm256_add_epi64(prod1, prod2), prod3);
            prod2 = _mm256_srli_epi64(prod1, 30);
            prod0 = _mm256_add_epi64(_mm256_add_epi64(prod0, prod3), prod2);
            prod1 = _mm256_and_si256(_mm256_add_epi64(prod1, prod2), mask);
            prod0 = _mm256_add_epi64(prod0, _mm256_slli_epi64(prod1, 30));
            _mm256_storeu_si256((__m256i*)(array+j), prod0);
        }
        prof_stop();
    }

    flint_free(array);
    flint_free(brray);
    FLINT_TEST_CLEAR(state);
}

// pretending modulus is special, like 2**60 - 2**30 - 1
#if defined(__AVX512F__)
void sample_avx512(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = flint_malloc(NB_ITER * sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randlimb(state);  // a is arbitrary

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, d);  // must be < d

        const __m512i ahi = _mm512_set1_epi64(__ll_highpart(a));
        const __m512i alo = _mm512_set1_epi64(__ll_lowpart(a));
        const __m512i mask = _mm512_set1_epi64((1L << 30) - 1L);
        prof_start();
        for (ulong j = 0; j+7 < NB_ITER; j+=8)
        {
            __m512i aj = _mm512_loadu_si512((__m512i*)(array+j));
            __m512i ajhi = _mm512_srli_epi64(aj, 30);
            __m512i ajlo = _mm512_and_si512(aj, mask);
            __m512i prod0 = _mm512_mul_epu32(ajlo, alo);
            __m512i prod1 = _mm512_mul_epu32(ajhi, alo);
            __m512i prod2 = _mm512_mul_epu32(ajlo, ahi);
            __m512i prod3 = _mm512_mul_epu32(ajhi, ahi);
            prod1 = _mm512_add_epi64(_mm512_add_epi64(prod1, prod2), prod3);
            prod2 = _mm512_srli_epi64(prod1, 30);
            prod0 = _mm512_add_epi64(_mm512_add_epi64(prod0, prod3), prod2);
            prod1 = _mm512_and_si512(_mm512_add_epi64(prod1, prod2), mask);
            prod0 = _mm512_add_epi64(prod0, _mm512_slli_epi64(prod1, 30));
            _mm512_storeu_si512((__m512i*)(array+j), prod0);
        }
        prof_stop();
    }

    flint_aligned_free(array);
    FLINT_TEST_CLEAR(state);
}
#endif






void sample_shoup(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            // unrolling helps: ~2.0c -> ~1.5c
            ulong p_hi, p_lo, res0, res1, res2, res3;
            ulong aj0 = array[j+0];
            ulong aj1 = array[j+1];
            ulong aj2 = array[j+2];
            ulong aj3 = array[j+3];
            umul_ppmm(p_hi, p_lo, a_pr, aj0);
            res0 = a * aj0 - p_hi * d;
            umul_ppmm(p_hi, p_lo, a_pr, aj1);
            res1 = a * aj1 - p_hi * d;
            umul_ppmm(p_hi, p_lo, a_pr, aj2);
            res2 = a * aj2 - p_hi * d;
            umul_ppmm(p_hi, p_lo, a_pr, aj3);
            res3 = a * aj3 - p_hi * d;
            array[j+0] = res0;
            array[j+1] = res1;
            array[j+2] = res2;
            array[j+3] = res3;
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

// __attribute__((optimize("no-tree-vectorize")))
void sample_shoup_autovec(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);
        const ulong a_pr_lo = __ll_lowpart(a_pr);
        const ulong a_pr_hi = __ll_highpart(a_pr);

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
        {
            ulong aj = array[j];

            ulong ajlo = __ll_lowpart(aj);
            ulong ajhi = __ll_highpart(aj);
            ulong prod1 = ajlo * a_pr_hi;
            ulong prod2 = ajhi * a_pr_lo;
            ulong p_hi = ajhi * a_pr_hi + ((prod1>>32) + (prod2>>32));
            aj = a * aj - p_hi * d;
            if (aj >= 2*d)
                aj -= 2*d;
            array[j] = aj;
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

#if defined(__AVX512DQ__) && defined(__AVX512VL__)

void sample_shoup_avx2_function(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 2) + 1;  // 1...62
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        const __m256i veca = _mm256_set1_epi64x(a);
        const __m256i vecn = _mm256_set1_epi64x(d);
        const __m256i vecn2 = _mm256_set1_epi64x(d);
        const __m256i veca_pr = _mm256_set1_epi64x(a_pr);
        //const __m256i a_pr_lo = _mm256_set1_epi64x(__ll_lowpart(a_pr));
        //const __m256i a_pr_hi = _mm256_set1_epi64x(__ll_highpart(a_pr));
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            const __m256i vecaj = _mm256_loadu_si256((__m256i *)(array+j));

            __m256i res = _mm256_mulmod_precomp_lazy(veca, vecaj, veca_pr, vecn, vecn2);
            _mm256_storeu_si256((__m256i *)(array+j), res);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_shoup_avx2_macro(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 2) + 1;  // 1...62
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        const __m256i veca = _mm256_set1_epi64x(a);
        const __m256i vecn = _mm256_set1_epi64x(d);
        const __m256i vecn2 = _mm256_set1_epi64x(d);
        const __m256i veca_pr = _mm256_set1_epi64x(a_pr);
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            __m256i vecaj = _mm256_loadu_si256((__m256i *)(array+j));
            //_MM256_MULMOD_PRECOMP_LAZY_V2(vecaj, veca, vecaj, veca_pr, veca_pr_hi, vecn, vecn2);
            const __m256i veca_pr_hi = _mm256_set1_epi64x(__ll_highpart(a_pr));
            _MM256_MULMOD_PRECOMP_LAZY_V2(vecaj, veca, vecaj, veca_pr, veca_pr_hi, vecn, vecn2);
            _mm256_storeu_si256((__m256i *)(array+j), vecaj);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_shoup_avx2(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 2) + 1;  // 1...62
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        const __m256i veca = _mm256_set1_epi64x(a);
        const __m256i vecn = _mm256_set1_epi64x(d);
        const __m256i vecn2 = _mm256_set1_epi64x(2*d);
        const __m256i veca_pr = _mm256_set1_epi64x(a_pr);
        const __m256i a_pr_hi = _mm256_set1_epi64x(__ll_highpart(a_pr));
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            __m256i aj = _mm256_loadu_si256((__m256i *)(array+j));

            __m256i p_hi = aj >> 32;
            __m256i prod1 = _mm256_mul_epu32(aj, a_pr_hi);
            aj = _mm256_mullo_epi64(veca, aj);

            __m256i prod2 = _mm256_mul_epu32(p_hi, veca_pr);
            prod1 = (prod1>>32) + (prod2>>32);
            p_hi = _mm256_mul_epu32(p_hi, a_pr_hi);
            p_hi = p_hi + prod1;

            p_hi = _mm256_mullo_epi64(p_hi, vecn);

            aj = _mm256_sub_epi64(aj, p_hi);
            aj = _mm256_subtract_if_cmpge(aj, vecn2);

            _mm256_storeu_si256((__m256i *)(array+j), aj);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_shoup2_avx2_macro(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 2) + 1;  // 1...62
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        const __m256i veca = _mm256_set1_epi64x(a);
        const __m256i vecn = _mm256_set1_epi64x(d);
        const __m256i vecn2 = _mm256_set1_epi64x(d);
        const __m256i veca_pr = _mm256_set1_epi64x(a_pr);
        for (ulong j = 0; j+7 < NB_ITER; j+=8)
        {
            __m256i vecaj = _mm256_loadu_si256((__m256i *)(array+j));
            __m256i vecaj2 = _mm256_loadu_si256((__m256i *)(array+j+4));
            const __m256i veca_pr_hi = _mm256_srli_epi64(veca_pr, 32);
            _MM256_MULMOD_PRECOMP_LAZY_V2(vecaj, veca, vecaj, veca_pr, veca_pr_hi, vecn, vecn2);
            _MM256_MULMOD_PRECOMP_LAZY_V2(vecaj2, veca, vecaj2, veca_pr, veca_pr_hi, vecn, vecn2);
            _mm256_storeu_si256((__m256i *)(array+j), vecaj);
            _mm256_storeu_si256((__m256i *)(array+j+4), vecaj2);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_shoup_avx512_macro(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 2) + 1;  // 1...62
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        const __m512i veca = _mm512_set1_epi64(a);
        const __m512i vecn = _mm512_set1_epi64(d);
        const __m512i vecn2 = _mm512_set1_epi64(d);
        const __m512i veca_pr = _mm512_set1_epi64(a_pr);
        for (ulong j = 0; j+7 < NB_ITER; j+=8)
        {
            __m512i vecaj = _mm512_loadu_si512((__m512i *)(array+j));
            const __m512i veca_pr_hi = _mm512_srli_epi64(veca_pr, 32);
            _MM512_MULMOD_PRECOMP_LAZY_V2(vecaj, veca, vecaj, veca_pr, veca_pr_hi, vecn, vecn2);
            _mm512_storeu_si512((__m512i *)(array+j), vecaj);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}
#endif

int main(void)
{
    double min, max;

    flint_printf("mulmod min time / max time is:\n");

//    prof_repeat(&min, &max, sample, NULL);
//    flint_printf("   - basic special red: %.3f cycles / %.3f cycles\n",
//            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);
//
//    prof_repeat(&min, &max, sample_avx2, NULL);
//    flint_printf("   - basic special red AVX2: %.3f cycles / %.3f cycles\n",
//            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);
//
//    prof_repeat(&min, &max, sample_avx2_2, NULL);
//    flint_printf("   - basic special red AVX2 (two arrays): %.3f cycles / %.3f cycles\n",
//            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);
//
//#if defined(__AVX512F__)
//    prof_repeat(&min, &max, sample_avx512, NULL);
//    flint_printf("   - basic special red AVX512: %.3f cycles / %.3f cycles\n",
//            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);
//#endif

    prof_repeat(&min, &max, sample_shoup, NULL);
    flint_printf("   - Shoup excluding precomputation: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_shoup_autovec, NULL);
    flint_printf("   - Shoup excluding precomputation, autovec: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

#if defined(__AVX512DQ__) && defined(__AVX512VL__)
    prof_repeat(&min, &max, sample_shoup_avx2, NULL);
    flint_printf("   - Shoup excluding precomputation, avx2: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_shoup_avx2_function, NULL);
    flint_printf("   - Shoup excluding precomputation lazy, avx2: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_shoup_avx2_macro, NULL);
    flint_printf("   - Shoup excluding precomputation macro, avx2: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_shoup2_avx2_macro, NULL);
    flint_printf("   - Shoup excluding precomputation macro, avx2: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_shoup_avx512_macro, NULL);
    flint_printf("   - Shoup excluding precomputation macro, avx512: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);
#endif

    return 0;
}
