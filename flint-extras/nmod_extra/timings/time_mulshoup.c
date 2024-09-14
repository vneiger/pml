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


#define NB_ITER 16384

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
            //array[j] = n_mulmod_shoup(a, array[j], a_pr, d);
            const ulong aj = array[j];

            ulong ajlo = __ll_lowpart(aj);
            ulong ajhi = __ll_highpart(aj);
            //ulong prod0 = ajlo * a_pr_lo;
            ulong prod1 = ajlo * a_pr_hi;
            ulong prod2 = ajhi * a_pr_lo;
            ulong p_hi = ajhi * a_pr_hi + ((prod1 + prod2) >> 32);
            ulong res = a * aj - p_hi * d;
            if (res >= 2*d)
                res -= 2*d;
            array[j] = res;
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

// __attribute__((optimize("no-tree-vectorize")))
#if defined(__AVX512DQ__) && defined(__AVX512VL__)
void sample_shoup_avx2(void * FLINT_UNUSED(arg), ulong count)
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
        const __m256i avec = _mm256_set1_epi64x(a);
        const __m256i dvec = _mm256_set1_epi64x(d);
        const __m256i a_pr_lo = _mm256_set1_epi64x(__ll_lowpart(a_pr));
        const __m256i a_pr_hi = _mm256_set1_epi64x(__ll_highpart(a_pr));
        const __m256i mask = _mm256_set1_epi64x(__ll_B - 1);
        for (ulong j = 0; j+3 < NB_ITER; j+=4)
        {
            const __m256i aj = _mm256_loadu_si256((__m256i *)(array+j));

            __m256i ajlo = _mm256_and_si256(aj, mask);
            __m256i ajhi = aj >> 32;
            // one of these two should be < 64 bits (well, is a_pr < 64 bits?)
            // if not, may need to first shift both right, then add?
            __m256i prod1 = _mm256_mul_epu32(ajlo, a_pr_hi);
            __m256i prod2 = _mm256_mul_epu32(ajhi, a_pr_lo);
            prod1 = _mm256_add_epi64(prod1, prod2);
            prod1 = _mm256_srli_epi64(prod1, 32);
            __m256i p_hi = _mm256_mul_epu32(ajhi, a_pr_hi);
            p_hi = _mm256_add_epi64(p_hi, prod1);

            // p_hi * d
            p_hi = _mm256_mullo_epi64(p_hi, dvec);
            __m256i prod = _mm256_mullo_epi64(avec, aj);

            __m256i res = _mm256_sub_epi64(prod, p_hi);
            //if (res >= 2*d)  // TODO
                //res -= 2*d;
            _mm256_storeu_si256((__m256i *)(array+j), res);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}
#endif


// __attribute__((optimize("no-tree-vectorize")))
void sample_shoup_autovec_2(void * FLINT_UNUSED(arg), ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr brray = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr brray_pr = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup

        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_randlimb(state);  // array[j] is arbitrary
            brray[j] = n_randint(state, d);  // array[j] is < d
            brray_pr[j] = n_mulmod_precomp_shoup(brray[j], d);
        }

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
        {
            const ulong aj = array[j];
            const ulong bj = brray[j];
            const ulong bj_pr = brray_pr[j];

            ulong ajlo = __ll_lowpart(aj);
            ulong ajhi = __ll_highpart(aj);
            ulong bj_pr_lo = __ll_lowpart(bj_pr);
            ulong bj_pr_hi = __ll_highpart(bj_pr);

            ulong prod1 = (ulong)ajlo * bj_pr_hi;
            ulong prod2 = (ulong)ajhi * bj_pr_lo;
            ulong p_hi = (ulong)ajhi * bj_pr_hi + ((prod1 + prod2) >> 32);

            ulong res = bj * aj - p_hi * d;
            /* if (res >= 2*d) */
            /*     res -= 2*d; */
            array[j] = res;
        }
        prof_stop();
    }

    flint_free(array);
    flint_free(brray);
    FLINT_TEST_CLEAR(state);
}

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
        for (ulong j = 0; j < NB_ITER; j++)
        {
            //array[j] = n_mulmod_shoup(a, array[j], a_pr, d);
            ulong p_hi, p_lo, res;
            ulong aj = array[j];
            umul_ppmm(p_hi, p_lo, a_pr, aj);
            res = a * aj - p_hi * d;
            //if (p_lo >= d)
            //    p_lo -= d;
            array[j] = res;
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;

    flint_printf("mulmod min time / max time is:\n");

    prof_repeat(&min, &max, sample, NULL);
    flint_printf("   - basic special red: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_avx2, NULL);
    flint_printf("   - basic special red AVX2: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_avx2_2, NULL);
    flint_printf("   - basic special red AVX2 (two arrays): %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

#if defined(__AVX512F__)
    prof_repeat(&min, &max, sample_avx512, NULL);
    flint_printf("   - basic special red AVX512: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);
#endif

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
#endif


    prof_repeat(&min, &max, sample_shoup_autovec_2, NULL);
    flint_printf("   - Shoup excluding precomputation, autovec (two arrays): %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    return 0;
}
