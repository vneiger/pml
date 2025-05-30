#include <immintrin.h>

#include "flint/nmod_vec.h"  // for DOT_SPLIT_MASK
#include "nmod32_vec.h"

// horizontal sum
FLINT_FORCE_INLINE ulong _mm256_hsum(__m256i a) {
    __m256i a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    __m256i sum_lo = _mm256_add_epi64(a, a_hi);
    __m128i sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (ulong) _mm_cvtsi128_si64(sum);
}

#if HAVE_AVX512   // TODO handle AVX flags
FLINT_FORCE_INLINE ulong _mm512_hsum(__m512i a) {
    return _mm512_reduce_add_epi64(a);
}
#endif

void _nmod32_vec_dot2_split(uint * res0, uint * res1, n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, slong len, nmod_t mod, ulong pow2_precomp)
{
    ulong dp_lo0 = 0;
    ulong dp_lo1 = 0;
    ulong dp_hi0 = 0;
    ulong dp_hi1 = 0;

    slong i = 0;

    for ( ; i+3 < len; i+=4)
    {
        dp_lo0 += (ulong)vec1[i+0] * (ulong)vec2_0[i+0];
        dp_lo1 += (ulong)vec1[i+0] * (ulong)vec2_1[i+0];

        dp_lo0 += (ulong)vec1[i+1] * (ulong)vec2_0[i+1];
        dp_lo1 += (ulong)vec1[i+1] * (ulong)vec2_1[i+1];

        dp_lo0 += (ulong)vec1[i+2] * (ulong)vec2_0[i+2];
        dp_lo1 += (ulong)vec1[i+2] * (ulong)vec2_1[i+2];

        dp_lo0 += (ulong)vec1[i+3] * (ulong)vec2_0[i+3];
        dp_lo1 += (ulong)vec1[i+3] * (ulong)vec2_1[i+3];

        dp_hi0 += (dp_lo0 >> DOT_SPLIT_BITS);
        dp_hi1 += (dp_lo1 >> DOT_SPLIT_BITS);
        dp_lo0 &= DOT_SPLIT_MASK;
        dp_lo1 &= DOT_SPLIT_MASK;
    }

    // less than 4 terms remaining, can be accumulated
    for ( ; i < len; i++)
    {
        dp_lo0 += (ulong)vec1[i] * (ulong)vec2_0[i];
        dp_lo1 += (ulong)vec1[i] * (ulong)vec2_1[i];
    }
    dp_hi0 += (dp_lo0 >> DOT_SPLIT_BITS);
    dp_hi1 += (dp_lo1 >> DOT_SPLIT_BITS);
    dp_lo0 &= DOT_SPLIT_MASK;
    dp_lo1 &= DOT_SPLIT_MASK;

    NMOD_RED(*res0, pow2_precomp * dp_hi0 + dp_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * dp_hi1 + dp_lo1, mod);
}

void _nmod32_vec_dot2_split_avx2(uint * res0, uint * res1,
                                 n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                 slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_lo1 = _mm256_setzero_si256();
    __m256i dp_hi0 = _mm256_setzero_si256();
    __m256i dp_hi1 = _mm256_setzero_si256();

    slong i = 0;

    for ( ; i+15 < len; i += 16)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i v1_1 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 8));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        __m256i v2_0_1 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 8));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        __m256i v2_1_1 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 8));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        __m256i dp_lo2 = _mm256_mul_epu32(v1_1, v2_0_1);
        __m256i dp_lo3 = _mm256_mul_epu32(v1_1, v2_1_1);

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm256_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm256_shuffle_epi32(v2_1_1, 0xB1);
        // alternative 2: vpsrlq
        //v1_0 = _mm256_srli_epi64(v1_0, 32);
        //v2_0 = _mm256_srli_epi64(v2_0, 32);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo2);
        dp_lo1 = _mm256_add_epi64(dp_lo1, dp_lo3);
        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
        dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);
    }

    for ( ; i+7 < len; i += 8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        // alternative 2: vpsrlq (tested)
        //v1_0 = _mm256_srli_epi64(v1_0, 32);
        //v2_0 = _mm256_srli_epi64(v2_0, 32);
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));

        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
        dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);
    }

    ulong hsum_lo0 = _mm256_hsum(dp_lo0);
    ulong hsum_hi0 = _mm256_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm256_hsum(dp_lo1);
    ulong hsum_hi1 = _mm256_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    // some accumulation could be done here (but not up to 7 terms)
    for (; i < len; i++)
    {
        hsum_lo0 += (ulong)vec1[i] * (ulong)vec2_0[i];
        hsum_hi0 += (hsum_lo0 >> DOT_SPLIT_BITS);
        hsum_lo0 &= DOT_SPLIT_MASK;

        hsum_lo1 += (ulong)vec1[i] * (ulong)vec2_1[i];
        hsum_hi1 += (hsum_lo1 >> DOT_SPLIT_BITS);
        hsum_lo1 &= DOT_SPLIT_MASK;
    }

    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}

#if HAVE_AVX512

void _nmod32_vec_dot2_split_avx512(uint * res0, uint * res1,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1,
                                   slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);
    __m512i dp_lo0 = _mm512_setzero_si512();
    __m512i dp_lo1 = _mm512_setzero_si512();
    __m512i dp_hi0 = _mm512_setzero_si512();
    __m512i dp_hi1 = _mm512_setzero_si512();

    slong i = 0;

    for ( ; i+31 < len; i += 32)
    {
        __m512i v1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i+ 0));
        __m512i v1_1 = _mm512_loadu_si512((const __m512i *) (vec1+i+16));
        __m512i v2_0_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+ 0));
        __m512i v2_0_1 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+16));
        __m512i v2_1_0 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+ 0));
        __m512i v2_1_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+16));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        __m512i dp_lo2 = _mm512_mul_epu32(v1_1, v2_0_1);
        __m512i dp_lo3 = _mm512_mul_epu32(v1_1, v2_1_1);

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm512_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm512_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm512_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm512_shuffle_epi32(v2_1_1, 0xB1);
        // alternative 2: vpsrlq
        //v1_0 = _mm512_srli_epi64(v1_0, 32);
        //v2_0 = _mm512_srli_epi64(v2_0, 32);
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo2);
        dp_lo1 = _mm512_add_epi64(dp_lo1, dp_lo3);
        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
        dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);
    }

    for ( ; i+15 < len; i += 16)
    {
        __m512i v1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i));
        __m512i v2_0_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i));
        __m512i v2_1_0 = _mm512_loadu_si512((const __m512i *) (vec2_1+i));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v2_0_0 = _mm512_shuffle_epi32(v2_0_0, 0xB1);
        v2_1_0 = _mm512_shuffle_epi32(v2_1_0, 0xB1);
        // alternative 2: vpsrlq (tested)
        //v1_0 = _mm512_srli_epi64(v1_0, 32);
        //v2_0 = _mm512_srli_epi64(v2_0, 32);
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));

        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
        dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);
    }

    ulong hsum_lo0 = _mm512_hsum(dp_lo0);
    ulong hsum_hi0 = _mm512_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm512_hsum(dp_lo1);
    ulong hsum_hi1 = _mm512_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    // some accumulation could be done here (but not up to 7 terms)
    for (; i < len; i++)
    {
        hsum_lo0 += (ulong)vec1[i] * (ulong)vec2_0[i];
        hsum_hi0 += (hsum_lo0 >> DOT_SPLIT_BITS);
        hsum_lo0 &= DOT_SPLIT_MASK;

        hsum_lo1 += (ulong)vec1[i] * (ulong)vec2_1[i];
        hsum_hi1 += (hsum_lo1 >> DOT_SPLIT_BITS);
        hsum_lo1 &= DOT_SPLIT_MASK;
    }

    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}

#endif

// dot3_split_avx2, tested ok, gain not big
void _nmod32_vec_dot3_split_avx2(uint * res0, uint * res1, uint * res2,
                                 n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                 slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_lo1 = _mm256_setzero_si256();
    __m256i dp_lo2 = _mm256_setzero_si256();
    __m256i dp_hi0 = _mm256_setzero_si256();
    __m256i dp_hi1 = _mm256_setzero_si256();
    __m256i dp_hi2 = _mm256_setzero_si256();

    slong i = 0;

    for ( ; i+15 < len; i += 16)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i v1_1 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 8));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        __m256i v2_0_1 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 8));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        __m256i v2_1_1 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 8));
        __m256i v2_2_0 = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 0));
        __m256i v2_2_1 = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 8));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_0, v2_2_0));
        __m256i dp_lo3 = _mm256_mul_epu32(v1_1, v2_0_1);
        __m256i dp_lo4 = _mm256_mul_epu32(v1_1, v2_1_1);
        __m256i dp_lo5 = _mm256_mul_epu32(v1_1, v2_2_1);

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm256_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm256_shuffle_epi32(v2_1_1, 0xB1);
        v2_2_0 = _mm256_shuffle_epi32(v2_2_0, 0xB1);
        v2_2_1 = _mm256_shuffle_epi32(v2_2_1, 0xB1);
        // alternative 2: vpsrlq
        //v1_0 = _mm256_srli_epi64(v1_0, 32);
        //v2_0 = _mm256_srli_epi64(v2_0, 32);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_0, v2_2_0));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo4 = _mm256_add_epi64(dp_lo4, _mm256_mul_epu32(v1_1, v2_1_1));
        dp_lo5 = _mm256_add_epi64(dp_lo5, _mm256_mul_epu32(v1_1, v2_2_1));

        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo3);
        dp_lo1 = _mm256_add_epi64(dp_lo1, dp_lo4);
        dp_lo2 = _mm256_add_epi64(dp_lo2, dp_lo5);
        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_hi2 = _mm256_add_epi64(dp_hi2, _mm256_srli_epi64(dp_lo2, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
        dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);
        dp_lo2 = _mm256_and_si256(dp_lo2, low_bits);
    }

    for ( ; i+7 < len; i += 8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i));
        __m256i v2_2_0 = _mm256_loadu_si256((const __m256i *) (vec2_2+i));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_0, v2_2_0));

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_2_0 = _mm256_shuffle_epi32(v2_2_0, 0xB1);
        // alternative 2: vpsrlq (tested)
        //v1_0 = _mm256_srli_epi64(v1_0, 32);
        //v2_0 = _mm256_srli_epi64(v2_0, 32);
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_0, v2_2_0));

        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
        dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);
        dp_hi2 = _mm256_add_epi64(dp_hi2, _mm256_srli_epi64(dp_lo2, DOT_SPLIT_BITS));
        dp_lo2 = _mm256_and_si256(dp_lo2, low_bits);
    }

    ulong hsum_lo0 = _mm256_hsum(dp_lo0);
    ulong hsum_hi0 = _mm256_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm256_hsum(dp_lo1);
    ulong hsum_hi1 = _mm256_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;
    ulong hsum_lo2 = _mm256_hsum(dp_lo2);
    ulong hsum_hi2 = _mm256_hsum(dp_hi2) + (hsum_lo2 >> DOT_SPLIT_BITS);
    hsum_lo2 &= DOT_SPLIT_MASK;

    // some accumulation could be done here (but not up to 7 terms)
    for (; i < len; i++)
    {
        hsum_lo0 += (ulong)vec1[i] * (ulong)vec2_0[i];
        hsum_hi0 += (hsum_lo0 >> DOT_SPLIT_BITS);
        hsum_lo0 &= DOT_SPLIT_MASK;

        hsum_lo1 += (ulong)vec1[i] * (ulong)vec2_1[i];
        hsum_hi1 += (hsum_lo1 >> DOT_SPLIT_BITS);
        hsum_lo1 &= DOT_SPLIT_MASK;

        hsum_lo2 += (ulong)vec1[i] * (ulong)vec2_2[i];
        hsum_hi2 += (hsum_lo2 >> DOT_SPLIT_BITS);
        hsum_lo2 &= DOT_SPLIT_MASK;
    }

    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
    NMOD_RED(*res2, pow2_precomp * hsum_hi2 + hsum_lo2, mod);
}

#if HAVE_AVX512

void _nmod32_vec_dot4_split_avx512(uint * res0, uint * res1, uint * res2, uint * res3,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2, n32_srcptr vec2_3,
                                   slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);
    __m512i dp_lo0 = _mm512_setzero_si512();
    __m512i dp_lo1 = _mm512_setzero_si512();
    __m512i dp_lo2 = _mm512_setzero_si512();
    __m512i dp_lo3 = _mm512_setzero_si512();
    __m512i dp_hi0 = _mm512_setzero_si512();
    __m512i dp_hi1 = _mm512_setzero_si512();
    __m512i dp_hi2 = _mm512_setzero_si512();
    __m512i dp_hi3 = _mm512_setzero_si512();

    slong i = 0;

    for ( ; i+15 < len; i += 16)
    {
        __m512i v1 = _mm512_loadu_si512((const __m512i *) (vec1+i));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i));
        __m512i v2_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i));
        __m512i v2_2 = _mm512_loadu_si512((const __m512i *) (vec2_2+i));
        __m512i v2_3 = _mm512_loadu_si512((const __m512i *) (vec2_3+i));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1, v2_2));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1, v2_3));

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1 = _mm512_shuffle_epi32(v1, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);
        v2_2 = _mm512_shuffle_epi32(v2_2, 0xB1);
        v2_3 = _mm512_shuffle_epi32(v2_3, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1, v2_2));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1, v2_3));

        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
        dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);
        dp_hi2 = _mm512_add_epi64(dp_hi2, _mm512_srli_epi64(dp_lo2, DOT_SPLIT_BITS));
        dp_lo2 = _mm512_and_si512(dp_lo2, low_bits);
        dp_hi3 = _mm512_add_epi64(dp_hi3, _mm512_srli_epi64(dp_lo3, DOT_SPLIT_BITS));
        dp_lo3 = _mm512_and_si512(dp_lo3, low_bits);
    }

    ulong hsum_lo0 = _mm512_hsum(dp_lo0);
    ulong hsum_hi0 = _mm512_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm512_hsum(dp_lo1);
    ulong hsum_hi1 = _mm512_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;
    ulong hsum_lo2 = _mm512_hsum(dp_lo2);
    ulong hsum_hi2 = _mm512_hsum(dp_hi2) + (hsum_lo2 >> DOT_SPLIT_BITS);
    hsum_lo2 &= DOT_SPLIT_MASK;
    ulong hsum_lo3 = _mm512_hsum(dp_lo3);
    ulong hsum_hi3 = _mm512_hsum(dp_hi3) + (hsum_lo3 >> DOT_SPLIT_BITS);
    hsum_lo3 &= DOT_SPLIT_MASK;

    // some accumulation could be done here (but not up to 7 terms)
    for (; i < len; i++)
    {
        hsum_lo0 += (ulong)vec1[i] * (ulong)vec2_0[i];
        hsum_hi0 += (hsum_lo0 >> DOT_SPLIT_BITS);
        hsum_lo0 &= DOT_SPLIT_MASK;

        hsum_lo1 += (ulong)vec1[i] * (ulong)vec2_1[i];
        hsum_hi1 += (hsum_lo1 >> DOT_SPLIT_BITS);
        hsum_lo1 &= DOT_SPLIT_MASK;

        hsum_lo2 += (ulong)vec1[i] * (ulong)vec2_2[i];
        hsum_hi2 += (hsum_lo2 >> DOT_SPLIT_BITS);
        hsum_lo2 &= DOT_SPLIT_MASK;

        hsum_lo3 += (ulong)vec1[i] * (ulong)vec2_3[i];
        hsum_hi3 += (hsum_lo3 >> DOT_SPLIT_BITS);
        hsum_lo3 &= DOT_SPLIT_MASK;
    }

    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
    NMOD_RED(*res2, pow2_precomp * hsum_hi2 + hsum_lo2, mod);
    NMOD_RED(*res3, pow2_precomp * hsum_hi3 + hsum_lo3, mod);
}
#endif
