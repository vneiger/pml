#include "flint/nmod_vec.h"  // for DOT_SPLIT_MASK
#include "nmod32_vec.h"
#include <immintrin.h>

// horizontal sum
FLINT_FORCE_INLINE ulong _mm256_hsum(__m256i a) {
    __m256i a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    __m256i sum_lo = _mm256_add_epi64(a, a_hi);
    __m128i sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (ulong) _mm_cvtsi128_si64(sum);
}

uint _nmod32_vec_dot2_split_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_hi = _mm256_setzero_si256();

    slong i = 0;
    for ( ; i+15 < len; i += 16)
    {
        __m256i dp_lo1 = _mm256_setzero_si256();
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i v1_1 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 8));
        __m256i v2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 0));
        __m256i v2_1 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 8));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_1, v2_1));

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0 = _mm256_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm256_shuffle_epi32(v2_1, 0xB1);
        // alternative 2: vpsrlq
        //v1_0 = _mm256_srli_epi64(v1_0, 32);
        //v2_0 = _mm256_srli_epi64(v2_0, 32);
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_1, v2_1));

        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo1);
        dp_hi = _mm256_add_epi64(dp_hi, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
    }

    for ( ; i+7 < len; i += 8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i v2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 0));

        // handle low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));

        // handle high 32 bit word of each 64 bit word
        // alternative 1: vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm256_shuffle_epi32(v2_0, 0xB1);
        // alternative 2: vpsrlq (tested)
        //v1_0 = _mm256_srli_epi64(v1_0, 32);
        //v2_0 = _mm256_srli_epi64(v2_0, 32);
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));

        dp_hi = _mm256_add_epi64(dp_hi, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
    }
    //dp_hi = _mm256_add_epi64(dp_hi, vec4n_bit_shift_right(dp_lo0, DOT_SPLIT_BITS));
    //dp_lo0 = vec4n_bit_and(dp_lo0, low_bits);

    // from here on, not really careful about operations that could be done 
    // with int rather than long
    ulong hsum_lo = _mm256_hsum(dp_lo0);
    ulong hsum_hi = _mm256_hsum(dp_hi) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // some accumulation could be done here (but not up to 7 terms)
    for (; i < len; i++)
    {
        hsum_lo += (ulong)vec1[i] * (ulong)vec2[i];
        hsum_hi += (hsum_lo >> DOT_SPLIT_BITS);
        hsum_lo &= DOT_SPLIT_MASK;
    }

    ulong res;
    NMOD_RED(res, pow2_precomp * hsum_hi + hsum_lo, mod);
    return (uint)res;
}

