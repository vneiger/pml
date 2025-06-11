/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h> // for NMOD_RED
#include <flint/nmod_vec.h>  // for DOT_SPLIT_MASK

#include "nmod32_vec.h"


void _nmod32_vec_dot2_split(uint * res0, uint * res1, n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, slong len, nmod_t mod, ulong pow2_precomp)
{
    ulong dp_lo0 = 0;
    ulong dp_lo1 = 0;
    ulong dp_hi0 = 0;
    ulong dp_hi1 = 0;

    slong i = 0;

    for ( ; i+7 < len; i+=8)
    {
        dp_lo0 += (ulong)vec1[i+0] * vec2_0[i+0]
                + (ulong)vec1[i+1] * vec2_0[i+1]
                + (ulong)vec1[i+2] * vec2_0[i+2]
                + (ulong)vec1[i+3] * vec2_0[i+3]
                + (ulong)vec1[i+4] * vec2_0[i+4]
                + (ulong)vec1[i+5] * vec2_0[i+5]
                + (ulong)vec1[i+6] * vec2_0[i+6]
                + (ulong)vec1[i+7] * vec2_0[i+7];

        dp_lo1 += (ulong)vec1[i+0] * vec2_1[i+0]
                + (ulong)vec1[i+1] * vec2_1[i+1]
                + (ulong)vec1[i+2] * vec2_1[i+2]
                + (ulong)vec1[i+3] * vec2_1[i+3]
                + (ulong)vec1[i+4] * vec2_1[i+4]
                + (ulong)vec1[i+5] * vec2_1[i+5]
                + (ulong)vec1[i+6] * vec2_1[i+6]
                + (ulong)vec1[i+7] * vec2_1[i+7];

        dp_hi0 += (dp_lo0 >> DOT_SPLIT_BITS);
        dp_hi1 += (dp_lo1 >> DOT_SPLIT_BITS);
        dp_lo0 &= DOT_SPLIT_MASK;
        dp_lo1 &= DOT_SPLIT_MASK;
    }

    // less than 8 terms remaining, can be accumulated
    for ( ; i < len; i++)
    {
        dp_lo0 += (ulong)vec1[i] * vec2_0[i];
        dp_lo1 += (ulong)vec1[i] * vec2_1[i];
    }
    dp_hi0 += (dp_lo0 >> DOT_SPLIT_BITS);
    dp_hi1 += (dp_lo1 >> DOT_SPLIT_BITS);
    dp_lo0 &= DOT_SPLIT_MASK;
    dp_lo1 &= DOT_SPLIT_MASK;

    NMOD_RED(*res0, pow2_precomp * dp_hi0 + dp_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * dp_hi1 + dp_lo1, mod);
}

#if PML_HAVE_AVX2
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

    for ( ; i+31 < len; i+=32)
    {
        __m256i v1_0   = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 0));
        __m256i v1_1   = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 8));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        __m256i v2_0_1 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 8));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        __m256i v2_1_1 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 8));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        __m256i dp_lo2 = _mm256_mul_epu32(v1_1, v2_0_1);
        __m256i dp_lo3 = _mm256_mul_epu32(v1_1, v2_1_1);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm256_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm256_shuffle_epi32(v2_1_1, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        v1_0   = _mm256_loadu_si256((const __m256i *) (vec1  +i+16));
        v1_1   = _mm256_loadu_si256((const __m256i *) (vec1  +i+24));
        v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+16));
        v2_0_1 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+24));
        v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+16));
        v2_1_1 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+24));

        // 3rd term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        // 4th term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm256_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm256_shuffle_epi32(v2_1_1, 0xB1);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        // gather results in dp_lo0 and dp_lo1
        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo2);
        dp_lo1 = _mm256_add_epi64(dp_lo1, dp_lo3);

        // split
        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
        dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo0 and dp_lo1
    for ( ; i+7 < len; i+=8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
    }

    dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
    dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
    dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);

    ulong hsum_lo0 = _mm256_hsum(dp_lo0);
    ulong hsum_hi0 = _mm256_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm256_hsum(dp_lo1);
    ulong hsum_hi1 = _mm256_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    // less than 8 terms remaining, can accumulate
    for (; i < len; i++)
    {
        hsum_lo0 += (ulong)vec1[i] * vec2_0[i];
        hsum_lo1 += (ulong)vec1[i] * vec2_1[i];
    }
    hsum_hi0 += (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    hsum_hi1 += (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    // the requirement on "len <= DOT2_ACC8_MAX_LEN"
    // ensures pow2_precomp * hsum_hi + hsum_lo fits in 64 bits
    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}

void _nmod32_vec_dot3_split_avx2(uint * res, n32_srcptr vec1,
                                 n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                 slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);

    __m256i dp_lo[3];
    __m256i dp_hi[3];
    __m256i v1;
    __m256i v2[3];

    dp_lo[0] = _mm256_setzero_si256();
    dp_lo[1] = _mm256_setzero_si256();
    dp_lo[2] = _mm256_setzero_si256();

    dp_hi[0] = _mm256_setzero_si256();
    dp_hi[1] = _mm256_setzero_si256();
    dp_hi[2] = _mm256_setzero_si256();

    slong i = 0;

    for ( ; i+31 < len; i+=32)
    {
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 0));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 0));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 3rd+4th term
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 8));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 8));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 8));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 8));

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 5th+6th term
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+16));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+16));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+16));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+16));

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 7th+8th term
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+24));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+24));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+24));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+24));

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // split
        dp_hi[0] = _mm256_add_epi64(dp_hi[0], _mm256_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
        dp_hi[1] = _mm256_add_epi64(dp_hi[1], _mm256_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
        dp_hi[2] = _mm256_add_epi64(dp_hi[2], _mm256_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
        dp_lo[0] = _mm256_and_si256(dp_lo[0], low_bits);
        dp_lo[1] = _mm256_and_si256(dp_lo[1], low_bits);
        dp_lo[2] = _mm256_and_si256(dp_lo[2], low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms
    for ( ; i+7 < len; i += 8)
    {
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 0));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 0));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));
    }

    // split
    dp_hi[0] = _mm256_add_epi64(dp_hi[0], _mm256_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
    dp_hi[1] = _mm256_add_epi64(dp_hi[1], _mm256_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
    dp_hi[2] = _mm256_add_epi64(dp_hi[2], _mm256_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
    dp_lo[0] = _mm256_and_si256(dp_lo[0], low_bits);
    dp_lo[1] = _mm256_and_si256(dp_lo[1], low_bits);
    dp_lo[2] = _mm256_and_si256(dp_lo[2], low_bits);

    ulong hsum_lo[3];
    ulong hsum_hi[3];
    hsum_lo[0] = _mm256_hsum(dp_lo[0]);
    hsum_hi[0] = _mm256_hsum(dp_hi[0]) + (hsum_lo[0] >> DOT_SPLIT_BITS);
    hsum_lo[0] &= DOT_SPLIT_MASK;
    hsum_lo[1] = _mm256_hsum(dp_lo[1]);
    hsum_hi[1] = _mm256_hsum(dp_hi[1]) + (hsum_lo[1] >> DOT_SPLIT_BITS);
    hsum_lo[1] &= DOT_SPLIT_MASK;
    hsum_lo[2] = _mm256_hsum(dp_lo[2]);
    hsum_hi[2] = _mm256_hsum(dp_hi[2]) + (hsum_lo[2] >> DOT_SPLIT_BITS);
    hsum_lo[2] &= DOT_SPLIT_MASK;

    // less than 8 terms remaining, can accumulate
    for (; i < len; i++)
    {
        hsum_lo[0] += (ulong)vec1[i] * vec2_0[i];
        hsum_lo[1] += (ulong)vec1[i] * vec2_1[i];
        hsum_lo[2] += (ulong)vec1[i] * vec2_2[i];
    }

    hsum_hi[0] += (hsum_lo[0] >> DOT_SPLIT_BITS);
    hsum_hi[1] += (hsum_lo[1] >> DOT_SPLIT_BITS);
    hsum_hi[2] += (hsum_lo[2] >> DOT_SPLIT_BITS);
    hsum_lo[0] &= DOT_SPLIT_MASK;
    hsum_lo[1] &= DOT_SPLIT_MASK;
    hsum_lo[2] &= DOT_SPLIT_MASK;

    NMOD_RED(res[0], pow2_precomp * hsum_hi[0] + hsum_lo[0], mod);
    NMOD_RED(res[1], pow2_precomp * hsum_hi[1] + hsum_lo[1], mod);
    NMOD_RED(res[2], pow2_precomp * hsum_hi[2] + hsum_lo[2], mod);
}
#endif  /* PML_HAVE_AVX2 */


#if PML_HAVE_AVX512
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

    for ( ; i+63 < len; i+=64)
    {
        __m512i v1_0   = _mm512_loadu_si512((const __m512i *) (vec1  +i+ 0));
        __m512i v1_1   = _mm512_loadu_si512((const __m512i *) (vec1  +i+16));
        __m512i v2_0_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+ 0));
        __m512i v2_0_1 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+16));
        __m512i v2_1_0 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+ 0));
        __m512i v2_1_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+16));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        __m512i dp_lo2 = _mm512_mul_epu32(v1_1, v2_0_1);
        __m512i dp_lo3 = _mm512_mul_epu32(v1_1, v2_1_1);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm512_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm512_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm512_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm512_shuffle_epi32(v2_1_1, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        v1_0   = _mm512_loadu_si512((const __m512i *) (vec1  +i+32));
        v1_1   = _mm512_loadu_si512((const __m512i *) (vec1  +i+48));
        v2_0_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+32));
        v2_0_1 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+48));
        v2_1_0 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+32));
        v2_1_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+48));

        // 3rd term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        // 4th term: high 32 bit word of each 64 bit word
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm512_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm512_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm512_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm512_shuffle_epi32(v2_1_1, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        // gather results in dp_lo0 and dp_lo1
        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo2);
        dp_lo1 = _mm512_add_epi64(dp_lo1, dp_lo3);

        // split
        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
        dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo0 and dp_lo1
    for ( ; i+15 < len; i+=16)
    {
        __m512i v1   = _mm512_loadu_si512((const __m512i *) (vec1  +i));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i));
        __m512i v2_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));

        // 2nd term: high 32 bit word of each 64 bit word
        v1 = _mm512_shuffle_epi32(v1, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));
    }

    // finally, do a last iteration which may be "incomplete": use mask load
    // (at each position we accumulate 0 or 2 terms, so including the above
    // loop this is a total of <= 8 terms)
    if (i < len)
    {
        // mask == 0b0...01...1 with number of 1's == number of remaining terms
        __mmask16 mask = (1 << (len-i)) - 1;
        __m512i v1   = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec1  +i));
        __m512i v2_0 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_0+i));
        __m512i v2_1 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_1+i));

        // 1st term
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));

        // 2nd term
        v1   = _mm512_shuffle_epi32(v1  , 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));
    }

    // split
    dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
    dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
    dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);

    // gather 8 terms in single ulong
    ulong hsum_lo0 = _mm512_hsum(dp_lo0);
    ulong hsum_hi0 = _mm512_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm512_hsum(dp_lo1);
    ulong hsum_hi1 = _mm512_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}

void _nmod32_vec_dot3_split_avx512(uint * res,
                                   n32_srcptr vec1, n32_srcptr vec2_0, n32_srcptr vec2_1, n32_srcptr vec2_2,
                                   slong len, nmod_t mod, ulong pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);

    __m512i dp_lo[3];
    __m512i dp_hi[3];
    __m512i v1[2];
    __m512i v2[6];

    dp_lo[0] = _mm512_setzero_si512();
    dp_lo[1] = _mm512_setzero_si512();
    dp_lo[2] = _mm512_setzero_si512();

    dp_hi[0] = _mm512_setzero_si512();
    dp_hi[1] = _mm512_setzero_si512();
    dp_hi[2] = _mm512_setzero_si512();

    slong i = 0;

    for ( ; i+63 < len; i+=64)
    {
        v1[0] = _mm512_loadu_si512((const __m512i *) (vec1  +i+ 0));
        v1[1] = _mm512_loadu_si512((const __m512i *) (vec1  +i+16));
        v2[0] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+ 0));
        v2[1] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+16));
        v2[2] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+ 0));
        v2[3] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+16));
        v2[4] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+ 0));
        v2[5] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+16));

        // 1st+3rd term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        // 2nd+4th term: high 32 bit word of each 64 bit word
        v1[0] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v1[1] = _mm512_shuffle_epi32(v1[1], 0xB1);
        v2[0] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm512_shuffle_epi32(v2[2], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[3], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[4], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[5], 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm512_srli_epi64(v1_0, 32)

        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        // 5th, 6th, 7th, 8th terms
        v1[0] = _mm512_loadu_si512((const __m512i *) (vec1  +i+32));
        v1[1] = _mm512_loadu_si512((const __m512i *) (vec1  +i+48));
        v2[0] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+32));
        v2[1] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+48));
        v2[2] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+32));
        v2[3] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+48));
        v2[4] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+32));
        v2[5] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+48));

        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        v1[0] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v1[1] = _mm512_shuffle_epi32(v1[1], 0xB1);
        v2[0] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm512_shuffle_epi32(v2[2], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[3], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[4], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[5], 0xB1);

        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        // split
        dp_hi[0] = _mm512_add_epi64(dp_hi[0], _mm512_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
        dp_lo[0] = _mm512_and_si512(dp_lo[0], low_bits);
        dp_hi[1] = _mm512_add_epi64(dp_hi[1], _mm512_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
        dp_lo[1] = _mm512_and_si512(dp_lo[1], low_bits);
        dp_hi[2] = _mm512_add_epi64(dp_hi[2], _mm512_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
        dp_lo[2] = _mm512_and_si512(dp_lo[2], low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo
    for ( ; i+15 < len; i+=16)
    {
        v1[0] = _mm512_loadu_si512((const __m512i *) (vec1  +i));
        v2[0] = _mm512_loadu_si512((const __m512i *) (vec2_0+i));
        v2[1] = _mm512_loadu_si512((const __m512i *) (vec2_1+i));
        v2[2] = _mm512_loadu_si512((const __m512i *) (vec2_2+i));
        v1[1] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[2], 0xB1);

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[1]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[4]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));
    }

    // finally, do a last iteration which may be "incomplete": use mask load
    // (at each position we accumulate 0 or 2 terms, so including the above
    // loop this is a total of <= 8 terms)
    if (i < len)
    {
        // mask == 0b0...01...1 with number of 1's == number of remaining terms
        __mmask16 mask = (1 << (len-i)) - 1;
        v1[0] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec1  +i));
        v2[0] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_0+i));
        v2[1] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_1+i));
        v2[2] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_2+i));
        v1[1] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[2], 0xB1);

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[1]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[4]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));
    }

    // split
    dp_hi[0] = _mm512_add_epi64(dp_hi[0], _mm512_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
    dp_lo[0] = _mm512_and_si512(dp_lo[0], low_bits);
    dp_hi[1] = _mm512_add_epi64(dp_hi[1], _mm512_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
    dp_lo[1] = _mm512_and_si512(dp_lo[1], low_bits);
    dp_hi[2] = _mm512_add_epi64(dp_hi[2], _mm512_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
    dp_lo[2] = _mm512_and_si512(dp_lo[2], low_bits);

    // gather 8 terms in single ulong
    ulong hsum_lo0 = _mm512_hsum(dp_lo[0]);
    ulong hsum_hi0 = _mm512_hsum(dp_hi[0]) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = _mm512_hsum(dp_lo[1]);
    ulong hsum_hi1 = _mm512_hsum(dp_hi[1]) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;
    ulong hsum_lo2 = _mm512_hsum(dp_lo[2]);
    ulong hsum_hi2 = _mm512_hsum(dp_hi[2]) + (hsum_lo2 >> DOT_SPLIT_BITS);
    hsum_lo2 &= DOT_SPLIT_MASK;

    NMOD_RED(res[0], pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(res[1], pow2_precomp * hsum_hi1 + hsum_lo1, mod);
    NMOD_RED(res[2], pow2_precomp * hsum_hi2 + hsum_lo2, mod);
}
#endif  /* PML_HAVE_AVX512 */
