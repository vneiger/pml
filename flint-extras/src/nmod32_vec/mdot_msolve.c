/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "nmod32_vec.h"

#if PML_HAVE_AVX2

#define MONE32 ((uint)0xFFFFFFFF)
/* Euclidean division of LENGTH by 8 */
/* LENGTH/8 */
#define LENGTHQ(l) (l>>3)
/* LENGTH%8 */
#define LENGTHR(l) (l&7)
#define LENGTHQ8(l) (l / 32)
#define LENGTHR8(l) (l%32)

#define LENGTHQ14(l) (l / 56)
#define LENGTHR14(l) (l%56)

#define AVX2LOAD(A) _mm256_load_si256((__m256i*)(A))
#define AVX2LOADU(A) _mm256_loadu_si256((__m256i*)(A))
#define AVX2STORE(res,A) _mm256_store_si256((__m256i*)(res),A);
#define AVX2STOREU(res,A) _mm256_storeu_si256((__m256i*)(res),A);
#define AVX2SETZERO() _mm256_setzero_si256()
#define AVX2SET1_64(A0) _mm256_set1_epi64x(A0)
#define AVX2AND_(A,B) _mm256_and_si256(A,B)
#define AVX2ADD_64(A,B) _mm256_add_epi64(A,B)
#define AVX2MUL(A,B) _mm256_mul_epu32(A,B)
#define AVX2SRLI_64(A,i) _mm256_srli_epi64(A,i)


void _nmod32_vec_mdot_msolve_via_dot_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                          slong nrows, slong len, slong stride, nmod_t mod)
{
    const uint PRIME = mod.n;
    for (slong i = 0; i < nrows; i++)
        mv[i] = _nmod32_vec_dot_msolve_avx2(mat + i*stride, vec, len, PRIME);
}


void _nmod32_vec_mdot_msolve_native_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                         slong nrows, slong len, slong stride, nmod_t mod)
{
    const uint PRIME = mod.n;
    uint RED_32 = ((ulong)2<<31) % PRIME;
    uint RED_64 = ((ulong)1<<63) % PRIME;
    RED_64 = (RED_64*2) % PRIME;

    //mask pour recuperer les parties basses
    __m256i mask = AVX2SET1_64(MONE32);
    const long quo = LENGTHQ8(len);
    const long rem = LENGTHR8(len);

    unsigned int i,j;
    ulong acc4x64[8];
    ulong acc64;
    __m256i acc_low,acc_high;
    __m256i prod1,prod2;
    __m256i res1;
    const uint *vec_cp;
    const uint *mat_cp;

    /* For each row of the matrix, we compute a dot product */
    for(j = 0; j < nrows; ++j)
    {
        vec_cp=vec;
        mat_cp=mat + j*stride;

        acc_low=AVX2SETZERO();
        acc_high=AVX2SETZERO();

        const long bs = quo;
        /* Accumulation */
        for(i = 0 ; i < bs; ++i){

            __m256i mat8=AVX2LOADU(mat_cp);
            __m256i vec8=AVX2LOADU(vec_cp);
            /* four 32-bits mul, lower parts */
            /* four 32-bits mul, higher parts */
            prod1=AVX2MUL(mat8,vec8);
            prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            res1 = prod1 + prod2;
            mat_cp+=8;
            vec_cp+=8;
            mat8=AVX2LOADU(mat_cp);
            vec8=AVX2LOADU(vec_cp);
            prod1=AVX2MUL(mat8,vec8);
            prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            res1 += prod1 + prod2;
            mat_cp+=8;
            vec_cp+=8;
            mat8=AVX2LOADU(mat_cp);
            vec8=AVX2LOADU(vec_cp);
            prod1=AVX2MUL(mat8,vec8);
            prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            res1 += prod1 + prod2;
            mat_cp+=8;
            vec_cp+=8;
            mat8=AVX2LOADU(mat_cp);
            vec8=AVX2LOADU(vec_cp);
            prod1=AVX2MUL(mat8,vec8);
            prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            res1 += prod1 + prod2;
            mat_cp+=8;
            vec_cp+=8;

            acc_low=AVX2ADD_64(acc_low,AVX2AND_(res1,mask));
            acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(res1,32));
        }

        AVX2STOREU(acc4x64,acc_low);
        AVX2STOREU(acc4x64+4,acc_high);

        /* Reduction */
        acc64=0;
        for(i=0;i<4;++i){
            // high part (2^64->2^95)
            acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
            // low part (2^32->2^63)
            acc4x64[i]+=((acc4x64[i+4]&((ulong)0xFFFFFFFF))*RED_32)%PRIME;
            acc64+=acc4x64[i]%PRIME;
        }

        mv[j]=acc64%PRIME;

        long tmp = 0;
        for(long k = 0; k < rem; k++){
            tmp += ((long)mat_cp[k] * (long)vec_cp[k]) % PRIME;
        }
        mv[j] = (mv[j] + tmp) % PRIME;
    }
}

#endif  /* PML_HAVE_AVX2 */
