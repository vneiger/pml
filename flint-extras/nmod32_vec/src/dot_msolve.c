#include "nmod32_vec.h"

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

uint _nmod32_vec_dot_msolve_avx2(n32_srcptr vec1, n32_srcptr vec2, slong len, uint PRIME)
{
    uint RED_32 = ((ulong)2<<31) % PRIME;
    uint RED_64 = ((ulong)1<<63) % PRIME;
    RED_64 = (RED_64*2) % PRIME;

    __m256i mask = AVX2SET1_64(MONE32);
    const long quo = LENGTHQ8(len);
    const long rem = LENGTHR8(len);

    ulong acc4x64[8];
    ulong acc64;
    const uint *vec1_cp;
    const uint *vec2_cp;

    vec1_cp=vec1;
    vec2_cp=vec2;

    __m256i acc_low = AVX2SETZERO();
    __m256i acc_high = AVX2SETZERO();

    long i = 0;
    for( ; i < quo; ++i)
    {
        __m256i v1 = AVX2LOADU(vec1_cp);
        __m256i v2 = AVX2LOADU(vec2_cp);
        vec1_cp+=8;
        vec2_cp+=8;
        __m256i prod1=AVX2MUL(v1,v2);
        __m256i prod2= AVX2MUL(AVX2SRLI_64(v1,32),AVX2SRLI_64(v2,32));
        __m256i res1 = AVX2ADD_64(prod1, prod2);

        v1=AVX2LOADU(vec1_cp);
        v2=AVX2LOADU(vec2_cp);
        vec1_cp+=8;
        vec2_cp+=8;
        prod1=AVX2MUL(v1,v2);
        prod2= AVX2MUL(AVX2SRLI_64(v1,32),AVX2SRLI_64(v2,32));
        res1 = AVX2ADD_64(res1, AVX2ADD_64(prod1, prod2));

        v1=AVX2LOADU(vec1_cp);
        v2=AVX2LOADU(vec2_cp);
        vec1_cp+=8;
        vec2_cp+=8;
        prod1=AVX2MUL(v1,v2);
        prod2= AVX2MUL(AVX2SRLI_64(v1,32),AVX2SRLI_64(v2,32));
        res1 = AVX2ADD_64(res1, AVX2ADD_64(prod1, prod2));

        v1=AVX2LOADU(vec1_cp);
        v2=AVX2LOADU(vec2_cp);
        vec1_cp+=8;
        vec2_cp+=8;
        prod1=AVX2MUL(v1,v2);
        prod2= AVX2MUL(AVX2SRLI_64(v1,32),AVX2SRLI_64(v2,32));
        res1 = AVX2ADD_64(res1, AVX2ADD_64(prod1, prod2));

        acc_low=AVX2ADD_64(acc_low,AVX2AND_(res1,mask));
        acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(res1,32));
    }

    AVX2STOREU(acc4x64,acc_low);
    AVX2STOREU(acc4x64+4,acc_high);

    /* Reduction */
    acc64=0;
    for(i=0;i<4;++i){
        //partie haute du registre haut (2^64 ->2^95)
        acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
        //partie basse du registre haut (2^32->2^63)
        acc4x64[i]+=((acc4x64[i+4]&((ulong)0xFFFFFFFF))*RED_32)%PRIME;
        acc64+=acc4x64[i]%PRIME;
    }

    /* *vec_res[j]=MODRED32(acc64, PRIME, preinv); */
    uint res = acc64 % PRIME;

    long tmp = 0;
    for(long k = 0; k < rem; k++)
        tmp += ((long)vec1_cp[k] * (long)vec2_cp[k]) % PRIME;

    res = (res + tmp) % PRIME;

    return res;
}

void nmod_mat_mul_nmod_vec_msolve32_avx2(n32_ptr vec_res, n32_srcptr mat, n32_srcptr vec,
                                         uint ncols, uint nrows, uint PRIME)
{
    uint RED_32 = ((ulong)2<<31) % PRIME;
    uint RED_64 = ((ulong)1<<63) % PRIME;
    RED_64 = (RED_64*2) % PRIME;

    //mask pour recuperer les parties basses
    __m256i mask = AVX2SET1_64(MONE32);
    const long quo = LENGTHQ8(ncols);
    const long rem = LENGTHR8(ncols);

    unsigned int i,j;
    ulong acc4x64[8];
    ulong acc64;
    __m256i acc_low,acc_high;
    __m256i prod1,prod2;
    __m256i res1;
    const uint *vec_cp;
    const uint *mat_cp;

    /* For each row of the matrix, we compute a dot product */
    for(j = 0; j < nrows; ++j){
        vec_cp=vec;
        mat_cp=mat + j*ncols;

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
            //partie haute du registre haut (2^64 ->2^95)
            acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
            //partie basse du registre haut (2^32->2^63)
            acc4x64[i]+=((acc4x64[i+4]&((ulong)0xFFFFFFFF))*RED_32)%PRIME;
            acc64+=acc4x64[i]%PRIME;
        }

        /* *vec_res[j]=MODRED32(acc64, PRIME, preinv); */
        vec_res[j]=acc64%PRIME;

        long tmp = 0;
        for(long k = 0; k < rem; k++){
            tmp += ((long)mat_cp[k] * (long)vec_cp[k]) % PRIME;
        }
        vec_res[j] = (vec_res[j] + tmp) % PRIME;
    }
}
