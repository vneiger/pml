#include "nmod_extra.h"
#include "nmod_vec_extra.h"
#include <flint/ulong_extras.h>

/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** returns dot(a, b)                                         */
/** power_two = 2^45 mod p, pinv = 1/p                        */
/*------------------------------------------------------------*/
ulong _nmod_vec_dot_small_modulus(nn_ptr a, nn_ptr b, ulong len,
                                      ulong power_two, vec1d p, vec1d pinv)
{
    const ulong num_full_blocks = len >> 5;
    // mask for 45 low bits
    const vec4n low_bits = vec4n_set_n((1L << 45) - 1);

    vec4n * a_4n = (vec4n *) a;
    vec4n * b_4n = (vec4n *) b;
    vec4n sum_low = vec4n_zero();
    vec4n sum_high = vec4n_zero();
    vec4n tmp;

    // all blocks of length 32
    // each pass in the loop reads 8x4 = 32 coefficients
    for (ulong k = 0; k < num_full_blocks; k++)
    {
        tmp = sum_low;
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));

        sum_low = vec4n_bit_and(tmp, low_bits);
        sum_high = vec4n_add(sum_high, vec4n_bit_shift_right_45(tmp));
    }

    // left with at most 31 coefficients
    // k is now the number of coefficients already read
    ulong k = num_full_blocks << 5;
    for (; k + 3 < len; k += 4)
        sum_low = vec4n_add(sum_low, vec4n_mul(vec4n_load_unaligned((nn_ptr) a_4n++), vec4n_load_unaligned((nn_ptr) b_4n++)));

    sum_high = vec4n_add(sum_high, vec4n_bit_shift_right_45(sum_low));
    sum_low = vec4n_bit_and(sum_low, low_bits);

    // left with at most 3 coefficients
    nn_ptr as = (nn_ptr) a_4n;
    nn_ptr bs = (nn_ptr) b_4n;
    ulong acc_last = 0;
    for (; k < len; k++)
        acc_last += (*as++) * (*bs++);

    const ulong total_low = sum_low[0] + sum_low[1] + sum_low[2] + sum_low[3] + (acc_last & ((UWORD(1) << 45) - 1));
    const ulong total_high = sum_high[0] + sum_high[1] + sum_high[2] + sum_high[3] + (acc_last >> 45);

    vec1d sum = total_low + power_two * total_high;
    return (ulong) vec1d_reduce_to_0n(sum, p, pinv);
}

// dot product using single limb, avx2
ulong _nmod_vec_dot_product_1_avx2(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod)
{
    // compute 4 vertical sub-dot products
    __m256i res_vec = _mm256_setzero_si256();
    ulong i;
    for (i = 0; i+3 < len; i += 4)
    {
        // load + multiplication + addition
        __m256i x = _mm256_loadu_si256((__m256i *) (vec1+i));
        __m256i y = _mm256_loadu_si256((__m256i *) (vec2+i));
        x = _mm256_mul_epu32(x, y);
        res_vec = _mm256_add_epi64(res_vec, x);
    }

    // horizontal add
    ulong res = res_vec[0] + res_vec[1] + res_vec[2] + res_vec[3];

    // scalar loop for leftover entries
    for (; i < len; ++i)
        res += vec1[i] * vec2[i];
    NMOD_RED(res, res, mod);

    return res;
}

// dot product using single limb, avx512
ulong _nmod_vec_dot_product_1_avx512(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod) {
    // compute 4 vertical sub-dot products
    __m512i res_vec = _mm512_setzero_si512();
    ulong i;
    for (i = 0; i+7 < len; i += 8)
    {
        // load + multiplication + addition
        __m512i x = _mm512_loadu_si512((__m512i *) (vec1+i));
        __m512i y = _mm512_loadu_si512((__m512i *) (vec2+i));
        x = _mm512_mul_epu32(x, y);
        res_vec = _mm512_add_epi64(res_vec, x);
    }

    // horizontal add
    //ulong res = res_vec[0] + res_vec[1] + res_vec[2] + res_vec[3] 
              //+ res_vec[4] + res_vec[5] + res_vec[6] + res_vec[7];
    ulong res = _mm512_reduce_add_epi64(res_vec);

    // scalar loop for leftover entries
    for (; i < len; ++i)
        res += vec1[i] * vec2[i];
    NMOD_RED(res, res, mod);

    return res;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
