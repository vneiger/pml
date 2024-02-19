#include "nmod_extra.h"
#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** returns dot(a1, b)                                        */
/** power_two = 2^45 mod p, pinv = 1/p                        */
/*------------------------------------------------------------*/
mp_limb_t _nmod_vec_dot_small_modulus(mp_ptr a, mp_ptr b, ulong len,
                                      mp_limb_t power_two, vec1d p, vec1d pinv)
{
    vec4n sum_low, sum_high, tmp, low_bits;
    vec4n *a_4n, *b_4n;
    mp_limb_t total_low, total_high;
    mp_ptr as, bs;
    ulong k, num_full_blocks, acc_last;
    vec1d sum;

    num_full_blocks = len >> 5;
    // mask for 45 low bits
    low_bits = vec4n_set_n((1L << 45) - 1);

    a_4n = (vec4n *) a;
    b_4n = (vec4n *) b;
    sum_low = vec4n_zero();
    sum_high = vec4n_zero();

    // all blocks of length 32
    // each pass in the loop reads 8x4 = 32 coefficients
    for(k = 0; k < num_full_blocks; k++)
    {
        tmp = sum_low;
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));
        tmp = vec4n_add(tmp, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));

        sum_low = vec4n_bit_and(tmp, low_bits);
        sum_high = vec4n_add(sum_high, vec4n_bit_shift_right_45(tmp));
    }

    // left with at most 31 coefficients
    // k is now the number of coefficients already read
    for (k = num_full_blocks << 5; k + 3 < len; k += 4)
        sum_low = vec4n_add(sum_low, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++)));

    sum_high = vec4n_add(sum_high, vec4n_bit_shift_right_45(sum_low));
    sum_low = vec4n_bit_and(sum_low, low_bits);

    // left with at most 3 coefficients
    as = (mp_ptr) a_4n;
    bs = (mp_ptr) b_4n;
    acc_last = 0;
    for (; k < len; k++)
        acc_last += (*as++) * (*bs++);

    total_low = sum_low[0] + sum_low[1] + sum_low[2] + sum_low[3] + (acc_last & ((1L << 45) - 1));
    total_high = sum_high[0] + sum_high[1] + sum_high[2] + sum_high[3] + (acc_last >> 45);

    sum = total_low + power_two * total_high;
    return (mp_limb_t) vec1d_reduce_to_0n(sum, p, pinv);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
