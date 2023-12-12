#include <flint/machine_vectors.h>
#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** res[0] = dot(a1, b), res[1] = dot(a2, b)                  */
/** power_two = 2^45 mod p, p2 = (p,p), pinv2 = (1/p,1/p)     */
/*------------------------------------------------------------*/
void _nmod_vec_dot2_small_modulus(mp_ptr res, mp_ptr a1, mp_ptr a2, mp_ptr b, ulong len,
                                  mp_limb_t power_two, vec2d p2, vec2d pinv2)
{
    vec4n sum_low1, sum_high1, low_bits;
    vec4n sum_low2, sum_high2;
    vec4n *a1_4n, *b_4n;
    vec4n *a2_4n;
    mp_limb_t total_low2, total_high2;
    mp_limb_t total_low1, total_high1;
    mp_ptr as1, as2, bs;
    ulong k, num_full_blocks, acc_last1, acc_last2;
    vec2d sum, redsum;
    
    num_full_blocks = len >> 5;
    // mask for 45 low bits
    low_bits = vec4n_set_n((1L << 45) - 1);
    
    a1_4n = (vec4n *) a1;
    a2_4n = (vec4n *) a2;
    b_4n = (vec4n *) b;
    sum_low1 = vec4n_zero();
    sum_high1 = vec4n_zero();
    sum_low2 = vec4n_zero();
    sum_high2 = vec4n_zero();

    // all blocks of length 32
    // each pass in the loop reads 8x4 = 32 coefficients
    for(k = 0; k < num_full_blocks; k++)
    {
        vec4n val_b;
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
       
        sum_high1 = vec4n_add(sum_high1, vec4n_bit_shift_right_45(sum_low1));
        sum_low1 = vec4n_bit_and(sum_low1, low_bits);
        sum_high2 = vec4n_add(sum_high2, vec4n_bit_shift_right_45(sum_low2));
        sum_low2 = vec4n_bit_and(sum_low2, low_bits);
    }
    
    // left with at most 31 coefficients
    // k is now the number of coefficients already read
    for (k = num_full_blocks << 5; k + 3 < len; k += 4)
    {
        vec4n val_b;
        val_b = vec4n_load_unaligned((mp_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((mp_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((mp_ptr) a2_4n++), val_b));
    }

    sum_high1 = vec4n_add(sum_high1, vec4n_bit_shift_right_45(sum_low1));
    sum_low1 = vec4n_bit_and(sum_low1, low_bits);
    sum_high2 = vec4n_add(sum_high2, vec4n_bit_shift_right_45(sum_low2));
    sum_low2 = vec4n_bit_and(sum_low2, low_bits);

    // left with at most 3 coefficients
    as1 = (mp_ptr) a1_4n;
    as2 = (mp_ptr) a2_4n;
    bs = (mp_ptr) b_4n;
    acc_last1 = 0;
    acc_last2 = 0;
    for (; k < len; k++)
    {
        mp_limb_t val_b;
        val_b = (*bs++);
        acc_last1 += (*as1++) * val_b;
        acc_last2 += (*as2++) * val_b;
    }

    total_low1 = sum_low1[0] + sum_low1[1] + sum_low1[2] + sum_low1[3] + (acc_last1 & ((1L << 45) - 1));
    total_high1 = sum_high1[0] + sum_high1[1] + sum_high1[2] + sum_high1[3] + (acc_last1 >> 45);
    total_low2 = sum_low2[0] + sum_low2[1] + sum_low2[2] + sum_low2[3] + (acc_last2 & ((1L << 45) - 1));
    total_high2 = sum_high2[0] + sum_high2[1] + sum_high2[2] + sum_high2[3] + (acc_last2 >> 45);

    sum[0] = total_low1 + power_two*total_high1;
    sum[1] = total_low2 + power_two*total_high2;
    redsum = vec2d_reduce_to_0n(sum, p2, pinv2);
    res[0] = (mp_limb_t) redsum[0];
    res[1] = (mp_limb_t) redsum[1];
}

