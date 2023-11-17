#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


mp_limb_t _nmod_vec_dot_small_modulus(mp_ptr a, mp_ptr b, ulong len, nmod_t mod)
{
    vec4n sum_low, sum_high, tmp, low_bits;
    vec4n *a_4n, *b_4n;
    mp_limb_t total_low, total_high, two_32;;
    mp_ptr as, bs;
    mp_limb_t values[4] __attribute__((aligned(32)));
    ulong k, num_full_blocks, acc_last;
    
    num_full_blocks = len >> 5;
    // mask for 32 low bits
    low_bits = vec4n_set_n(4294967295);
    
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
        sum_high = vec4n_add(sum_high, vec4n_bit_shift_right_32(tmp));
    }

    // left with at most 31 coefficients
    // k is now the number of coefficients already read
    for (k = num_full_blocks << 5; k + 4 < len; k += 4)
        sum_low = vec4n_add(sum_low, vec4n_mul(vec4n_load_unaligned((mp_ptr) a_4n++), vec4n_load_unaligned((mp_ptr) b_4n++))); 

    sum_high = vec4n_add(sum_high, vec4n_bit_shift_right_32(sum_low));
    sum_low = vec4n_bit_and(sum_low, low_bits);

    // left with at most 3 coefficients
    as = (mp_ptr) a_4n;
    bs = (mp_ptr) b_4n;
    acc_last = 0;
    for (; k < len; k++)
        acc_last += (*as++) * (*bs++);
    
    vec4n_store_aligned(values, sum_low);
    total_low = values[0] + values[1] + values[2] + values[3] + (acc_last & 4294967295);

    vec4n_store_aligned(values, sum_high);
    total_high = values[0] + values[1] + values[2] + values[3] + (acc_last >> 32);

    NMOD_RED(total_low, total_low, mod);
    NMOD_RED(total_high, total_high, mod);
    NMOD_RED(two_32, 1L << 32, mod);
    NMOD_ADDMUL(total_low, total_high, two_32, mod);
    return total_low;
}
