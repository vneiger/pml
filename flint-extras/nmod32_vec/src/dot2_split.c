#include "nmod32_vec.h"

uint _nmod_vec_dot2_split(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
{
    const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    slong i = 0;
    for ( ; i+31 < len; i += 32)
    {
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_load_unaligned(vec2+i+ 0)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_load_unaligned(vec2+i+ 4)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_load_unaligned(vec2+i+ 8)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_load_unaligned(vec2+i+12)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_load_unaligned(vec2+i+16)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_load_unaligned(vec2+i+20)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_load_unaligned(vec2+i+24)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_load_unaligned(vec2+i+28)));

        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    for ( ; i + 3 < len; i += 4)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_load_unaligned(vec2+i)));

    dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
    dp_lo = vec4n_bit_and(dp_lo, low_bits);

    ulong hsum_lo = vec4n_horizontal_sum(dp_lo);
    const ulong hsum_hi = vec4n_horizontal_sum(dp_hi) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    for (; i < len; i++)
        hsum_lo += vec1[i] * vec2[i];

    ulong res;
    NMOD_RED(res, pow2_precomp * hsum_hi + hsum_lo, mod);
    return res;
}

