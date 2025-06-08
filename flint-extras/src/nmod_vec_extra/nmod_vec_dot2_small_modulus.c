#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

#if FLINT_BITS == 64
void _nmod_vec_2dot2_split(nn_ptr res, nn_srcptr vec10, nn_srcptr vec11, nn_srcptr vec2,
                           slong len, nmod_t mod, ulong pow2_precomp)
#if PML_HAVE_AVX2
{
    const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
    vec4n dp_lo0 = vec4n_zero();
    vec4n dp_hi0 = vec4n_zero();
    vec4n dp_lo1 = vec4n_zero();
    vec4n dp_hi1 = vec4n_zero();

    slong i = 0;
    for ( ; i+31 < len; i += 32)
    {
        vec4n v2;

        v2 = vec4n_load_unaligned(vec2+i+ 0);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+ 0), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+ 0), v2));

        v2 = vec4n_load_unaligned(vec2+i+ 4);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+ 4), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+ 4), v2));

        v2 = vec4n_load_unaligned(vec2+i+ 8);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+ 8), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+ 8), v2));

        v2 = vec4n_load_unaligned(vec2+i+12);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+12), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+12), v2));

        v2 = vec4n_load_unaligned(vec2+i+16);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+16), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+16), v2));

        v2 = vec4n_load_unaligned(vec2+i+20);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+20), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+20), v2));

        v2 = vec4n_load_unaligned(vec2+i+24);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+24), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+24), v2));

        v2 = vec4n_load_unaligned(vec2+i+28);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i+28), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i+28), v2));

        dp_hi0 = vec4n_add(dp_hi0, vec4n_bit_shift_right(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = vec4n_add(dp_hi1, vec4n_bit_shift_right(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = vec4n_bit_and(dp_lo0, low_bits);
        dp_lo1 = vec4n_bit_and(dp_lo1, low_bits);
    }

    for ( ; i + 3 < len; i += 4)
    {
        vec4n v2 = vec4n_load_unaligned(vec2+i);
        dp_lo0 = vec4n_add(dp_lo0, vec4n_mul(vec4n_load_unaligned(vec10+i), v2));
        dp_lo1 = vec4n_add(dp_lo1, vec4n_mul(vec4n_load_unaligned(vec11+i), v2));
    }

    dp_hi0 = vec4n_add(dp_hi0, vec4n_bit_shift_right(dp_lo0, DOT_SPLIT_BITS));
    dp_hi1 = vec4n_add(dp_hi1, vec4n_bit_shift_right(dp_lo1, DOT_SPLIT_BITS));
    dp_lo0 = vec4n_bit_and(dp_lo0, low_bits);
    dp_lo1 = vec4n_bit_and(dp_lo1, low_bits);

    ulong hsum_lo0 = vec4n_horizontal_sum(dp_lo0);
    const ulong hsum_hi0 = vec4n_horizontal_sum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    ulong hsum_lo1 = vec4n_horizontal_sum(dp_lo1);
    const ulong hsum_hi1 = vec4n_horizontal_sum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    for (; i < len; i++)
    {
        ulong v2 = vec2[i];
        hsum_lo0 += vec10[i] * v2;
        hsum_lo1 += vec11[i] * v2;
    }

    NMOD_RED(res[0], pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(res[1], pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}
#else  // PML_HAVE_AVX2
{
    ulong dp_lo0 = 0;
    ulong dp_hi0 = 0;
    ulong dp_lo1 = 0;
    ulong dp_hi1 = 0;

    slong i = 0;
    for ( ; i+7 < len; i += 8)
    {
        ulong v20, v21, v22, v23;

        v20 = vec2[i+0];
        v21 = vec2[i+1];
        v22 = vec2[i+2];
        v23 = vec2[i+3];
        dp_lo0 += vec10[i+0] * v20;
        dp_lo1 += vec11[i+0] * v20;
        dp_lo0 += vec10[i+1] * v21;
        dp_lo1 += vec11[i+1] * v21;
        dp_lo0 += vec10[i+2] * v22;
        dp_lo1 += vec11[i+2] * v22;
        dp_lo0 += vec10[i+3] * v23;
        dp_lo1 += vec11[i+3] * v23;

        v20 = vec2[i+4];
        v21 = vec2[i+5];
        v22 = vec2[i+6];
        v23 = vec2[i+7];
        dp_lo0 += vec10[i+4] * v20;
        dp_lo1 += vec11[i+4] * v20;
        dp_lo0 += vec10[i+5] * v21;
        dp_lo1 += vec11[i+5] * v21;
        dp_lo0 += vec10[i+6] * v22;
        dp_lo1 += vec11[i+6] * v22;
        dp_lo0 += vec10[i+7] * v23;
        dp_lo1 += vec11[i+7] * v23;

        dp_hi0 += dp_lo0 >> DOT_SPLIT_BITS;
        dp_lo0 &= DOT_SPLIT_MASK;
    }

    for ( ; i < len; i++)
    {
        ulong v2 = vec2[i];
        dp_lo0 += vec10[i] * v2;
        dp_lo1 += vec11[i] * v2;
    }

    NMOD_RED(res[0], pow2_precomp * dp_hi0 + dp_lo0, mod);
    NMOD_RED(res[1], pow2_precomp * dp_hi1 + dp_lo1, mod);
}
#endif  // PML_HAVE_AVX2
#endif  // FLINT_BITS == 64
