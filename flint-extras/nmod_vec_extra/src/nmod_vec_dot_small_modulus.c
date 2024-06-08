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

ulong _nmod_vec_dot_small_modulus_v2(nn_ptr a, nn_ptr b, ulong len,
                                      ulong power_two, double p, double pinv)
{
    // the dot product without modular reduction is
    //       dp  =  dp_lo + 2**sp_nb * dp_hi
    // where sp_nb is a fixed number of bits where we split (below, this is 45)
    const vec4n low_bits = vec4n_set_n((1L << 45) - 1);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    // NOTES:
    // -> constraint 1:
    // the above representation of dp requires len * (p-1)**2 <= 2**sp_nb * (2**64-1)
    // since p is < 2**32, a sufficient condition is len * 2**64 <= 2**(sp_nb + 63),
    // i.e. len <= 2**(sp_nb-1)
    // -> constraint 2:
    // having dp_lo and dp_hi, we will actually compute dp_lo + power_two * dp_hi
    // so we need power_two * dp_hi

    // each pass in the loop handles 8x4 = 32 coefficients
    ulong k = 0;
    for (; k+31 < len; k += 32)
    {
        // no overflow in next 8 lines if 8*(p-1)**2 < 2**64, i.e. p <= 2**(30.5)    (*)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 0), vec4n_load_unaligned(b+k+ 0)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 4), vec4n_load_unaligned(b+k+ 4)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 8), vec4n_load_unaligned(b+k+ 8)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+12), vec4n_load_unaligned(b+k+12)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+16), vec4n_load_unaligned(b+k+16)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+20), vec4n_load_unaligned(b+k+20)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+24), vec4n_load_unaligned(b+k+24)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+28), vec4n_load_unaligned(b+k+28)));

        // add 19 bits 45...63 to dp_hi
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right_45(dp_lo));
        // keep only 45 lower bits in dp_lo
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    // left with len-k < 32 coefficients,
    // this handles (len-k)//4 < 8 blocks of 4 coefficients, no overflow in this loop if (*) satisfied
    for (; k + 3 < len; k += 4)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k), vec4n_load_unaligned(b+k)));

    dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right_45(dp_lo));
    dp_lo = vec4n_bit_and(dp_lo, low_bits);

    // left with < 4 coefficients
    // no overflow in dp_last if 3*(p-1)**2 < 2**64, i.e. p <= 2**32 / sqrt(3)  (up to ~31.207 bits)
    ulong dp_last = 0;
    for (; k < len; k++)
        dp_last += a[k] * b[k];

    const ulong total_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3] + (dp_last & ((1L << 45) - 1));
    const ulong total_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (dp_last >> 45);

    double total = total_lo + power_two * total_hi;
    double dp = fma(-rint(total*pinv), p, total);
    //double dp = total - rint(total*pinv) * p;
    return (dp >= 0) ? dp : (dp+p);
}

ulong _nmod_vec_dot_small_modulus_v3(nn_ptr a, nn_ptr b, ulong len, ulong power_two, nmod_t mod)
{
    // the dot product without modular reduction is
    //       dp  =  dp_lo + 2**sp_nb * dp_hi
    // where sp_nb is a fixed number of bits where we split (below, this is 45)
    const vec4n low_bits = vec4n_set_n((1L << 45) - 1);
    vec4n dp_lo = vec4n_zero();
    vec4n dp_hi = vec4n_zero();

    // NOTES:
    // -> constraint 1:
    // the above representation of dp requires len * (p-1)**2 <= 2**sp_nb * (2**64-1)
    // since p is < 2**32, a sufficient condition is len * 2**64 <= 2**(sp_nb + 63),
    // i.e. len <= 2**(sp_nb-1)
    // -> constraint 2:
    // having dp_lo and dp_hi, we will actually compute dp_lo + power_two * dp_hi
    // so we need power_two * dp_hi

    // each pass in the loop handles 8x4 = 32 coefficients
    ulong k = 0;
    for (; k+31 < len; k += 32)
    {
        // no overflow in next 8 lines if 8*(p-1)**2 < 2**64, i.e. p <= 2**(30.5)    (*)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 0), vec4n_load_unaligned(b+k+ 0)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 4), vec4n_load_unaligned(b+k+ 4)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+ 8), vec4n_load_unaligned(b+k+ 8)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+12), vec4n_load_unaligned(b+k+12)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+16), vec4n_load_unaligned(b+k+16)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+20), vec4n_load_unaligned(b+k+20)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+24), vec4n_load_unaligned(b+k+24)));
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k+28), vec4n_load_unaligned(b+k+28)));

        // add 19 bits 45...63 to dp_hi
        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right_45(dp_lo));
        // keep only 45 lower bits in dp_lo
        dp_lo = vec4n_bit_and(dp_lo, low_bits);
    }

    // left with len-k < 32 coefficients,
    // this handles (len-k)//4 < 8 blocks of 4 coefficients, no overflow in this loop if (*) satisfied
    for (; k + 3 < len; k += 4)
        dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(a+k), vec4n_load_unaligned(b+k)));

    dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right_45(dp_lo));
    dp_lo = vec4n_bit_and(dp_lo, low_bits);

    // left with < 4 coefficients
    // no overflow in dp_last if 3*(p-1)**2 < 2**64, i.e. p <= 2**32 / sqrt(3)  (up to ~31.207 bits)
    ulong dp_last = 0;
    for (; k < len; k++)
        dp_last += a[k] * b[k];


    const ulong total_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3] + (dp_last & ((1L << 45) - 1));
    const ulong total_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (dp_last >> 45);
    const ulong dp = power_two * total_hi + total_lo;
    ulong res;
    NMOD_RED(res, dp, mod);
    //return n_mod2_preinv(dp, mod.n, mod.ninv);
    return res;
}



/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
