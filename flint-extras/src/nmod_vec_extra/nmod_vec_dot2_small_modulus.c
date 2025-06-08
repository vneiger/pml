#include "machine_vectors.h"
#include "nmod_extra.h"
#include "nmod_vec_extra.h"

#ifdef PML_HAVE_MACHINE_VECTORS

/*------------------------------------------------------------*/
/** dot product for moduli less than 2^30                     */
/** reduction works if (p-1)^3*len < 2^96                     */
/** res[0] = dot(a1, b), res[1] = dot(a2, b)                  */
/** power_two = 2^45 mod p, p2 = (p,p), pinv2 = (1/p,1/p)     */
/*------------------------------------------------------------*/
void _nmod_vec_dot2_small_modulus(nn_ptr res, nn_ptr a1, nn_ptr a2, nn_ptr b, ulong len,
                                  ulong power_two, vec2d p2, vec2d pinv2)
{
    const ulong num_full_blocks = len >> 5;
    // mask for 45 low bits
    const vec4n low_bits = vec4n_set_n((UWORD(1) << 45) - 1);

    vec4n * a1_4n = (vec4n *) a1;
    vec4n * a2_4n = (vec4n *) a2;
    vec4n * b_4n = (vec4n *) b;
    vec4n sum_low1 = vec4n_zero();
    vec4n sum_high1 = vec4n_zero();
    vec4n sum_low2 = vec4n_zero();
    vec4n sum_high2 = vec4n_zero();

    // all blocks of length 32
    // each pass in the loop reads 8x4 = 32 coefficients
    for(ulong k = 0; k < num_full_blocks; k++)
    {
        vec4n val_b = vec4n_load_unaligned((nn_ptr) b_4n++);

        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
        val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));

        sum_high1 = vec4n_add(sum_high1, vec4n_bit_shift_right(sum_low1, 45));
        sum_low1 = vec4n_bit_and(sum_low1, low_bits);
        sum_high2 = vec4n_add(sum_high2, vec4n_bit_shift_right(sum_low2, 45));
        sum_low2 = vec4n_bit_and(sum_low2, low_bits);
    }

    // left with at most 31 coefficients
    // k is now the number of coefficients already read
    ulong k = num_full_blocks << 5;
    for (; k + 3 < len; k += 4)
    {
        vec4n val_b = vec4n_load_unaligned((nn_ptr) b_4n++);
        sum_low1 = vec4n_add(sum_low1, vec4n_mul(vec4n_load_unaligned((nn_ptr) a1_4n++), val_b));
        sum_low2 = vec4n_add(sum_low2, vec4n_mul(vec4n_load_unaligned((nn_ptr) a2_4n++), val_b));
    }

    sum_high1 = vec4n_add(sum_high1, vec4n_bit_shift_right(sum_low1, 45));
    sum_low1 = vec4n_bit_and(sum_low1, low_bits);
    sum_high2 = vec4n_add(sum_high2, vec4n_bit_shift_right(sum_low2, 45));
    sum_low2 = vec4n_bit_and(sum_low2, low_bits);

    // left with at most 3 coefficients
    nn_ptr as1 = (nn_ptr) a1_4n;
    nn_ptr as2 = (nn_ptr) a2_4n;
    nn_ptr bs = (nn_ptr) b_4n;
    ulong acc_last1 = UWORD(0);
    ulong acc_last2 = UWORD(0);
    for (; k < len; k++)
    {
        ulong val_b = (*bs++);
        acc_last1 += (*as1++) * val_b;
        acc_last2 += (*as2++) * val_b;
    }

    const ulong total_low1 = sum_low1[0] + sum_low1[1] + sum_low1[2] + sum_low1[3] + (acc_last1 & ((1L << 45) - 1));
    const ulong total_high1 = sum_high1[0] + sum_high1[1] + sum_high1[2] + sum_high1[3] + (acc_last1 >> 45);
    const ulong total_low2 = sum_low2[0] + sum_low2[1] + sum_low2[2] + sum_low2[3] + (acc_last2 & ((1L << 45) - 1));
    const ulong total_high2 = sum_high2[0] + sum_high2[1] + sum_high2[2] + sum_high2[3] + (acc_last2 >> 45);

    /* modified to avoid avx2-only functions and rely on avx2+neon from FLINT */
    vec1d sum0, sum1;
    sum0 = total_low1 + power_two*total_high1;
    sum1 = total_low2 + power_two*total_high2;
    sum0 = vec1d_reduce_to_0n(sum0, p2[0], pinv2[0]);
    sum1 = vec1d_reduce_to_0n(sum1, p2[1], pinv2[1]);
    /* vec2d redsum = vec2d_reduce_to_0n(sum, p2, pinv2); */
    /* res[0] = (ulong) redsum[0]; */
    /* res[1] = (ulong) redsum[1]; */
    res[0] = (ulong) sum0;
    res[1] = (ulong) sum1;
}

#endif  /* PML_HAVE_MACHINE_VECTORS */
