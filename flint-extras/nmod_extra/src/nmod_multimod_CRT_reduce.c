#include <flint/machine_vectors.h>
#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* residues[j][i] = input[i] mod prime[j]                     */
/* for i < nb, j < num_primes                                 */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_reduce(mp_ptr *residues, mp_ptr input, ulong nb, nmod_multimod_CRT_t C)
{
    ulong i, j;
    
    // if p < pi (and in particular if p < 2^49), nothing to do
    // else, use the 2^32 trick from nmod_poly_mul
    // TODO: if pi < p < 2^50, then since pi >= 2^49, we have p < 2pi, so we could use reduce_2n_to_n
    for (j = 0; j < C->num_primes; j++)
    {
        mp_ptr res = residues[j];
        if (C->p < C->primes[j])
            for (i = 0; i < nb; i++)
                    res[i] = input[i];
        else
        {
            vec1d n, ninv;
            vec4d n4, ninv4;
            n = C->primes[j];
            ninv = C->primes_inv[j];
            n4 = vec4d_set_d(n);
            ninv4 = vec4d_set_d(ninv);
            
            i = 0;
            if (nb >= 4)
            {
                for (; i + 3 < nb; i += 4)
                    {
                        vec4n t = vec4n_load_unaligned(input + i);
                        vec4d tlo = vec4n_convert_limited_vec4d(vec4n_bit_and(t, vec4n_set_n(4294967295)));
                        vec4d thi = vec4n_convert_limited_vec4d(vec4n_bit_shift_right(t, 32));
                        vec4d_store_unaligned_mp_ptr(res + i, vec4d_addmod(tlo,
                                                                           vec4d_reduce_pm1no_to_0n(
                                                                               vec4d_mulmod(thi, vec4d_set_d(1L << 32), n4, ninv4),
                                                                               n4),
                                                                           n4));
                    }
                }
            
            for (; i < nb; i++)
            {
                vec1n t = input[i];
                vec1d tlo = (vec1d) (t & 4294967295);
                vec1d thi = (vec1d) (t >> 32);
                res[i] = vec1d_addmod(tlo,
                                      vec1d_reduce_pm1no_to_0n(
                                          vec1d_mulmod(thi, (vec1d)(1L << 32), n, ninv),
                                          n),
                                          n);
            }
        }
    }
}
