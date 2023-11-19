#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* initializes all data in C                                  */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_init(nmod_multimod_CRT_t C, mp_limb_t modulus, ulong num_primes)
{
    FLINT_ASSERT(num_primes <= 4);

    C->num_primes = num_primes;
    C->p = modulus;
    nmod_init(&C->mod, modulus);

    C->primes[0] = PRIME0;
    C->primes[1] = PRIME1;
    C->primes[2] = PRIME2;
    C->primes[3] = PRIME3;
    
    C->primes_inv[0] = 1 / (double)PRIME0;
    C->primes_inv[1] = 1 / (double)PRIME1;
    C->primes_inv[2] = 1 / (double)PRIME2;
    C->primes_inv[3] = 1 / (double)PRIME3;

    nmod_init(&C->mod_primes[0], PRIME0);
    nmod_init(&C->mod_primes[1], PRIME1);
    nmod_init(&C->mod_primes[2], PRIME2);
    nmod_init(&C->mod_primes[3], PRIME3);

// default value set to 4
    if (num_primes == 0)
        num_primes = 4;
    

    if (modulus < (1L << 50)) // small modulus: case use SIMD floating-point representation 
    {
        C->pinv = 1 / (double)modulus;

        C->p0_red = (double) (PRIME0 % C->p);
        C->p1_red = (double) (PRIME1 % C->p);
        C->p0p1_red = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(C->p0_red, C->p1_red, C->p, C->pinv), C->p);
        C->p0p1p2_red = vec1d_reduce_pm1no_to_0n(vec1d_mulmod(C->p0p1_red, (double) (PRIME2 % C->p), C->p, C->pinv), C->p);
        C->p0_red2 = (double) (PRIME0 % PRIME2);
        
        nmod_t mod1;
        nmod_t mod2;
        nmod_t mod3;
        
        nmod_init(&mod1, PRIME1);
        C->invp0_p1 = nmod_inv(PRIME0, mod1);

        nmod_init(&mod2, PRIME2);
        mp_limb_t p0p1_red2 = nmod_mul(PRIME0, PRIME1, mod2);
        C->invp0p1_p2 = (double) nmod_inv(p0p1_red2, mod2);

        nmod_init(&mod3, PRIME3);
        C->p0p1_red3 = nmod_mul(PRIME0, PRIME1, mod3);
        mp_limb_t p0p1p2_red3 = nmod_mul(C->p0p1_red3, PRIME2, mod3);
        C->invp0p1p2_p3 = (double) nmod_inv(p0p1p2_red3, mod3);
    }
    else // large modulus. this is inspired by multimod and CRT in fft_small
    {
        ulong i, old_len;
        mp_ptr old_data;
        
        if (num_primes == 1)
            return;
        
        old_data = FLINT_ARRAY_ALLOC(1*1 + 1 + 1, ulong);
        old_len = 1;
        
        for (i = 0; i < num_primes; i++)
        {
            mp_limb_t p;
            p = C->mod_primes[i].n;
            
            if (i == 0)
            {
                old_data[0] = 1; // coprime
                old_data[1] = p; // prod
                old_data[2] = 1; // coprime_red
            }
            else
            {
                ulong len;
                ulong *t,  *tt;
                
                len = old_len;
                
                t = FLINT_ARRAY_ALLOC(2*(len + 2), ulong);
                tt = t + (len + 2);
                
                t[len + 1] = 0;
                t[len] = mpn_mul_1(t, old_data + len*i, len, p);
                
                /* leave enough room for (product of primes)*(number of primes) */
                len += 2;
                mpn_mul_1(tt, t, len, i + 1);
                while (tt[len - 1] == 0)
                    len--;
                
                flint_free(old_data);
                
                /* set product of primes */
                old_data = FLINT_ARRAY_ALLOC((i+1)*len + len + (i+1), ulong);
                flint_mpn_copyi(old_data + (i+1)*len, t, len);
                
                flint_free(t);
                old_len = len;
            }
        }
        
        /* set cofactors */
        for (i = 0; i < num_primes; i++)
        {
            mpn_divexact_1(old_data + i*old_len, old_data + num_primes*old_len, old_len, C->mod_primes[i].n);
            old_data[num_primes*old_len + old_len + i] = mpn_mod_1(old_data + i*old_len, old_len, C->mod_primes[i].n);
            C->inverse_cofactors[i] = nmod_inv(old_data[num_primes*old_len + old_len + i], C->mod_primes[i]);
        }
        
        C->data = FLINT_ARRAY_ALLOC(num_primes*old_len + num_primes + old_len, ulong);
        flint_mpn_copyi(C->data, old_data, num_primes*old_len + num_primes + old_len);
        flint_free(old_data);
    }
}
