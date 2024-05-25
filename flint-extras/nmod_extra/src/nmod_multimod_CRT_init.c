#include <flint/fmpz.h>

#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* initializes all data in C                                  */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_init(nmod_multimod_CRT_t C, ulong modulus, ulong num_primes)
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
        ulong p0p1_red2 = nmod_mul(PRIME0, PRIME1, mod2);
        C->invp0p1_p2 = (double) nmod_inv(p0p1_red2, mod2);

        nmod_init(&mod3, PRIME3);
        C->p0p1_red3 = nmod_mul(PRIME0, PRIME1, mod3);
        ulong p0p1p2_red3 = nmod_mul(C->p0p1_red3, PRIME2, mod3);
        C->invp0p1p2_p3 = (double) nmod_inv(p0p1p2_red3, mod3);
    }
    else // large modulus. this is inspired by multimod and CRT in fft_small
    {
        ulong i, len;
        fmpz_t prod;
        nn_ptr coeffs;
        __mpz_struct *ptr;
       
        if (num_primes == 1)
            return;
        
        fmpz_init(prod);
        fmpz_set_ui(prod, 1);

        for (i = 0; i < num_primes; i++)
            fmpz_mul_ui(prod, prod, C->mod_primes[i].n);

        len = num_primes * num_primes + num_primes + num_primes;
        C->data = FLINT_ARRAY_ALLOC(len, ulong);
        for (i = 0; i < len; i++)
            C->data[i] = 0;

        ptr = COEFF_TO_PTR(*prod);
	coeffs = ptr->_mp_d;
        
        flint_mpn_copyi(C->data + num_primes * num_primes, coeffs, num_primes);
        for (i = 0; i < num_primes; i++)
        {
            mpn_divexact_1(C->data + i*num_primes, coeffs, num_primes, C->mod_primes[i].n);
            C->data[num_primes*num_primes + num_primes + i] = mpn_mod_1(C->data + i*num_primes, num_primes, C->mod_primes[i].n);
            C->inverse_cofactors[i] = nmod_inv(C->data[num_primes*num_primes + num_primes + i], C->mod_primes[i]);
        }

        fmpz_clear(prod);
    }
    
}
