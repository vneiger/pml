#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* copies primes in mmod-> primes (length is num_primes)        */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_init(fmpz_multimod_naive_t mmod, mp_srcptr primes, ulong num_primes)
{
    ulong i, u;
    mp_limb_t max_prime;
    
    mmod->num_primes = num_primes;
    mmod->primes = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t));
    mmod->mod = (nmod_t *) flint_malloc(num_primes * sizeof(nmod_t));
    mmod->powers_of_two = (mp_ptr *) flint_malloc(num_primes * sizeof(mp_ptr));

    fmpz_init(mmod->prod);
    fmpz_set_ui(mmod->prod, 1);

    max_prime = 0;
    for (i = 0; i < num_primes; i++)
    {
        mmod->primes[i] = primes[i];
        nmod_init(&mmod->mod[i], primes[i]);
    	fmpz_mul_ui(mmod->prod, mmod->prod, primes[i]);
        if (primes[i] > max_prime)
            max_prime = primes[i];
    }

    mmod->num_limbs = fmpz_size(mmod->prod);

    for (i = 0; i < num_primes; i++)
    {
        mp_limb_t two_FLINT_BITS;
        nmod_t mod;
        mp_ptr v;

        mmod->small_moduli = 1;
        
        mod = mmod->mod[i];
        two_FLINT_BITS = ( UWORD(1) << (FLINT_BITS-1) ) % primes[i];
        two_FLINT_BITS = nmod_add(two_FLINT_BITS, two_FLINT_BITS, mod);

        // small primes:
        // every 64-bit word in the input will be sliced into 3x24 bits 
        // NB: could do better, we don't use all 3x24 = 72 bits.
        // instead, could slice 3x64 bits into 8x24 (instead of 9x24 here)
        if (max_prime < (1L << 30))
        {
            ulong len;
            mp_limb_t two_24, two_48;

            two_24 = (1L << 24) % primes[i];
            two_48 = (1L << 48) % primes[i];
            
            len = 3 * mmod->num_limbs; 
            len = ((len + 3) >> 2) << 2; // must be a multiple of 4
            mmod->powers_of_two[i] = (mp_ptr) aligned_alloc(32, len * sizeof(mp_limb_t));
            v = mmod->powers_of_two[i];

            /* Powers of 2^FLINT_BITS, 2^24*2^FLINT_BITS, 2^48*2^FLINT_BITS modulo p */
            v[0] = 1;
            v[1] = two_24;
            v[2] = two_48;
            for (u = 3; u < len; u += 3)
            {
                v[u] = nmod_mul(v[u - 3], two_FLINT_BITS, mod);
                if (u + 1 < len)
                {
                    v[u + 1] = nmod_mul(v[u], two_24, mod);
                    if (u + 2 < len)
                        v[u + 2] = nmod_mul(v[u + 1], two_24, mod);
                }
            }
        }
        else{
            mmod->small_moduli = 0;
            
            mmod->powers_of_two[i] = (mp_ptr) flint_malloc(mmod->num_limbs * sizeof(mp_limb_t));
            /* Powers of 2^FLINT_BITS modulo p */
            v = mmod->powers_of_two[i];
            v[0] = 1;
            for (u = 1; u < mmod->num_limbs; u++)
                v[u] = nmod_mul(v[u - 1], two_FLINT_BITS, mod);
        }
    }
}
/* #include <flint/flint.h> */
/* #include <flint/fmpz.h> */

/* #include "fmpz_extra.h" */

/* /\* ------------------------------------------------------------ *\/ */
/* /\* copies primes in mmod-> primes (length is num_primes)        *\/ */
/* /\* ------------------------------------------------------------ *\/ */
/* void fmpz_multimod_naive_init(fmpz_multimod_naive_t mmod, mp_srcptr primes, ulong num_primes) */
/* { */
/*     ulong i, u; */

/*     mmod->num_primes = num_primes; */
/*     mmod->primes = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t)); */
/*     mmod->mod = (nmod_t *) flint_malloc(num_primes * sizeof(nmod_t)); */
/*     mmod->powers_of_two = (mp_ptr *) flint_malloc(num_primes * sizeof(mp_ptr)); */

/*     fmpz_init(mmod->prod); */
/*     fmpz_set_ui(mmod->prod, 1); */

/*     for (i = 0; i < num_primes; i++) */
/*     { */
/*         mmod->primes[i] = primes[i]; */
/*         nmod_init(&mmod->mod[i], primes[i]); */
/*     	fmpz_mul_ui(mmod->prod, mmod->prod, primes[i]); */
/*     } */

/*     mmod->num_limbs = fmpz_size(mmod->prod); */
/*     for (i = 0; i < num_primes; i++) */
/*     { */
/*         mp_limb_t two_FLINT_BITS; */
/*         nmod_t mod; */
/*         mp_ptr v; */
        
/*         mmod->powers_of_two[i] = (mp_ptr) flint_malloc(mmod->num_limbs * sizeof(mp_limb_t)); */
/*         mod = mmod->mod[i]; */
/*         two_FLINT_BITS = ( UWORD(1) << (FLINT_BITS-1) ) % primes[i]; */
/*         two_FLINT_BITS = nmod_add(two_FLINT_BITS, two_FLINT_BITS, mod); */
        
/*         /\* Powers of 2^FLINT_BITS modulo p *\/ */
/*         v = mmod->powers_of_two[i]; */
/*         v[0] = 1; */
/*         for (u = 1; u < mmod->num_limbs; u++) */
/*             v[u] = nmod_mul(v[u - 1], two_FLINT_BITS, mod); */
/*     } */
/* } */
