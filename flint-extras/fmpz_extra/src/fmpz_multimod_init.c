#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* copies primes in mmod-> primes (length is num_primes)        */
/* max_bit_length = maximum bit-length of inputs we will reduce */
/* ------------------------------------------------------------ */

void fmpz_multimod_init(fmpz_multimod_t mmod, mp_srcptr primes, ulong num_primes, slong max_bit_length)
{
    ulong i, u;

    mmod->num_primes = num_primes;
    mmod->primes = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t));
    mmod->mod = (nmod_t *) flint_malloc(num_primes * sizeof(nmod_t));
    mmod->powers_of_two = (mp_ptr *) flint_malloc(num_primes * sizeof(mp_ptr));
    mmod->prime_bit_length = _nmod_vec_max_bits(primes, num_primes); /* finds the bit-length of the primes */
    mmod->num_limbs = (max_bit_length + FLINT_BITS - 1) / FLINT_BITS; /* round up to get the number of limbs  */

    /* printf("num_primes %ld\n", num_primes); */
    /* printf("mbl %ld\n", max_bit_length); */
    /* printf("num limbs %ld\n", mmod->num_limbs); */

    for (i = 0; i < num_primes; i++)
    {
        mp_limb_t two_FLINT_BITS;
        mp_ptr v;
        nmod_t mod;

        mmod->primes[i] = primes[i];
        nmod_init(&mmod->mod[i], primes[i]);
        
        mmod->powers_of_two[i] = (mp_ptr) flint_malloc(mmod->num_limbs * sizeof(mp_limb_t));
        mod = mmod->mod[i];
        two_FLINT_BITS = ( UWORD(1) << (FLINT_BITS-1) ) % primes[i];
        two_FLINT_BITS = nmod_add(two_FLINT_BITS, two_FLINT_BITS, mod);
        
        /* Powers of 2^FLINT_BITS modulo p */
        v = mmod->powers_of_two[i];
        v[0] = 1;
        for (u = 1; u < mmod->num_limbs; u++)
            v[u] = nmod_mul(v[u - 1], two_FLINT_BITS, mod);
                
    }
}
