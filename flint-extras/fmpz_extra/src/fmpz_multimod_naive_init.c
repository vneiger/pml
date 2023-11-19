#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* copies primes in mmod-> primes (length is num_primes)        */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_init(fmpz_multimod_naive_t mmod, mp_srcptr primes, ulong num_primes)
{
    ulong i, u;

    mmod->num_primes = num_primes;
    mmod->primes = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t));
    mmod->mod = (nmod_t *) flint_malloc(num_primes * sizeof(nmod_t));
    mmod->powers_of_two = (mp_ptr *) flint_malloc(num_primes * sizeof(mp_ptr));

    fmpz_init(mmod->prod);
    fmpz_set_ui(mmod->prod, 1);

    for (i = 0; i < num_primes; i++)
    {
        mmod->primes[i] = primes[i];
        nmod_init(&mmod->mod[i], primes[i]);
    	fmpz_mul_ui(mmod->prod, mmod->prod, primes[i]);
    }

    mmod->num_limbs = fmpz_size(mmod->prod);
    for (i = 0; i < num_primes; i++)
    {
        mp_limb_t two_FLINT_BITS;
        nmod_t mod;
        mp_ptr v;
        
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
