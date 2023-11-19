#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a CRT with num_primes of FLINT_BITS-8 size          */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_CRT_CRT(ulong num_primes)
{
    flint_rand_t state;
    fmpz_t comb;
    fmpz_multimod_CRT_t mCRT; 
    mp_ptr primes, residues;
    ulong i;

    flint_randinit(state);
    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    fmpz_init(comb);

    nmod_vec_primes(primes, num_primes, FLINT_BITS-8);
    for (i = 0; i < num_primes; i++)
    	residues[i] = n_randlimb(state) % primes[i];

    fmpz_multimod_CRT_init(mCRT, primes, num_primes);
    fmpz_multimod_CRT_CRT(comb, residues, mCRT);

    for (i = 0; i < num_primes; i++)
        assert (fmpz_fdiv_ui(comb, primes[i]) == residues[i]);
    
    fmpz_multimod_CRT_clear(mCRT);

    fmpz_clear(comb);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    ulong i;
    for (i = 1; i < 1000; i += 1)
    {
	check_fmpz_multimod_CRT_CRT(i);
    }

    return 0;
}