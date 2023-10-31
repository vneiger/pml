#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* uses 62 bits primes                                          */
/* uses num_primes, total bit length is max_bit_length          */
/* reduces 100 random integer modulo all primes                   */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_reduce(ulong num_primes, ulong max_bit_length)
{
    flint_rand_t state;
    fmpz_multimod_t mmod; 
    mp_ptr primes, residues;
    ulong i, j;

    flint_randinit(state);

    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, FLINT_BITS-2);
    fmpz_multimod_init(mmod, primes, num_primes, max_bit_length);
    
    for (j = 0; j < 100; j++)
    {
    	fmpz_t A;
    	fmpz_init(A);
    	fmpz_randtest(A, state, max_bit_length);
    	fmpz_multimod_reduce(residues, A, mmod);
    	for (i = 0; i < num_primes; i++)
    	    assert (residues[i] == fmpz_fdiv_ui(A, primes[i]));
    	fmpz_clear(A);
    }

    fmpz_multimod_clear(mmod);
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
    for (i = 1; i < 1000; i += 100)
    {
	check_fmpz_multimod_reduce(i, 100);
	check_fmpz_multimod_reduce(i, 100000);
    }
    return 0;
}
