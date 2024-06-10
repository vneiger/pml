#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"

/*--------------------------------------------------------------*/
/* computes a random combination with num_primes terms          */
/* check against the result of the naive algorithm              */
/*--------------------------------------------------------------*/
void check_fmpz_CRT_naive_combine(ulong num_primes, ulong n_bits)
{
    flint_rand_t state;
    fmpz_t comb, check, temp;
    fmpz_CRT_naive_t mCRT; 
    nn_ptr primes, residues;
    ulong i, j;

    flint_rand_init(state);
    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    fmpz_init(comb);
    fmpz_init(check);
    fmpz_init(temp);

    nmod_vec_primes(primes, num_primes, n_bits);
    for (i = 0; i < num_primes; i++)
    	residues[i] = n_randlimb(state) % primes[i];

    fmpz_CRT_naive_init(mCRT, primes, num_primes);
    fmpz_CRT_naive_combine(comb, residues, mCRT);

    fmpz_set_ui(check, 0);
    for (i = 0; i < num_primes; i++)
    {
    	fmpz_set_ui(temp, residues[i]);
    	for (j = 0; j < num_primes; j++)
    	    if (i != j)
    		fmpz_mul_ui(temp, temp, primes[j]);
    	fmpz_add(check, check, temp);
    }
    assert(fmpz_equal(comb, check));

    fmpz_CRT_naive_clear(mCRT);

    fmpz_clear(temp);
    fmpz_clear(check);
    fmpz_clear(comb);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    ulong i;
    for (i = 1; i < 1000; i += 1)
    {
	check_fmpz_CRT_naive_combine(i, 60);
	check_fmpz_CRT_naive_combine(i, 50);
	check_fmpz_CRT_naive_combine(i, 29);
    }

    return 0;
}
