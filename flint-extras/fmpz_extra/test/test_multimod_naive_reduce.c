#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* uses num_primes of size n_bits                               */
/* reduces 100 random integer modulo all primes                 */
/* input initially reduced modulo the product of primes         */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_naive_reduce(ulong num_primes, ulong n_bits)
{
    flint_rand_t state;
    fmpz_multimod_naive_t mmod; 
    nn_ptr primes, residues;
    ulong i, j;

    flint_rand_init(state);

    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, n_bits);
    fmpz_multimod_naive_init(mmod, primes, num_primes);
    
    for (j = 0; j < 100; j++)
    {
    	fmpz_t A;
    	fmpz_init(A);

    	fmpz_randtest(A, state, fmpz_bits(mmod->prod) - 1);
   	fmpz_multimod_naive_reduce(residues, A, mmod);
    	for (i = 0; i < num_primes; i++)
    	    if (residues[i] != fmpz_fdiv_ui(A, primes[i]))
            {
                printf("A=");
                fmpz_print(A);
                printf("\n");
                printf("error with i=%lu, num_primes=%lu and n_bits=%lu\n", i, num_primes, n_bits);
                printf("%lu %lu\n", residues[i], fmpz_fdiv_ui(A, primes[i]));
                printf("p=%lu\n", primes[i]);
                exit(-1);
            }
    	fmpz_clear(A);
    }

   fmpz_multimod_naive_clear(mmod);
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
	check_fmpz_multimod_naive_reduce(i, 60);
	check_fmpz_multimod_naive_reduce(i, 50);
	check_fmpz_multimod_naive_reduce(i, 29);
    }
    return 0;
}
