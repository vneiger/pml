#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* uses 55 bits primes                                          */
/* uses num_primes, total bit length is max_bit_length          */
/* reduces 100 random integer modulo all primes                 */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_CRT_reduce(ulong num_primes)
{
    flint_rand_t state;
    fmpz_multimod_CRT_t mmod; 
    mp_ptr primes, residues;
    ulong i, j;

    flint_randinit(state);

    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, 55);
    fmpz_multimod_CRT_init(mmod, primes, num_primes);
    
    for (j = 0; j < 100; j++)
    {
	fmpz_t A;
	fmpz_init(A);
	fmpz_randtest(A, state, 54 * num_primes);
        fmpz_multimod_CRT_reduce(residues, A, mmod);
	for (i = 0; i < num_primes; i++)
	    if (residues[i] != fmpz_fdiv_ui(A, primes[i]))
            {
                printf("A=");
                fmpz_print(A);
                printf("\n");
                printf("error with i=%lu, num_primes=%lu\n", i, num_primes);
                printf("%lu %lu\n", residues[i], fmpz_fdiv_ui(A, primes[i]));
                printf("p=%lu\n", primes[i]);
                exit(-1);
            }
	fmpz_clear(A);
    }

   fmpz_multimod_CRT_clear(mmod);
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
    for (i = 1; i < 10000; i += 10)
    {
	check_fmpz_multimod_CRT_reduce(i);
    }
    
    return 0;
}
