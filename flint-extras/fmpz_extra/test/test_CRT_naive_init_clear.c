#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a CRT                                    */
/* uses num_primes of FLINT_BITS-8 bits                         */
/*--------------------------------------------------------------*/
void check_fmpz_CRT_naive_init_clear(slong num_primes)
{
    fmpz_CRT_naive_t mCRT; 
    mp_ptr primes;

    primes = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, FLINT_BITS-8);
    fmpz_CRT_naive_init(mCRT, primes, num_primes);
    _nmod_vec_clear(primes);
    fmpz_CRT_naive_clear(mCRT);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    slong i;
    for (i = 1; i < 1000; i += 10)
	check_fmpz_CRT_naive_init_clear(i);

    return 0;
}
