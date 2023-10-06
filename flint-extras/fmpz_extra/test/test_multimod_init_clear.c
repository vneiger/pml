#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* uses num_primes of FLINT_BITS-8 bits                         */
/* max bit-size of operands is max_bit_length                   */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_init_clear(slong num_primes, slong max_bit_length)
{
    fmpz_multimod_t mmod; 
    mp_ptr primes;

    primes = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, FLINT_BITS-8);
    fmpz_multimod_init(mmod, primes, num_primes, max_bit_length);
    _nmod_vec_clear(primes);
    fmpz_multimod_clear(mmod);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    slong i;
    for (i = 1; i < 1000; i += 50)
    {
	check_fmpz_multimod_init_clear(i, 59);
	check_fmpz_multimod_init_clear(i, 10000);
    }
    return 0;
}
