#include <assert.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>

#include "nmod_vec_extra.h"
#include "fmpz_extra.h"

/*--------------------------------------------------------------*/
/* computes a random combination with num_primes terms          */
/* check against the result of the naive algorithm              */
/*--------------------------------------------------------------*/
int check_fmpz_CRT_naive_combine(ulong num_primes, ulong n_bits, flint_rand_t state)
{
    fmpz_t comb, check, temp;
    fmpz_CRT_naive_t mCRT; 
    nn_ptr primes, residues;
    ulong i, j;

    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    fmpz_init(comb);
    fmpz_init(check);
    fmpz_init(temp);

    nmod_vec_primes(primes, num_primes, n_bits);
    for (i = 0; i < num_primes; i++)
    	residues[i] = n_randint(state, primes[i]);

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
    int res = fmpz_equal(comb, check);

    fmpz_CRT_naive_clear(mCRT);

    fmpz_clear(temp);
    fmpz_clear(check);
    fmpz_clear(comb);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);

    return res;
}

TEST_FUNCTION_START(fmpz_CRT_naive_combine, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong num_primes = 1;  /* corner case for i == 0 */
        if (i == 1)
            num_primes = 1000;  /* corner case for i == 1 */
        if (i > 1)
            num_primes = 1 + n_randint(state, 200);

        result = check_fmpz_CRT_naive_combine(num_primes, 60, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 60\n", num_primes);
        result = check_fmpz_CRT_naive_combine(num_primes, 50, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 60\n", num_primes);
        result = check_fmpz_CRT_naive_combine(num_primes, 29, state);
    }

    TEST_FUNCTION_END(state);
}
