#include <flint/flint.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"
#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* compute A = sum_i m[i] \prod_{j \ne i} primes[j]             */
/* ------------------------------------------------------------ */
void 
fmpz_CRT_naive_combine(fmpz_t A, nn_srcptr m, const fmpz_CRT_naive_t mCRT)
{
    mp_size_t size;
    ulong tmp[3];
    nn_ptr dat;
    mpz_t temp;
    ulong i;

    mpz_init2(temp, FLINT_BITS * (mCRT->num_limbs + 2));
    dat = temp->_mp_d;

    for (i = 0; i < mCRT->num_limbs + 2; i++)
    {
	dat[i] = 0;
    }

    /* computes a 3-limb dot product */
    /* adds it to dat                */
    for (i = 0; i < mCRT->num_limbs; i++)
    {
	nmod_vec_integer_dot_product(tmp, mCRT->coefficients[i], m, mCRT->num_primes, FLINT_BITS, mCRT->prime_bit_length);
    	add_sssaaaaaa(dat[i+2], dat[i+1], dat[i], dat[i+2], dat[i+1], dat[i], tmp[2], tmp[1], tmp[0]);
    }

    /* final adjustment of the mpz */
    size = mCRT->num_limbs + 2;
    while (size >= 1 && dat[size-1] == 0)
    {
    	size--;
    }
    if (size == 1 && dat[0] == 0)
    {
    	size = 0;
    }

    temp->_mp_size = size;
    fmpz_set_mpz(A, temp); // TODO !! to be fixed !!

    mpz_clear(temp);
}
