#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* prepares the vectors of coefficients and inverses            */
/* ------------------------------------------------------------ */
void fmpz_CRT_naive_init(fmpz_CRT_naive_t mCRT, mp_srcptr primes, ulong num_primes)
{
    unsigned long i, k;
    int j;

    mCRT->num_primes = num_primes;
    mCRT->primes = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t));
    mCRT->inverses = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t));
    mCRT->mod = (nmod_t *) flint_malloc(num_primes * sizeof(nmod_t));
    mCRT->prime_bit_length = _nmod_vec_max_bits(primes, num_primes); /* finds the bit-length of the primes */
    fmpz_init(mCRT->prod);
    fmpz_set_ui(mCRT->prod, 1);

    for (i = 0; i < num_primes; i++)
    {
        nmod_init(&mCRT->mod[i], primes[i]);
	fmpz_mul_ui(mCRT->prod, mCRT->prod, primes[i]);
    }

    mCRT->num_limbs = (fmpz_bits(mCRT->prod) + FLINT_BITS - 1) / FLINT_BITS;   /* number of limbs in the product */
    mCRT->coefficients = (mp_ptr *) flint_malloc(mCRT->num_limbs * sizeof(mp_ptr));
    for (k = 0; k < mCRT->num_limbs; k++)
    {
	mCRT->coefficients[k] = (mp_ptr) flint_malloc(num_primes * sizeof(mp_limb_t));
    }

    for (i = 0; i < num_primes; i++)
    {
	fmpz_t cofactor;
	fmpz_init(cofactor);
	fmpz_divexact_ui(cofactor, mCRT->prod, primes[i]);

	mCRT->inverses[i] = n_invmod(fmpz_fdiv_ui(cofactor, primes[i]), primes[i]);

	if (!COEFF_IS_MPZ(*cofactor))
	{
	    fmpz value = *cofactor;
	    mCRT->coefficients[0][i] = value;
	    for (k = 1; k < mCRT->num_limbs; k++)
	    {
		mCRT->coefficients[k][i] = 0;
	    }
	}
	else
	{
	    __mpz_struct *ptr = COEFF_TO_PTR(*cofactor);
	    mp_ptr coeffs = ptr->_mp_d;
	    for (j = 0; j < ptr->_mp_size; j++)
	    {
		mCRT->coefficients[j][i] = coeffs[j];
	    }
	    for (k = ptr->_mp_size; k < mCRT->num_limbs; k++)
	    {
		mCRT->coefficients[k][i] = 0;
	    }
	}
	fmpz_clear(cofactor);
    }

}
