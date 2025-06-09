#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* clears all memory                                            */
/* ------------------------------------------------------------ */
void 
fmpz_CRT_naive_clear(fmpz_CRT_naive_t mCRT)
{
    ulong i;

    flint_free(mCRT->primes);
    flint_free(mCRT->inverses);
    for (i = 0; i < mCRT->num_limbs; i++)
	flint_free(mCRT->coefficients[i]);

    flint_free(mCRT->coefficients);
    flint_free(mCRT->mod);
    fmpz_clear(mCRT->prod);
}
