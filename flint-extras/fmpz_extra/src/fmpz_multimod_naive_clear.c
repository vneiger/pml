#include <flint/flint.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* clears all memory used by mmod                               */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_clear(fmpz_multimod_naive_t mmod)
{
    ulong i;

    flint_free(mmod->primes);
    flint_free(mmod->mod);
    
    for (i = 0; i < mmod->num_primes; i++)
    {
        flint_free(mmod->powers_of_two[i]);
    }
    
    flint_free(mmod->powers_of_two);
    fmpz_clear(mmod->prod);
}
