#include <stdlib.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/*                                                              */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_clear(fmpz_multimod_CRT_t mmod)
{
    ulong i;

    fmpz_clear(mmod->product_primes);
    fmpz_multi_mod_clear(mmod->top_mod);
    fmpz_multi_CRT_clear(mmod->top_CRT);
    
    for (i = 0; i < mmod->num_leaves; i++)
    {
        fmpz_clear(mmod->products_leaves + i);
        fmpz_multimod_naive_clear(mmod->leaves_mod[i]);
        fmpz_CRT_naive_clear(mmod->leaves_CRT[i]);
    }

    free(mmod->leaves_mod);
    free(mmod->leaves_CRT);
    free(mmod->products_leaves);
    _nmod_vec_clear(mmod->inverse_cofactors);
}
