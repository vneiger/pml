#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* compute A s.t. A mod primes[j] = m[j] for all j              */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_CRT(fmpz_t A, mp_srcptr m, const fmpz_multimod_CRT_t mmod)
{
    ulong i, cpt;
    fmpz * tmp;
    mp_ptr precomp_m;

    precomp_m = _nmod_vec_init(mmod->num_primes);
    tmp = (fmpz *) malloc(mmod->num_leaves * sizeof(fmpz));

    cpt = 0;
    for (i = 0; i < mmod->num_leaves; i++)
    {
        ulong j, nb;
        mp_ptr inv_cof;
        
        nb = mmod->leaves_CRT[i]->num_primes;
        inv_cof = mmod->inverse_cofactors;
        for (j = 0; j < nb; j++)
        {
            precomp_m[cpt] = nmod_mul(inv_cof[cpt], m[cpt], mmod->leaves_CRT[i]->mod[j]);
            cpt++;
        }

        fmpz_init(tmp + i);
        fmpz_CRT_naive_combine(tmp + i, precomp_m + mmod->size_leaves * i, mmod->leaves_CRT[i]);
    }
    
    if (mmod->num_leaves == 1)
        fmpz_set(A, tmp);
    else
        fmpz_multi_CRT_combine(A, mmod->top_CRT, tmp);
    
    for (i = 0; i < mmod->num_leaves; i++)
        fmpz_clear(tmp + i);
    free(tmp);
    _nmod_vec_clear(precomp_m);
}
