#include "nmod_poly_extra.h"

#ifdef HAS_INT128

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_64_fft_clear(nmod_64_fft_t F)
{
    ulong i;
    
    for (i = 0; i <= F->order; i++)
    {
        flint_free(F->powers_w[i]);
        flint_free(F->i_powers_w[i]);
        flint_free(F->powers_inv_w[i]);
        flint_free(F->i_powers_inv_w[i]);
    }
    for (i = 0; i < F->order; i++)
    {
        flint_free(F->powers_inv_w_over_2[i]);
        flint_free(F->i_powers_inv_w_over_2[i]);
    }

    flint_free(F->powers_inv_w_over_2);
    flint_free(F->i_powers_inv_w_over_2);
    flint_free(F->powers_w);
    flint_free(F->powers_inv_w);
    flint_free(F->powers_inv_2);
    flint_free(F->i_powers_inv_2);
}

#endif
