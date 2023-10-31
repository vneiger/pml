#include <flint/nmod_poly.h>

#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_sd_fft_clear(nmod_sd_fft_t F)
{
    ulong i;
    
    for (i = 0; i <= F->order; i++)
    {
        flint_free(F->powers_w[i]);
        flint_free(F->powers_inv_w_t[i]);
    }
    
    for (i = 0; i < F->order; i++)
    {
        flint_free(F->powers_inv_w_over_2[i]);
    }

    flint_free(F->powers_w);
    flint_free(F->powers_inv_w_t);
    flint_free(F->powers_inv_2);
    flint_free(F->powers_inv_w_over_2);
}
