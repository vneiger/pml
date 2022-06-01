#include <flint/nmod_poly.h>

#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_plain_fft_clear(nmod_plain_fft_t F)
{
    slong i;
    
    for (i = 0; i <= F->order; i++)
    {
        flint_free(F->powers_w[i]);
        flint_free(F->powers_inv_w[i]);
    }
    flint_free(F->powers_w);
    flint_free(F->powers_inv_w);
    flint_free(F->powers_inv_2);
}
