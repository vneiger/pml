#include "nmod_poly_extra.h"

#ifdef HAS_AVX2

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_clear(nmod_avx2_32_fft_t F)
{
    slong i;
    
    for (i = 0; i <= F->order; i++)
    {
        free(F->powers_w[i]);
        free(F->i_powers_w[i]);
        free(F->powers_inv_w[i]);
        free(F->i_powers_inv_w[i]);
    }
    flint_free(F->powers_w);
    flint_free(F->powers_inv_w);
    flint_free(F->powers_inv_2);
    flint_free(F->i_powers_inv_2);
}

#endif
