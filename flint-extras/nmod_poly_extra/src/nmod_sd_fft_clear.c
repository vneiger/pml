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

    if (F->order > 0)
        flint_free(F->powers_inv_w_over_2);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
