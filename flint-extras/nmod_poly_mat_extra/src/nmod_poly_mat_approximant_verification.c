#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include <flint/nmod_poly.h>

int nmod_poly_mat_is_approximant_basis(const nmod_poly_mat_t appbas,
                                       const nmod_poly_mat_t pmat,
                                       slong order,
                                       const slong * shift,
                                       orientation_t orient)
{
    // context
    const slong rdim = pmat->r, cdim = pmat->c;
    const ulong prime = pmat->modulus;

    // check appbas is square with the right dimension
    if (appbas->r != rdim || appbas->c != rdim)
    {
        printf("basis has wrong row dimension or column dimension\n");
        return 0;
    }

    // check appbas has form at least "form"
    if (!nmod_poly_mat_is_ordered_weak_popov(appbas, shift, orient))
    {
        printf("basis is not shifted-reduced\n");
        return 0;
    }

    // compute residual
    nmod_poly_mat_t residual;
    nmod_poly_mat_init(residual, rdim, cdim, prime);
    // FIXME next line could be a multiplication truncated at order order+1
    nmod_poly_mat_mul(residual, appbas, pmat);

    // check rows of appbas are approximants
    nmod_poly_t pol;
    nmod_poly_init(pol, prime);
    for(slong i = 0; i < rdim; i++)
        for(slong j = 0; j < cdim; j++)
        {
            nmod_poly_set_trunc(pol, nmod_poly_mat_entry(residual, i, j), order);
            if (!nmod_poly_is_zero(pol))
            {
                printf("entry %ld, %ld of residual has valuation less than the order\n",i,j);
                return 0;
            }
        }
    nmod_poly_clear(pol);
    nmod_poly_mat_clear(residual);

    // TODO check generation!

    // 
    //slong lead_pos[rdim];
    //nmod_poly_mat_pivot_index_shifted_rowwise(lead_pos, appbas, shift);
    //printf("\nleading positions\n");
    //for (slong i = 0; i < rdim; i++)
    //    printf("%lu ", lead_pos[i]);

    return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
