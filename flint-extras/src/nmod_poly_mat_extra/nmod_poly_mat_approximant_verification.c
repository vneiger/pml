#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

/* TODO currently specialized to ROW_LOWER (or at least ROW_stuff) */
int nmod_poly_mat_is_approximant_basis(const nmod_poly_mat_t appbas,
                                       const nmod_poly_mat_t pmat,
                                       slong order,
                                       const slong * shift,
                                       orientation_t orient)
{
    // context
    const slong rdim = pmat->r;
    const slong cdim = pmat->c;
    const ulong prime = pmat->modulus;

    nmod_poly_mat_t residual;
    nmod_poly_t pol;
    nmod_mat_t CP0;

    nmod_poly_init(pol, prime);
    nmod_poly_mat_init(residual, rdim, cdim, prime);
    nmod_mat_init(CP0, rdim, rdim+cdim, prime);

    int success = 1;

    // check appbas is square with the right dimension
    if (appbas->r != rdim || appbas->c != rdim)
    {
        printf("basis has wrong row dimension or column dimension\n");
        success = 0;
    }

    // check appbas has form at least "form"
    if (!nmod_poly_mat_is_ordered_weak_popov(appbas, shift, orient))
    {
        printf("basis is not shifted-reduced\n");
        success = 0;
    }

    // compute residual, check rows of appbas are approximants
    nmod_poly_mat_mul(residual, appbas, pmat);

    for(slong i = 0; i < rdim; i++)
    {
        for(slong j = 0; j < cdim; j++)
        {
            nmod_poly_set_trunc(pol, nmod_poly_mat_entry(residual, i, j), order);
            if (!nmod_poly_is_zero(pol))
            {
                printf("entry %ld, %ld of residual has valuation less than the order\n",i,j);
                success = 0;
            }
        }
    }

    // check generation: follows ideas from Algorithm 1 in Giorgi-Neiger, ISSAC 2018

    // generation, test 1: check determinant of appbas is lambda * x**D
    // since ordered weak Popov, deg-det is sum of diagonal degrees
    slong D = 0;
    for (slong i = 0; i < rdim; i++)
        D += nmod_poly_degree(nmod_poly_mat_entry(appbas, i, i));
    nmod_poly_mat_det(pol, appbas);
    if (nmod_poly_degree(pol) != D)
    {
        printf("determinant is not lambda * x**(sum(diag-deg))");
        success = 0;
    }

    // generation, test 2: check that [P(0)  C] has full rank
    // where C = (appbas * pmat * X^{-order})  mod X
    // (coefficient "C" of degree order of the residual)
    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < rdim; j++)
        {
            ulong c = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(appbas, i, j), 0);
            nmod_mat_set_entry(CP0, i, j, c);
        }
        for (slong j = 0; j < cdim; j++)
        {
            ulong c = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(residual, i, j), order);
            nmod_mat_set_entry(CP0, i, rdim+j, c);
        }
    }

    slong rank = nmod_mat_rank(CP0);
    if (rank != rdim)
    {
        printf("generation test (step: [C P(0)] has full rank) failed");
        success = 0;
    }

    nmod_poly_clear(pol);
    nmod_poly_mat_clear(residual);
    nmod_mat_clear(CP0);

    return success;
}
