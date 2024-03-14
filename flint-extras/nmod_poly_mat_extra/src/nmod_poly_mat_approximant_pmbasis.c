#include "nmod_poly_mat_multiply.h"  // for middle product
#include "nmod_poly_mat_approximant.h"

void nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           slong order)
{
    if (order <= PMBASIS_THRES)
    {
        nmod_poly_mat_mbasis(appbas, shift, pmat, order);
        return;
    }

    const long order1 = order>>1;
    const long order2 = order - order1;
    nmod_poly_mat_t appbas2, residual;

    nmod_poly_mat_init(appbas2, pmat->r, pmat->r, pmat->modulus);
    nmod_poly_mat_init(residual, pmat->r, pmat->c, pmat->modulus);

    nmod_poly_mat_pmbasis(appbas, shift, pmat, order1);

    nmod_poly_mat_middle_product_naive(residual, appbas, pmat, order1, order2-1);

    nmod_poly_mat_pmbasis(appbas2, shift, residual, order2);

    nmod_poly_mat_mul(appbas, appbas2, appbas);

    nmod_poly_mat_clear(appbas2);
    nmod_poly_mat_clear(residual);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
