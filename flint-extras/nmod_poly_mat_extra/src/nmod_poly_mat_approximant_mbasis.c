#include "nmod_mat_extra.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

void nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                          slong * shift,
                          const nmod_poly_mat_t pmat,
                          ulong order)
{
    nmod_mat_poly_t app, matp;
    nmod_mat_poly_init(matp, pmat->r, pmat->c, pmat->modulus);
    // TODO improve: set init
    nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, order);
    nmod_mat_poly_init(app, pmat->r, pmat->r, pmat->modulus);
    nmod_mat_poly_mbasis(app, shift, matp, order);
    // TODO improve: set init
    nmod_poly_mat_set_from_mat_poly(appbas, app);
    nmod_mat_poly_clear(matp);
    nmod_mat_poly_clear(app);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
