#include <flint/nmod_poly.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* RANDOM - UNIFORM                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_rand(nmod_poly_t poly,
                    flint_rand_t state,
                    slong len)
{
    nmod_poly_fit_length(poly, len);
    _nmod_vec_rand(poly->coeffs, state, len, poly->mod);
    poly->length = len;
    _nmod_poly_normalise(poly);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
