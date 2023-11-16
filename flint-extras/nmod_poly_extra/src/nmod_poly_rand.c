#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"
#include <flint/nmod_poly.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* RANDOM - UNIFORM                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_rand(nmod_poly_t pol,
                    flint_rand_t state,
                    slong len)
{
    if (len <= 0)
        nmod_poly_zero(pol);
    else
    {
        nmod_poly_fit_length(pol, len);
        _nmod_vec_rand(pol->coeffs, state, len, pol->mod);
        pol->length = len;
        _nmod_poly_normalise(pol);
    }
}

void nmod_poly_rand_monic(nmod_poly_t pol,
                          flint_rand_t state,
                          slong len)
{
    if (len <= 0)
        nmod_poly_zero(pol);
    else
    {
        nmod_poly_fit_length(pol, len);
        _nmod_vec_rand(pol->coeffs, state, len - 1, pol->mod);
        pol->coeffs[len-1] = 1;
        pol->length = len;
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
