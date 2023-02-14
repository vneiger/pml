#include <flint/nmod_mat.h> // for nmod_mat_randtest()
#include "nmod_mat_extra.h" // for nmod_mat_rand()
#include "nmod_mat_poly.h"

void nmod_mat_poly_randtest(nmod_mat_poly_t matp,
                            flint_rand_t state, 
                            slong len)
{
    nmod_mat_poly_fit_length(matp, len);
    _nmod_mat_poly_set_length(matp, len);
    for (slong i = 0; i < len; ++i)
        nmod_mat_randtest(matp->coeffs + i, state);
    _nmod_mat_poly_normalise(matp);
}

void nmod_mat_poly_rand(nmod_mat_poly_t matp,
                        flint_rand_t state, 
                        slong len)
{
    nmod_mat_poly_fit_length(matp, len);
    _nmod_mat_poly_set_length(matp, len);
    for (slong i = 0; i < len; ++i)
        nmod_mat_rand(matp->coeffs + i, state);
    _nmod_mat_poly_normalise(matp);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
