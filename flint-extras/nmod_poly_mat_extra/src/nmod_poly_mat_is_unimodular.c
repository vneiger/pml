#include "nmod_poly_mat_utils.h"

int nmod_poly_mat_is_unimodular(const nmod_poly_mat_t pmat)
{
    if (pmat->r != pmat->c)
        return 0;
    nmod_poly_t det;
    nmod_poly_init(det, pmat->modulus);
    nmod_poly_mat_det(det, pmat);
    int is_uni = (det->length == 1);
    nmod_poly_clear(det);
    return is_uni;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
