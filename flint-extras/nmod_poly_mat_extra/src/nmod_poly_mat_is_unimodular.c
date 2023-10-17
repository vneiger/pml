#include "nmod_poly_mat_utils.h"
#include <flint/flint.h>
#include <flint/nmod_poly_mat.h>

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

int nmod_poly_mat_is_unimodular_randomized(const nmod_poly_mat_t pmat, flint_rand_t state)
{
    if (pmat->r != pmat->c)
        return 0;

    if (pmat->modulus == 1)
        return 0;

    // pick two random distinct points
    mp_limb_t pt1 = 0;
    mp_limb_t pt2 = 0;
    while (pt1 == pt2)
    {
        pt1 = n_randint(state, pmat->modulus);
        pt2 = n_randint(state, pmat->modulus);
    }

    // det1 = det(pmat(pt1))
    nmod_mat_t ev1;
    nmod_mat_init(ev1, pmat->r, pmat->c, pmat->modulus);
    nmod_poly_mat_evaluate_nmod(ev1, pmat, pt1);
    mp_limb_t det1 = nmod_mat_det(ev1);

    // det2 = det(pmat(pt2))
    nmod_mat_t ev2;
    nmod_mat_init(ev2, pmat->r, pmat->c, pmat->modulus);
    nmod_poly_mat_evaluate_nmod(ev2, pmat, pt2);
    mp_limb_t det2 = nmod_mat_det(ev1);

    nmod_mat_clear(ev1);
    nmod_mat_clear(ev2);

    return (det1 != 0 && det1 == det2);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
