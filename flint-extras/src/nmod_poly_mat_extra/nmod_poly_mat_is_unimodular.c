/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

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
    ulong pt1 = 0;
    ulong pt2 = 0;
    while (pt1 == pt2)
    {
        pt1 = n_randint(state, pmat->modulus);
        pt2 = n_randint(state, pmat->modulus);
    }

    // det1 = det(pmat(pt1))
    nmod_mat_t ev1;
    nmod_mat_init(ev1, pmat->r, pmat->c, pmat->modulus);
    nmod_poly_mat_evaluate_nmod(ev1, pmat, pt1);
    ulong det1 = nmod_mat_det(ev1);

    // det2 = det(pmat(pt2))
    nmod_mat_t ev2;
    nmod_mat_init(ev2, pmat->r, pmat->c, pmat->modulus);
    nmod_poly_mat_evaluate_nmod(ev2, pmat, pt2);
    ulong det2 = nmod_mat_det(ev1);

    nmod_mat_clear(ev1);
    nmod_mat_clear(ev2);

    return (det1 != 0 && det1 == det2);
}
