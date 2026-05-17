/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_extra.h"  /* for NMOD_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_multiply.h"

void nmod_poly_mat_multiply(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2)
{
    slong len1 = nmod_poly_mat_max_length(pmat1);
    slong len2 = nmod_poly_mat_max_length(pmat2);

    /* TODO once available, call multiply with constant */
    /* TODO once available, call mat-vec | vec-mat multiply */

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (pmat1->r == 1 || pmat2->c == 1  /* vec-mat or mat-vec */
        || pmat1->c == 1                /* independent products */
        || len1 == 1 || len2 == 1)      /* constant lhs or rhs */
    {
        nmod_poly_mat_mul(res, pmat1, pmat2);
        return;
    }

    /* FIXME decide where to handle this: currently done at multiply levels... */
    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, pmat1->r, pmat2->c, pmat1->modulus);
        nmod_poly_mat_multiply(tmp, pmat1, pmat2);
        nmod_poly_mat_swap_entrywise(res, tmp);
        nmod_poly_mat_clear(tmp);
        return;
    }

    /* from here, all dimensions and lengths are >= 2 */
    /* TODO thresholds essentially tuned from square cases
     * -> would be safe to check that they are ok for fat vectors */
    /* TODO see impact of FLINT+BLAS on thresholds */

    const slong dim = n_cbrt(pmat1->r * pmat1->c * pmat2->c);
    const slong len = len1 + len2 - 1;
    const ulong modn = pmat1->modulus;

    if (dim > 12)
    {
        if (NMOD_POLY_CAN_USE_GEOMETRIC(modn, len) && len > 300)
            nmod_poly_mat_mul_geometric(res, pmat1, pmat2);
        else if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len))
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
        /* FIXME small fields (cannot use geom/vdm): should probably call waksman in some cases... */
    }

    else if (dim > 10)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && len < 200)
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
            nmod_poly_mat_mul_waksman(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
    }

    else if (dim > 8)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && len < 120)
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
            nmod_poly_mat_mul_waksman(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
    }

    else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
        nmod_poly_mat_mul_waksman(res, pmat1, pmat2);

    else
        nmod_poly_mat_mul(res, pmat1, pmat2);
}

