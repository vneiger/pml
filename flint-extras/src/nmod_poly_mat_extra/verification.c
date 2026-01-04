/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h"
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

    /* check appbas is square with the right dimension */
    if (appbas->r != rdim || appbas->c != rdim)
    {
        printf("basis has wrong row dimension or column dimension\n");
        success = 0;
    }

    /* check appbas is shifted reduced */
    if (!nmod_poly_mat_is_ordered_weak_popov(appbas, shift, orient))
    {
        printf("basis is not shifted-weak Popov\n");
        success = 0;
    }

    /* compute residual, check rows of appbas are approximants */
    nmod_poly_mat_mul(residual, appbas, pmat);

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            nmod_poly_set_trunc(pol, nmod_poly_mat_entry(residual, i, j), order);
            if (!nmod_poly_is_zero(pol))
            {
                printf("entry %ld, %ld of residual has valuation less than the order\n",i,j);
                success = 0;
            }
        }
    }

    /* check generation: follows ideas from Algorithm 1 in Giorgi-Neiger, ISSAC 2018 */

    /* generation, test 1: check determinant of appbas is lambda * x**D      *
     * since ordered weak Popov, deg(det(appbas)) is sum of diagonal degrees *
     */
    slong D = 0;
    for (slong i = 0; i < rdim; i++)
        D += nmod_poly_degree(nmod_poly_mat_entry(appbas, i, i));
    nmod_poly_mat_det(pol, appbas);
    if (nmod_poly_degree(pol) != D)
    {
        printf("determinant is not lambda * x**(sum(diag-deg))");
        success = 0;
    }

    /* generation, test 2: check that [P(0)  C] has full rank *
     * where C = (appbas * pmat * X^{-order})  mod X          *
     * (coefficient "C" of degree order of the residual)      *
     **/
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

/* TODO currently specialized to ROW_LOWER (or at least ROW_stuff) */
/* -> checks whether ker is a shift-ordered weak Popov left kernel basis of pmat */
int nmod_poly_mat_is_kernel(const nmod_poly_mat_t ker,
                            slong nz,
                            const slong * shift,
                            const nmod_poly_mat_t pmat,
                            orientation_t orient)
{
    const slong rdim = pmat->r;
    const slong cdim = pmat->c;
    const ulong prime = pmat->modulus;

    nmod_poly_mat_t residual;
    nmod_poly_mat_init(residual, nz, cdim, prime);

    nmod_poly_mat_t kernz;
    nmod_poly_mat_window_init(kernz, ker, 0, 0, nz, rdim);

    int success = 1;

    /* check kernel has the right column dimension */
    if (ker->c != rdim)
    {
        printf("basis has wrong column dimension\n");
        success = 0;
    }

    if (ker->r < nz)
    {
        printf("basis has wrong row dimension\n");
        success = 0;
    }

    /* check kernel is shifted reduced */
    if (!nmod_poly_mat_is_ordered_weak_popov(kernz, shift, orient))
    {
        /* printf("basis is not shifted-weak Popov\n"); */
        /* TODO temporarily allow reduced but not s-weak Popov */
        if (!nmod_poly_mat_is_reduced(kernz, shift, orient))
        {
            printf("basis is not shifted-reduced\n");
            success = 0;
        }
    }

    /* compute residual, check rows of ker are in the kernel */
    /* TODO enhancement: offer randomized option, using Freivalds-like randomized strategy (see ntl-extras) */
    nmod_poly_mat_mul(residual, kernz, pmat);
    /* nmod_poly_mat_degree_matrix_print_pretty(residual); */
    if (!nmod_poly_mat_is_zero(residual))
    {
        printf("not all rows are in the kernel\n");
        success = 0;
    }

    /* check rank */
    slong rk = nmod_poly_mat_rank(pmat);
    if (nz != rdim - rk)
    {
        printf("number of rows does not equal nullity\n");
        success = 0;
    }

    /* !! TODO check generation !! */

    nmod_poly_mat_clear(residual);

    return success;
}
