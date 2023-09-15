#include <flint/nmod_poly.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_mat_poly.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SETTING AND GETTING COEFFICIENTS                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, slong degree)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            nmod_mat_set_entry(res, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j), degree));
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TRANSPOSE, TRUNC, SHIFT, REVERSE                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_set_trunc(nmod_poly_mat_t tmat, const nmod_poly_mat_t pmat, long len)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_set_trunc(tmat->rows[i] + j, pmat->rows[i] + j, len);
}

void nmod_poly_mat_truncate(nmod_poly_mat_t pmat, long len)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_truncate(pmat->rows[i] + j, len);
}

void nmod_poly_mat_shift_left(nmod_poly_mat_t smat, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < smat->r; i++)
        for (slong j = 0; j < smat->c; j++)
            nmod_poly_shift_left(nmod_poly_mat_entry(smat, i, j), nmod_poly_mat_entry(pmat, i, j), k);
}

void nmod_poly_mat_shift_right(nmod_poly_mat_t smat, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < smat->r; i++)
        for (slong j = 0; j < smat->c; j++)
            nmod_poly_shift_right(smat->rows[i] + j, pmat->rows[i] + j, k);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM CONSTANT                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_set_from_nmod_mat(nmod_poly_mat_t pmat, const nmod_mat_t cmat)
{
    for (slong i = 0; i < cmat->r; ++i)
        for (slong j = 0; j < cmat->c; ++j)
        {
            if (nmod_mat_entry(cmat, i, j) == 0)
                nmod_poly_zero(nmod_poly_mat_entry(pmat, i, j));
            else
            {
                nmod_poly_realloc(nmod_poly_mat_entry(pmat, i, j), 1);
                nmod_poly_mat_entry(pmat, i, j)->coeffs[0]
                                    = nmod_mat_entry(cmat, i, j);
            }
        }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM MATRIX POLYNOMIAL                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_set_from_mat_poly0(nmod_poly_mat_t pmat, const nmod_mat_poly0_t matp)
{
	nmod_poly_mat_zero(pmat);
	for (slong k = 0; k <= matp->degree; k++)
        for (slong i = 0; i < matp->r; ++i)
            for (slong j = 0; j < matp->c; ++j)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k, nmod_mat_get_entry(matp->mat + k, i, j));
}

void nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                           const nmod_mat_poly_t matp,
                                           slong order)
{
    if (order > matp->length)
        order = matp->length;

    // prepare memory
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_fit_length(nmod_poly_mat_entry(pmat, i, j), order);

    // fill data
    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < pmat->r; i++)
            for (slong j = 0; j < pmat->c; j++)
                nmod_poly_mat_entry(pmat, i, j)->coeffs[k] = nmod_mat_poly_entry(matp, k, i, j);

    // normalize
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), order);
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, i, j));
        }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
