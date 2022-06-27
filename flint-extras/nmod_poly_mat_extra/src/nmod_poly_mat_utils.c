#include <flint/nmod_poly.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SETTING AND GETTING COEFFICIENTS                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, slong degree)
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

void nmod_poly_mat_shift_left(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < res->r; i++)
        for (slong j = 0; j < res->c; j++)
            if (!nmod_poly_is_zero(nmod_poly_mat_entry(res, i, j)))
                nmod_poly_shift_left(nmod_poly_mat_entry(res, i, j), nmod_poly_mat_entry(pmat, i, j), k);
}

void nmod_poly_mat_shift_right(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < res->r; i++)
        for (slong j = 0; j < res->c; j++)
            if (!nmod_poly_is_zero(nmod_poly_mat_entry(res, i, j)))
                nmod_poly_shift_right(nmod_poly_mat_entry(res, i, j), nmod_poly_mat_entry(pmat, i, j), k);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CONVERT                                                    */
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

void nmod_poly_mat_set_from_mat_poly(nmod_poly_mat_t pmat, const nmod_mat_poly_t matp)
{
    // TODO reinsert this kind of things, with flag "debug mode" or something?
	//if (pmat->modulus != matp->mod)
	//{
	//	printf("\nERROR! Wrong modulus: nmod_mat_poly_to_poly_mat\n");
	//	return;
	//}
	//if (pmat->r != r || pmat->c != c)
	//{
	//	printf("\nERROR! Wrong shape: nmod_mat_poly_to_poly_mat\n");
	//	printf("shape pmat = (%ld, %ld), shape matp = (%ld, %ld)\n", pmat->r, pmat->c, r, c);
	//	return;
	//}
    // TODO start with the right length,
    // then fill, then normalize
    // (zero is not necessary if all entries are touched, simply reallocate to
    // right length)
	nmod_poly_mat_zero(pmat);
	for (slong k = 0; k <= matp->degree; k++)
        for (slong i = 0; i < matp->r; ++i)
            for (slong j = 0; j < matp->c; ++j)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k, nmod_mat_get_entry(matp->mat + k, i, j));
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
