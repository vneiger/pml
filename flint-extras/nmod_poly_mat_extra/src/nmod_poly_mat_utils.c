#include <flint/nmod_poly.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_mat_poly.h"


// rotate rows of mat from i to j (requirement: 0 <= i <= j < mat->r)
// and apply the corresponding transformation to vec (requirement: j < len(vec))
// If i == j, then nothing happens.
// vec can be NULL, in case it is omitted
// More precisely this performs simultaneously:
//      mat[i,:]     <--    mat[j,:]
//      mat[i+1,:]   <--    mat[i,:]
//      mat[i+2,:]   <--    mat[i+1,:]
//         ...       <--       ...
//      mat[j-1,:]   <--    mat[j-2,:]
//      mat[j,:]     <--    mat[j-1,:]
// as well as
//      vec[i]     <--    vec[j]
//      vec[i+1]   <--    vec[i]
//      vec[i+2]   <--    vec[i+1]
//        ...      <--      ...
//      vec[j-1]   <--    vec[j-2]
//      vec[j]     <--    vec[j-1]
void _nmod_poly_mat_rotate_rows(nmod_poly_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

        nmod_poly_struct * tmp_mat = mat->rows[j];
        for (slong ii = j; ii > i; ii--)
            mat->rows[ii] = mat->rows[ii-1];
        mat->rows[i] = tmp_mat;
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM CONSTANT                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

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
