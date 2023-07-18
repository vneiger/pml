#include "nmod_poly_extra.h" // for is_monic
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include <flint/nmod_poly.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - POPOV                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

int nmod_poly_mat_is_popov_rowwise(const nmod_poly_mat_t mat,
                                   const slong * shift)
{
    // matrix must be ordered weak Popov
    if (!nmod_poly_mat_is_ordered_weak_popov_rowwise(mat, shift))
        return 0;

    // retrieve pivot profile
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    slong * pivdeg = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_pivot_profile_rowwise(pivind, pivdeg, mat, shift);

    // pivots must be monic
    for (slong i = 0; i < mat->r; i++)
        if (!nmod_poly_is_monic(nmod_poly_mat_entry(mat, i, pivind[i])))
            return 0;

    // pivots must have degree strictly larger than other entries in same column
    for (slong i = 0; i < mat->r; i++)
        for (slong ii = 0; ii < mat->r; ii++)
            if (i != ii && nmod_poly_degree(nmod_poly_mat_entry(mat, ii, pivind[i])) >= pivdeg[i])
                return 0;

    flint_free(pivind);
    flint_free(pivdeg);
    return 1;
}

int nmod_poly_mat_is_popov_columnwise(const nmod_poly_mat_t mat,
                                      const slong * shift)
{
    // matrix must be ordered weak Popov
    if (!nmod_poly_mat_is_ordered_weak_popov_columnwise(mat, shift))
        return 0;

    // retrieve pivot profile
    slong * pivind = flint_malloc(mat->c * sizeof(slong));
    slong * pivdeg = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_pivot_profile_columnwise(pivind, pivdeg, mat, shift);

    // pivots must be monic
    for (slong j = 0; j < mat->c; j++)
        if (!nmod_poly_is_monic(nmod_poly_mat_entry(mat, pivind[j], j)))
            return 0;

    // pivots must have degree strictly larger than other entries in same column
    for (slong j = 0; j < mat->c; j++)
        for (slong jj = 0; jj < mat->c; jj++)
            if (j != jj && nmod_poly_degree(nmod_poly_mat_entry(mat, pivind[j], jj)) >= pivdeg[j])
                return 0;

    flint_free(pivind);
    flint_free(pivdeg);
    return 1;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - HERMITE                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

//int is_hermite(const nmod_poly_mat_t mat, orientation_t row_wise)
//{
//    slong rdim = mat->r, cdim = mat->c;
//    slong deg_mat = nmod_poly_mat_degree(mat);
//
//    if (row_wise)
//    {
//        slong shifts[cdim];
//        for (slong i = 0; i < cdim; i++)
//            shifts[i] = (cdim - i) * (deg_mat + 1);
//        return is_popov(mat, shifts, row_wise, 0);
//    }
//
//    slong shifts[rdim];
//    for (slong i = 0; i < rdim; i++)
//        shifts[i] = i * (deg_mat + 1);
//    return is_popov(mat, shifts, row_wise, 0);
//
//}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
