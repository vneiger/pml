#include "nmod_poly_mat_forms.h"

int nmod_poly_mat_is_ordered_weak_popov_rowwise(const nmod_poly_mat_t mat,
                                                const slong * shift)
{
    slong pivind[mat->r];
    nmod_poly_mat_pivot_index_rowwise(pivind, mat, shift);

    // first row must be nonzero
    if (mat->r > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->r - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_ordered_weak_popov_columnwise(const nmod_poly_mat_t mat,
                                                   const slong * shift)
{
    slong pivind[mat->c];
    nmod_poly_mat_pivot_index_columnwise(pivind, mat, shift);

    // first column must be nonzero
    if (mat->c > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->c - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
