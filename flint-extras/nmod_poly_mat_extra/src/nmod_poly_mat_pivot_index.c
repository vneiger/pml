#include "nmod_poly_mat_forms.h"

void nmod_poly_mat_pivot_index_rowwise(slong *pivind, const nmod_poly_mat_t mat)
{
    slong max, piv, d;
    for (slong i = 0; i < mat->r; i++)
    {
        max = -1;
        piv = -1;
        for (slong j = 0; j < mat->c; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            if (0 <= d && max <= d)
            {
                max = d;
                piv = j;
            }
        }
        pivind[i] = piv;
    }
}

void nmod_poly_mat_pivot_index_columnwise(slong *pivind, const nmod_poly_mat_t mat)
{
    slong max, piv, d;
    for (slong j = 0; j < mat->c; j++)
    {
        max = -1;
        piv = -1;
        for (slong i = 0; i < mat->r; i++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            if (0 <= d && max <= d)
            {
                max = d;
                piv = i;
            }
        }
        pivind[j] = piv;
    }
}

void nmod_poly_mat_pivot_index_shifted_rowwise(slong *pivind, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;
    slong min_shift = (mat->c > 0) ? shift[0] : 0;

    // find minimum of shift
    for (slong j = 0; j < mat->c; ++j)
        if (shift[j] < min_shift)
            min_shift = shift[j];

    for (slong i = 0; i < mat->r; i++)
    {
        max = min_shift-1; // zero rows will have this as rdeg
        piv = -1;
        for (slong j = 0; j < mat->c; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            if (0 <= d && max <= d+shift[j])
            {
                max = d;
                piv = j;
            }
        }
        pivind[i] = piv;
    }
}

void nmod_poly_mat_pivot_index_shifted_columnwise(slong *pivind, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;
    slong min_shift = (mat->r > 0) ? shift[0] : 0;

    // find minimum of shift
    for (slong i = 0; i < mat->r; ++i)
        if (shift[i] < min_shift)
            min_shift = shift[i];

    for (slong j = 0; j < mat->c; j++)
    {
        max = min_shift-1;
        piv = -1;
        for (slong i = 0; i < mat->r; i++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            if (0 <= d && max <= d+shift[i])
            {
                max = d;
                piv = i;
            }
        }
        pivind[j] = piv;
    }
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
