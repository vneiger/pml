#include "nmod_poly_mat_forms.h"

void column_degrees_shifted(slong *res, const nmod_poly_mat_t mat, const slong *shifts)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max, d;
    for (slong i = 0; i < cdim; i++)
    {
        max = -1;
        for (slong j = 0; j < rdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
            if (max < d)
                max = d;
        }
        res[i] = max;
    }
}

void row_degrees_shifted(slong *res, const nmod_poly_mat_t mat, const slong *shifts)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max, d;
    for (slong i = 0; i < rdim; i++)
    {
        max = -1;
        for (slong j = 0; j < cdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
            if (max < d)
                max = d;
        }
        res[i] = max;
    }
}

void degree_matrix(slong *res, const nmod_poly_mat_t mat, const slong *shifts,
                   orientation_t row_wise)
{

    slong rdim = mat->r, cdim = mat->c;
    slong d;

    if (row_wise)
    {
        for(slong i = 0; i < rdim; i++)
        {
            for(slong j = 0; j < cdim; j++)
            {

                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
                *(res + (i * rdim) + j) = d;
            }
        }
        return;
    }
    for(slong i = 0; i < cdim; i++)
        for(slong j = 0; j < rdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
            *(res + (i * cdim) + j) = d;
        }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
