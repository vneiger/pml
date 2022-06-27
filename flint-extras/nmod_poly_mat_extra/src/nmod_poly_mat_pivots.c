#include "nmod_poly_mat_forms.h"

void leading_matrix(nmod_mat_t res, const nmod_poly_mat_t mat,
                    const slong *shifts, orientation_t row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    nmod_poly_struct *P;
    if (row_wise)
    {
        slong rdeg[rdim];
        row_degrees_shifted(rdeg, mat, shifts);
        for(slong i = 0; i < rdim; i++)
            for(slong j = 0; j < cdim; j++)
            {
                P = nmod_poly_mat_entry(mat, i, j);
                nmod_mat_set_entry(res, i, j, nmod_poly_get_coeff_ui(P, rdeg[i] - shifts[j]));
            }
        return;
    }
    {
        slong cdeg[cdim];
        column_degrees_shifted(cdeg, mat, shifts);
        for(slong i = 0; i < cdim; i++)
            for(slong j = 0; j < rdim; j++)
            {
                P = nmod_poly_mat_entry(mat, j, i);
                nmod_mat_set_entry(res, j, i, nmod_poly_get_coeff_ui(P, cdeg[i] - shifts[j]));
            }
    }
}

void leading_positions(slong *res, const nmod_poly_mat_t mat,
                       const slong *shifts, orientation_t row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max;
    slong d;
    slong ind; // TODO initial value
    if (row_wise)
    {
        for (slong i = 0; i < rdim; i++)
        {
            max = 0;
            for (slong j = 0; j < cdim; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
                if (max <= d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
        return;
    }
    {
        for (slong i = 0; i < cdim; i++)
        {
            max = 0;
            for (slong j = 0; j < rdim; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
                if (max <= d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
    }
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
