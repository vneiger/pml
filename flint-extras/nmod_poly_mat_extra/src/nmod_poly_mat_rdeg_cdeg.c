#include <flint/nmod_poly.h>

#include "nmod_poly_mat_forms.h"

void nmod_poly_mat_row_degree(slong *rdeg, const nmod_poly_mat_t mat, const slong *shift)
{
    if (shift == NULL)
    {
        slong max, d;
        for (slong i = 0; i < mat->r; i++)
        {
            max = -1;
            for (slong j = 0; j < mat->c; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
                if (max < d)
                    max = d;
            }
            rdeg[i] = max;
        }
    }
    else
    {
        slong max, d;
        slong min_shift = (mat->c > 0) ? shift[0] : 0;
        
        // find minimum of shift
        for (slong j = 0; j < mat->c; ++j)
            if (shift[j] < min_shift)
                min_shift = shift[j];
        
        for (slong i = 0; i < mat->r; i++)
        {
            max = min_shift-1; // zero rows will have this as rdeg
            for (slong j = 0; j < mat->c; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shift[j];
                // if new maximum at a nonzero entry, update
                if (shift[j] <= d && max < d)
                    max = d;
            }
            rdeg[i] = max;
        }
    }
}

void nmod_poly_mat_column_degree(slong *cdeg, const nmod_poly_mat_t mat, const slong *shift)
{
    if (shift == NULL)
    {
        slong max, d;
        for (slong j = 0; j < mat->c; j++)
        {
            max = -1;
            for (slong i = 0; i < mat->r; i++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
                if (max < d)
                    max = d;
            }
            cdeg[j] = max;
        }
    }
    else
    {
        slong max, d;
        slong min_shift = (mat->r > 0) ? shift[0] : 0;

        // find minimum of shift
        for (slong i = 0; i < mat->r; ++i)
            if (shift[i] < min_shift)
                min_shift = shift[i];

        for (slong j = 0; j < mat->c; j++)
        {
            max = min_shift-1;
            for (slong i = 0; i < mat->r; i++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shift[i];
                // if new maximum at a nonzero entry, update
                if (shift[i] <= d && max < d)
                    max = d;
            }
            cdeg[j] = max;
        }
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
