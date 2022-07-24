#include "nmod_poly_mat_forms.h"

void nmod_poly_mat_pivot_profile_rowwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat)
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
        pivdeg[i] = max;
        pivind[i] = piv;
    }
}

void nmod_poly_mat_pivot_profile_columnwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat)
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
        pivdeg[j] = max;
        pivind[j] = piv;
    }
}

void nmod_poly_mat_pivot_profile_shifted_rowwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;

    for (slong i = 0; i < mat->r; i++)
    {
        piv = -1; // zero rows will have this as pivind
        // the next starting value for max plays no role: in the next for loop,
        // either piv==-1 until the end or max has been set correctly at first
        // nonzero entry in the row
        max = 0;
        for (slong j = 0; j < mat->c; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            // if encountering first nonzero entry in the row, update as well
            if (0 <= d && (piv == -1 || max <= d+shift[j]))
            {
                piv = j;
                max = d+shift[piv];
            }
        }
        pivind[i] = piv;
        pivdeg[i] = (piv == -1) ? -1 : (max - shift[piv]);
    }
}

void nmod_poly_mat_pivot_profile_shifted_columnwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;

    for (slong j = 0; j < mat->c; j++)
    {
        piv = -1; // zero columns will have this as pivind
        // the next starting value for max plays no role: in the next for loop,
        // either piv==-1 until the end or max has been set correctly at first
        // nonzero entry in the column
        max = 0;
        for (slong i = 0; i < mat->r; i++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            // if encountering first nonzero entry in the column, update as well
            if (0 <= d && (piv == -1 || max <= d+shift[i]))
            {
                piv = i;
                max = d+shift[piv];
            }
        }
        pivind[j] = piv;
        pivdeg[j] = (piv == -1) ? -1 : (max - shift[piv]);
    }
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
