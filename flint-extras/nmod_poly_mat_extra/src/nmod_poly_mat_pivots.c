#include "nmod_poly_mat_forms.h"

/**********************************************************************
*                        shifted pivot index                         *
**********************************************************************/

void nmod_poly_mat_pivot_index_rowwise(slong *pivind, const nmod_poly_mat_t mat, const slong *shift)
{
    if (shift == NULL)
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
 
    else
    {
        slong max, piv, d;
        
        for (slong i = 0; i < mat->r; i++)
        {
            piv = -1;
            // the next starting value for max plays no role: in the next for loop,
            // either piv==-1 until the end or max has been set correctly at first
            // nonzero entry in the row
            max = 0;
            for (slong j = 0; j < mat->c; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
                // if new maximum (or equal maximum) reached at a nonzero entry, update
                // if encountering first nonzero entry in the row, update as well
                if (0 <= d && (piv==-1 || max <= d+shift[j]))
                {
                    max = d + shift[j];
                    piv = j;
                }
            }
            pivind[i] = piv;
        }
    }
}

void nmod_poly_mat_pivot_index_columnwise(slong *pivind, const nmod_poly_mat_t mat, const slong *shift)
{
    if (shift == NULL)
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
    else
    {
        slong max, piv, d;

        for (slong j = 0; j < mat->c; j++)
        {
            piv = -1;
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
                    max = d + shift[i];
                    piv = i;
                }
            }
            pivind[j] = piv;
        }
    }
}

/**********************************************************************
*                       shifted pivot profile                        *
**********************************************************************/
void nmod_poly_mat_pivot_profile_rowwise(slong *pivind,
                                         slong *pivdeg,
                                         const nmod_poly_mat_t mat,
                                         const slong *shift)
{
    if (shift == NULL)
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
    else
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
}

void nmod_poly_mat_pivot_profile_columnwise(slong *pivind,
                                            slong *pivdeg,
                                            const nmod_poly_mat_t mat,
                                            const slong *shift)
{
    if (shift == NULL)
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
    else
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
}

/**********************************************************************
*                        echelon pivot index                         *
**********************************************************************/

void nmod_poly_mat_lechelon_pivot_index_rowwise(slong * pivind, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong i = 0; i < mat->r; i++)
    {
        // rightmost nonzero entry in row i if nonzero, otherwise -1
        piv = mat->c - 1;
        while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
            piv--;
        pivind[i] = piv;
    }
}

void nmod_poly_mat_uechelon_pivot_index_rowwise(slong * pivind, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong i = 0; i < mat->r; i++)
    {
        // leftmost nonzero entry in row i if nonzero, otherwise mat->c
        piv = 0;
        while (piv < mat->c && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
            piv++;
        pivind[i] = piv;
    }
}

void nmod_poly_mat_lechelon_pivot_index_columnwise(slong * pivind, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong j = 0; j < mat->c; j++)
    {
        // topmost nonzero entry in column j if nonzero, otherwise mat->r
        piv = 0;
        while (piv < mat->r && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
            piv++;
        pivind[j] = piv;
    }
}

void nmod_poly_mat_uechelon_pivot_index_columnwise(slong * pivind, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong j = 0; j < mat->c; j++)
    {
        // bottommost nonzero entry in column j if nonzero, otherwise -1
        piv = mat->r - 1;
        while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
            piv--;
        pivind[j] = piv;
    }
}

/**********************************************************************
*                       echelon pivot profile                        *
**********************************************************************/

void nmod_poly_mat_lechelon_pivot_profile_rowwise(slong * pivind, slong * pivdeg, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong i = 0; i < mat->r; i++)
    {
        // rightmost nonzero entry in row i if nonzero, otherwise -1
        piv = mat->c - 1;
        while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
            piv--;
        pivind[i] = piv;
        pivdeg[i] = (piv==-1) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, i, piv));
    }
}

void nmod_poly_mat_uechelon_pivot_profile_rowwise(slong * pivind, slong * pivdeg, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong i = 0; i < mat->r; i++)
    {
        // leftmost nonzero entry in row i if nonzero, otherwise mat->c
        piv = 0;
        while (piv < mat->c && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
            piv++;
        pivind[i] = piv;
        pivdeg[i] = (piv==mat->c) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, i, piv));
    }
}

void nmod_poly_mat_lechelon_pivot_profile_columnwise(slong * pivind, slong * pivdeg, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong j = 0; j < mat->c; j++)
    {
        // topmost nonzero entry in column j if nonzero, otherwise mat->r
        piv = 0;
        while (piv < mat->r && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
            piv++;
        pivind[j] = piv;
        pivdeg[j] = (piv == mat->r) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, piv, j));
    }
}

void nmod_poly_mat_uechelon_pivot_profile_columnwise(slong * pivind, slong * pivdeg, const nmod_poly_mat_t mat)
{
    slong piv;
    for (slong j = 0; j < mat->c; j++)
    {
        // bottommost nonzero entry in column j if nonzero, otherwise -1
        piv = mat->r - 1;
        while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
            piv--;
        pivind[j] = piv;
        pivdeg[j] = (piv==-1) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, piv, j));
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
