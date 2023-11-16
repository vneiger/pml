#include "nmod_poly_mat_forms.h"


/**********************************************************************
*                    shifted pivot profile, vector                    *
**********************************************************************/

void _nmod_poly_vec_pivot_profile(slong * pivind,
                                  slong * pivdeg,
                                  const nmod_poly_struct * vec,
                                  const slong * shift,
                                  slong len,
                                  orientation_t orient)
{
    if (orient == ROW_LOWER || orient == COL_UPPER)
    {
        if (shift == NULL)
        {
            slong max = -1;
            slong piv = -1;
            slong d;
            for (slong j = 0; j < len; j++)
            {
                d = nmod_poly_degree(vec+j);
                if (0 <= d && max <= d)
                {
                    max = d;
                    piv = j;
                }
            }
            *pivdeg = max;
            *pivind = piv;
        }
        else
        {
            slong piv = -1; // zero vec will have this as pivind
            slong max = 0; // the starting value for max plays no role: in the next
                        // for loop, either piv==-1 until the end or max has
                        // been set correctly at first nonzero entry in the vec
            slong d;
            for (slong j = 0; j < len; j++)
            {
                d = nmod_poly_degree(vec+j);
                // if new maximum (or equal maximum) reached at a nonzero entry, update
                // if encountering first nonzero entry in the row, update as well
                if (0 <= d && (piv == -1 || max <= d+shift[j]))
                {
                    piv = j;
                    max = d+shift[piv];
                }
            }
            *pivind = piv;
            *pivdeg = (piv == -1) ? -1 : (max - shift[piv]);
        }
    }
    else // orient == COL_LOWER || orient == ROW_UPPER
    {
        if (shift == NULL)
        {
            slong max = -1;
            slong piv = -1;
            slong d;
            for (slong j = 0; j < len; j++)
            {
                d = nmod_poly_degree(vec+j);
                if (0 <= d && max < d)
                {
                    max = d;
                    piv = j;
                }
            }
            *pivdeg = max;
            *pivind = (piv == -1) ? len : piv;
        }
        else
        {
            slong piv = -1;
            slong max = 0; // the starting value for max plays no role: in the next
                        // for loop, either piv==-1 until the end or max has
                        // been set correctly at first nonzero entry in vec
            slong d;
            for (slong j = 0; j < len; j++)
            {
                d = nmod_poly_degree(vec+j);
                // if new strict maximum reached at a nonzero entry, update
                // if encountering first nonzero entry in vec, update as well
                if (0 <= d && (piv == -1 || max < d+shift[j]))
                {
                    piv = j;
                    max = d+shift[piv];
                }
            }
            *pivind = piv;
            *pivdeg = (piv == -1) ? len : (max - shift[piv]);
        }
    }
}

/**********************************************************************
*                        shifted pivot index                         *
**********************************************************************/

void nmod_poly_mat_pivot_index(slong *pivind,
                               const nmod_poly_mat_t mat,
                               const slong * shift,
                               orientation_t orient)
{
    if (orient == ROW_LOWER || orient == ROW_UPPER)
    {
        slong buf;
        for (slong i = 0; i < mat->r; i++)
            _nmod_poly_vec_pivot_profile(pivind+i, &buf, mat->rows[i], shift, mat->c, orient);
    }
    else
    {
        slong buf;
        // copy column into vec
        nmod_poly_struct * vec = flint_malloc(mat->r * sizeof(nmod_poly_struct));
        for (slong j = 0; j < mat->c; j++)
        {
            for (slong i = 0; i < mat->r; i++)
                vec[i] = mat->rows[i][j];
            _nmod_poly_vec_pivot_profile(pivind+j, &buf, vec, shift, mat->r, orient);
        }
        flint_free(vec);
    }
}

/**********************************************************************
*                       shifted pivot profile                        *
**********************************************************************/

void nmod_poly_mat_pivot_profile(slong * pivind,
                                 slong * pivdeg,
                                 const nmod_poly_mat_t mat,
                                 const slong * shift,
                                 orientation_t orient)
{
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        for (slong i = 0; i < mat->r; i++)
            _nmod_poly_vec_pivot_profile(pivind+i, pivdeg+i, mat->rows[i], shift, mat->c, orient);
    else
    {
        // force access to column as a vec
        nmod_poly_struct * vec = flint_malloc(mat->r * sizeof(nmod_poly_struct));
        for (slong j = 0; j < mat->c; j++)
        {
            for (slong i = 0; i < mat->r; i++)
                vec[i] = mat->rows[i][j];
            _nmod_poly_vec_pivot_profile(pivind+j, pivdeg+j, vec, shift, mat->r, orient);
        }
        flint_free(vec);
    }
}


/**********************************************************************
*                        echelon pivot index                         *
**********************************************************************/

void nmod_poly_mat_echelon_pivot_index(slong * pivind, const nmod_poly_mat_t mat, orientation_t orient)
{
    switch (orient)
    {
        slong piv;
        case ROW_LOWER:
            for (slong i = 0; i < mat->r; i++)
            {
                // rightmost nonzero entry in row i if nonzero, otherwise -1
                piv = mat->c - 1;
                while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
                    piv--;
                pivind[i] = piv;
            }
            break;

        case ROW_UPPER:
            for (slong i = 0; i < mat->r; i++)
            {
                // leftmost nonzero entry in row i if nonzero, otherwise mat->c
                piv = 0;
                while (piv < mat->c && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
                    piv++;
                pivind[i] = piv;
            }
            break;

        case COL_LOWER:
            for (slong j = 0; j < mat->c; j++)
            {
                // topmost nonzero entry in column j if nonzero, otherwise mat->r
                piv = 0;
                while (piv < mat->r && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
                    piv++;
                pivind[j] = piv;
            }
            break;

        case COL_UPPER:
            for (slong j = 0; j < mat->c; j++)
            {
                // bottommost nonzero entry in column j if nonzero, otherwise -1
                piv = mat->r - 1;
                while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
                    piv--;
                pivind[j] = piv;
            }
            break;
    }
}

/**********************************************************************
*                       echelon pivot profile                        *
**********************************************************************/

void nmod_poly_mat_echelon_pivot_profile(slong * pivind, slong * pivdeg, const nmod_poly_mat_t mat, orientation_t orient)
{
    switch (orient)
    {
        slong piv;

    case ROW_LOWER:
        for (slong i = 0; i < mat->r; i++)
        {
            // rightmost nonzero entry in row i if nonzero, otherwise -1
            piv = mat->c - 1;
            while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
                piv--;
            pivind[i] = piv;
            pivdeg[i] = (piv==-1) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, i, piv));
        }
        break;

    case ROW_UPPER:
        for (slong i = 0; i < mat->r; i++)
        {
            // leftmost nonzero entry in row i if nonzero, otherwise mat->c
            piv = 0;
            while (piv < mat->c && nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, piv)))
                piv++;
            pivind[i] = piv;
            pivdeg[i] = (piv==mat->c) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, i, piv));
        }
        break;

    case COL_LOWER:
        for (slong j = 0; j < mat->c; j++)
        {
            // topmost nonzero entry in column j if nonzero, otherwise mat->r
            piv = 0;
            while (piv < mat->r && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
                piv++;
            pivind[j] = piv;
            pivdeg[j] = (piv == mat->r) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, piv, j));
        }
        break;

    case COL_UPPER:
        for (slong j = 0; j < mat->c; j++)
        {
            // bottommost nonzero entry in column j if nonzero, otherwise -1
            piv = mat->r - 1;
            while (piv >= 0 && nmod_poly_is_zero(nmod_poly_mat_entry(mat, piv, j)))
                piv--;
            pivind[j] = piv;
            pivdeg[j] = (piv==-1) ? -1 : nmod_poly_degree(nmod_poly_mat_entry(mat, piv, j));
        }
        break;
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
