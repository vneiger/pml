/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

#include "nmod_poly_mat_forms.h"

/*------------------------------------------------------------*/
/*        row degree - column degree - degree matrix          */
/*------------------------------------------------------------*/

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

void nmod_poly_mat_degree_matrix(fmpz_mat_t dmat, const nmod_poly_mat_t mat)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
}

void nmod_poly_mat_degree_matrix_shifted(fmpz_mat_t dmat,
                                         const nmod_poly_mat_t mat,
                                         const slong * shift,
                                         orientation_t orient)
{
    if (shift)
    {
        if (orient == ROW_LOWER || orient == ROW_UPPER)
            for(slong i = 0; i < mat->r; i++)
                for(slong j = 0; j < mat->c; j++)
                    *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shift[j];
        else // orient == COL_*
            for(slong i = 0; i < mat->r; i++)
                for(slong j = 0; j < mat->c; j++)
                    *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shift[i];
    }
    else
        for(slong i = 0; i < mat->r; i++)
            for(slong j = 0; j < mat->c; j++)
                *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
}


/*------------------------------------------------------------*/
/*          shifted pivot profile (index + degree)            */
/*------------------------------------------------------------*/

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

void nmod_poly_mat_pivot_index(slong *pivind,
                               const nmod_poly_mat_t mat,
                               const slong * shift,
                               orientation_t orient)
{
    if (orient == ROW_LOWER || orient == ROW_UPPER)
    {
        slong buf;
        for (slong i = 0; i < mat->r; i++)
            _nmod_poly_vec_pivot_profile(pivind+i, &buf, nmod_poly_mat_entry(mat, i, 0), shift, mat->c, orient);
    }
    else
    {
        slong buf;
        // copy column into vec
        nmod_poly_struct * vec = flint_malloc(mat->r * sizeof(nmod_poly_struct));
        for (slong j = 0; j < mat->c; j++)
        {
            for (slong i = 0; i < mat->r; i++)
                vec[i] = *nmod_poly_mat_entry(mat, i, j);
            _nmod_poly_vec_pivot_profile(pivind+j, &buf, vec, shift, mat->r, orient);
        }
        flint_free(vec);
    }
}

void nmod_poly_mat_pivot_profile(slong * pivind,
                                 slong * pivdeg,
                                 const nmod_poly_mat_t mat,
                                 const slong * shift,
                                 orientation_t orient)
{
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        for (slong i = 0; i < mat->r; i++)
            _nmod_poly_vec_pivot_profile(pivind+i, pivdeg+i, nmod_poly_mat_entry(mat, i, 0), shift, mat->c, orient);
    else
    {
        // force access to column as a vec
        nmod_poly_struct * vec = flint_malloc(mat->r * sizeof(nmod_poly_struct));
        for (slong j = 0; j < mat->c; j++)
        {
            for (slong i = 0; i < mat->r; i++)
                vec[i] = *nmod_poly_mat_entry(mat, i, j);
            _nmod_poly_vec_pivot_profile(pivind+j, pivdeg+j, vec, shift, mat->r, orient);
        }
        flint_free(vec);
    }
}

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


/*------------------------------------------------------------*/
/*                      leading matrix                        */
/*------------------------------------------------------------*/

void nmod_poly_mat_leading_matrix(nmod_mat_t lmat,
                                  const nmod_poly_mat_t mat,
                                  const slong * shift,
                                  orientation_t orient)
{
    if (orient == ROW_LOWER || orient == ROW_UPPER)
    {
        // find row degrees
        slong rdeg[mat->r];
        nmod_poly_mat_row_degree(rdeg, mat, shift);

        // deduce leading matrix
        if (shift == NULL)
        {
            for (slong i = 0; i < mat->r; i++)
                if (rdeg[i] >= 0)
                    for (slong j = 0; j < mat->c; j++)
                        nmod_mat_set_entry(lmat, i, j,
                                           nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                                  rdeg[i]));
                else
                    for (slong j = 0; j < mat->c; j++)
                        nmod_mat_set_entry(lmat, i, j, 0);
        }
        else
        {
            for (slong i = 0; i < mat->r; i++)
                for (slong j = 0; j < mat->c; j++)
                    if (rdeg[i] >= shift[j])
                        nmod_mat_set_entry(lmat, i, j,
                                           nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                                  rdeg[i] - shift[j]));
                    else
                        nmod_mat_set_entry(lmat, i, j, 0);
        }
    }
    else // orient == COL_*
    {
        // find column degrees
        slong cdeg[mat->c];
        nmod_poly_mat_column_degree(cdeg, mat, shift);

        // deduce leading matrix
        if (shift == NULL)
        {
            for (slong j = 0; j < mat->c; j++)
                if (cdeg[j] >= 0)
                    for (slong i = 0; i < mat->r; i++)
                        nmod_mat_set_entry(lmat, i, j,
                                           nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                                  cdeg[j]));
                else
                    for (slong i = 0; i < mat->r; i++)
                        nmod_mat_set_entry(lmat, i, j, 0);
        }
        else
        {
            for (slong j = 0; j < mat->c; j++)
                for (slong i = 0; i < mat->r; i++)
                    if (cdeg[j] >= shift[i])
                        nmod_mat_set_entry(lmat, i, j,
                                           nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                                  cdeg[j] - shift[i]));
                    else
                        nmod_mat_set_entry(lmat, i, j, 0);
        }
    }
}
