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
#include <flint/nmod_poly.h>
#include "nmod_poly_mat_forms.h"

void nmod_poly_mat_degree_matrix(fmpz_mat_t dmat,
                                 const nmod_poly_mat_t mat)
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
