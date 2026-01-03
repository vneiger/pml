/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include "nmod_poly_mat_forms.h"

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
