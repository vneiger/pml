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
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_forms.h" // for nmod_poly_mat_degree_matrix, orientation

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PRETTY PRINTING THE MATRIX                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_print_pretty(const nmod_poly_mat_t mat, const char * var)
{
    slong rdim = mat->r, cdim = mat->c;

    flint_printf("<%wd x %wd matrix over Z/nZ[%s]>\n", mat->r, mat->c, var);
    flint_printf("[");
    for (slong i = 0; i < rdim; i++)
    {
        flint_printf("[");
        for (slong j = 0; j < cdim; j++)
        {
            nmod_poly_print_pretty(nmod_poly_mat_entry(mat, i, j), var);
            if (j+1 < cdim)
                flint_printf(", ");
        }
        if (i != rdim -1)
            flint_printf("],\n");
        else
            flint_printf("]");
    }
    flint_printf("]\n");
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PRETTY PRINTING DEGREE MATRIX                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_degree_matrix_print_pretty(const nmod_poly_mat_t mat)
{
    fmpz_mat_t dmat;
    fmpz_mat_init(dmat, mat->r, mat->c);
    nmod_poly_mat_degree_matrix(dmat, mat); 
    fmpz_mat_print_pretty(dmat);
    printf("\n");
    fmpz_mat_clear(dmat);
}

void nmod_poly_mat_degree_matrix_shifted_print_pretty(const nmod_poly_mat_t mat,
                                                      const slong *shift,
                                                      orientation_t orient)
{
    fmpz_mat_t dmat;
    fmpz_mat_init(dmat, mat->r, mat->c);
    nmod_poly_mat_degree_matrix_shifted(dmat, mat, shift, orient); 
    fmpz_mat_print_pretty(dmat);
    printf("\n");
    fmpz_mat_clear(dmat);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PRETTY PRINTING LEADING MATRIX                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_leading_matrix_print_pretty(const nmod_poly_mat_t mat,
                                               const slong * shift,
                                               orientation_t orient)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix(lmat, mat, shift, orient);
    nmod_mat_print_pretty(lmat);
    nmod_mat_clear(lmat);
}
