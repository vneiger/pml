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


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
