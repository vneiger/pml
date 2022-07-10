#include "nmod_poly_mat_forms.h"

void nmod_poly_mat_leading_matrix_rowwise(nmod_mat_t lmat,
                            const nmod_poly_mat_t mat)
{
    // find row degrees
    slong rdeg[mat->r];
    nmod_poly_mat_row_degree(rdeg, mat);

    // deduce leading matrix
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

void nmod_poly_mat_leading_matrix_columnwise(nmod_mat_t lmat,
                               const nmod_poly_mat_t mat)
{
    // find column degrees
    slong cdeg[mat->c];
    nmod_poly_mat_column_degree(cdeg, mat);

    // deduce leading matrix
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

void nmod_poly_mat_leading_matrix_shifted_rowwise(nmod_mat_t lmat,
                                    const nmod_poly_mat_t mat,
                                    const slong *shifts)
{
    // find row degrees
    slong rdeg[mat->r];
    nmod_poly_mat_row_degree_shifted(rdeg, mat, shifts);

    // deduce leading matrix
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            if (rdeg[i] >= shifts[j])
                nmod_mat_set_entry(lmat, i, j,
                                   nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                          rdeg[i] - shifts[j]));
            else
                nmod_mat_set_entry(lmat, i, j, 0);
}

void nmod_poly_mat_leading_matrix_shifted_columnwise(nmod_mat_t lmat,
                                       const nmod_poly_mat_t mat,
                                       const slong *shifts)
{
    // find column degrees
    slong cdeg[mat->c];
    nmod_poly_mat_column_degree_shifted(cdeg, mat, shifts);

    // deduce leading matrix
    for (slong j = 0; j < mat->c; j++)
        for (slong i = 0; i < mat->r; i++)
            if (cdeg[j] >= shifts[i])
                nmod_mat_set_entry(lmat, i, j,
                                   nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                          cdeg[j] - shifts[i]));
            else
                nmod_mat_set_entry(lmat, i, j, 0);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
