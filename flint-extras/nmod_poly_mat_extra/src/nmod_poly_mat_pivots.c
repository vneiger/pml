#include "nmod_poly_mat_forms.h"

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT INDEX/DEGREE                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void leading_positions(slong *res, const nmod_poly_mat_t mat,
                       const slong *shifts, orientation_t row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max;
    slong d;
    slong ind; // TODO initial value
    if (row_wise)
    {
        for (slong i = 0; i < rdim; i++)
        {
            max = 0;
            for (slong j = 0; j < cdim; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
                if (max <= d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
        return;
    }
    {
        for (slong i = 0; i < cdim; i++)
        {
            max = 0;
            for (slong j = 0; j < rdim; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
                if (max <= d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
    }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*    LEADING MATRIX                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void leading_matrix(nmod_mat_t lmat, const nmod_poly_mat_t mat, orientation_t orient)
{
    if (orient == ROW_WISE)
        leading_matrix_rowwise(lmat, mat);
    else if (orient == COLUMN_WISE)
        leading_matrix_columnwise(lmat, mat);
    // TODO add failure when not one of these two?
}

void leading_matrix_rowwise(nmod_mat_t lmat,
                            const nmod_poly_mat_t mat)
{
    // find row degrees
    slong rdeg[mat->r];
    row_degrees(rdeg, mat);

    // deduce leading matrix
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_mat_set_entry(lmat, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                      rdeg[i]));
    // TODO check whether bug when reaching zero row (rdeg[i] should be -1)
}

void leading_matrix_columnwise(nmod_mat_t lmat,
                               const nmod_poly_mat_t mat)
{
    // find column degrees
    slong cdeg[mat->c];
    column_degrees(cdeg, mat);

    // deduce leading matrix
    for (slong j = 0; j < mat->c; j++)
        for (slong i = 0; i < mat->r; i++)
            nmod_mat_set_entry(lmat, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                      cdeg[j]));
    // TODO check whether bug when reaching zero row (rdeg[i] should be -1)
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*   SHIFTED LEADING MATRIX                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void leading_matrix_shifted(nmod_mat_t lmat, const nmod_poly_mat_t mat,
                    const slong *shifts, orientation_t orient)
{
    if (orient == ROW_WISE)
        leading_matrix_shifted_rowwise(lmat, mat, shifts);
    else if (orient == COLUMN_WISE)
        leading_matrix_shifted_columnwise(lmat, mat, shifts);
    // TODO add failure when not one of these two?
}


void leading_matrix_shifted_rowwise(nmod_mat_t lmat,
                                    const nmod_poly_mat_t mat,
                                    const slong *shifts)
{
    // find row degrees
    slong rdeg[mat->r];
    row_degrees_shifted(rdeg, mat, shifts);

    // deduce leading matrix
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_mat_set_entry(lmat, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                      rdeg[i] - shifts[j]));
    // TODO check whether bug when reaching min_shift-1 (rdeg[i] - shifts[j] might be negative?)
}

void leading_matrix_shifted_columnwise(nmod_mat_t lmat,
                                       const nmod_poly_mat_t mat,
                                       const slong *shifts)
{
    // find column degrees
    slong cdeg[mat->c];
    column_degrees_shifted(cdeg, mat, shifts);

    // deduce leading matrix
    for (slong j = 0; j < mat->c; j++)
        for (slong i = 0; i < mat->r; i++)
            nmod_mat_set_entry(lmat, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j),
                                                      cdeg[j] - shifts[i]));
    // TODO check whether bug when reaching min_shift-1 (rdeg[i] - shifts[j] might be negative?)
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
