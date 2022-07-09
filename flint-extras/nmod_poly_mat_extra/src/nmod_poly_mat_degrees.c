#include "nmod_poly_mat_forms.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) ROW/COLUMN DEGREE                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


void row_degrees(slong *rdeg, const nmod_poly_mat_t mat)
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

void column_degrees(slong *cdeg, const nmod_poly_mat_t mat)
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

void row_degrees_shifted(slong *rdeg, const nmod_poly_mat_t mat, const slong *shift)
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

void column_degrees_shifted(slong *cdeg, const nmod_poly_mat_t mat, const slong *shift)
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


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT INDEX                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void pivot_index_rowwise(slong *pivind, const nmod_poly_mat_t mat)
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

void pivot_index_columnwise(slong *pivind, const nmod_poly_mat_t mat)
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

void pivot_index_shifted_rowwise(slong *pivind, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;
    slong min_shift = (mat->c > 0) ? shift[0] : 0;

    // find minimum of shift
    for (slong j = 0; j < mat->c; ++j)
        if (shift[j] < min_shift)
            min_shift = shift[j];

    for (slong i = 0; i < mat->r; i++)
    {
        max = min_shift-1; // zero rows will have this as rdeg
        piv = -1;
        for (slong j = 0; j < mat->c; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            if (0 <= d && max <= d+shift[j])
            {
                max = d;
                piv = j;
            }
        }
        pivind[i] = piv;
    }
}

void pivot_index_shifted_columnwise(slong *pivind, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;
    slong min_shift = (mat->r > 0) ? shift[0] : 0;

    // find minimum of shift
    for (slong i = 0; i < mat->r; ++i)
        if (shift[i] < min_shift)
            min_shift = shift[i];

    for (slong j = 0; j < mat->c; j++)
    {
        max = min_shift-1;
        piv = -1;
        for (slong i = 0; i < mat->r; i++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            if (0 <= d && max <= d+shift[i])
            {
                max = d;
                piv = i;
            }
        }
        pivind[j] = piv;
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT PROFILE                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void pivot_profile_rowwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat)
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

void pivot_profile_columnwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat)
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

void pivot_profile_shifted_rowwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;
    slong min_shift = (mat->c > 0) ? shift[0] : 0;

    // find minimum of shift
    for (slong j = 0; j < mat->c; ++j)
        if (shift[j] < min_shift)
            min_shift = shift[j];

    for (slong i = 0; i < mat->r; i++)
    {
        max = min_shift-1; // zero rows will have this as rdeg
        piv = -1;
        for (slong j = 0; j < mat->c; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            if (0 <= d && max <= d+shift[j])
            {
                max = d;
                piv = j;
            }
        }
        pivdeg[i] = max;
        pivind[i] = piv;
    }
}

void pivot_profile_shifted_columnwise(slong *pivind, slong *pivdeg, const nmod_poly_mat_t mat, const slong *shift)
{
    slong max, piv, d;
    slong min_shift = (mat->r > 0) ? shift[0] : 0;

    // find minimum of shift
    for (slong i = 0; i < mat->r; ++i)
        if (shift[i] < min_shift)
            min_shift = shift[i];

    for (slong j = 0; j < mat->c; j++)
    {
        max = min_shift-1;
        piv = -1;
        for (slong i = 0; i < mat->r; i++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            // if new maximum (or equal maximum) reached at a nonzero entry, update
            if (0 <= d && max <= d+shift[i])
            {
                max = d;
                piv = i;
            }
        }
        pivdeg[j] = max;
        pivind[j] = piv;
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) DEGREE MATRIX                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void degree_matrix(fmpz_mat_t dmat,
                   const nmod_poly_mat_t mat)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
}

void degree_matrix_row_shifted(fmpz_mat_t dmat,
                               const nmod_poly_mat_t mat,
                               const slong * shift)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shift[j];
}

void degree_matrix_column_shifted(fmpz_mat_t dmat,
                                  const nmod_poly_mat_t mat,
                                  const slong * shift)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shift[i];
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) LEADING MATRIX                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

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
