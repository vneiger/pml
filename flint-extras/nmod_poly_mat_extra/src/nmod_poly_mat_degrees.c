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

void row_degrees_shifted(slong *rdeg, const nmod_poly_mat_t mat, const slong *shifts)
{
    slong max, d;
    slong min_shift = (mat->c > 0) ? shifts[0] : 0;

    // find minimum of shift
    for (slong j = 0; j < mat->c; ++j)
        if (shifts[j] < min_shift)
            min_shift = shifts[j];

    for (slong i = 0; i < mat->r; i++)
    {
        max = min_shift-1; // zero rows will have this as rdeg
        for (slong j = 0; j < mat->c; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
            // if new maximum at a nonzero entry, update
            if (shifts[j] <= d && max < d)
                max = d;
        }
        rdeg[i] = max;
    }
}

void column_degrees_shifted(slong *cdeg, const nmod_poly_mat_t mat, const slong *shifts)
{
    slong max, d;
    slong min_shift = (mat->r > 0) ? shifts[0] : 0;

    // find minimum of shift
    for (slong i = 0; i < mat->r; ++i)
        if (shifts[i] < min_shift)
            min_shift = shifts[i];

    for (slong j = 0; j < mat->c; j++)
    {
        max = min_shift-1;
        for (slong i = 0; i < mat->r; i++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[i];
            // if new maximum at a nonzero entry, update
            if (shifts[i] <= d && max < d)
                max = d;
        }
        cdeg[j] = max;
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT INDEX/DEGREE                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) DEGREE MATRIX                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void degree_matrix(slong *res,
                   const nmod_poly_mat_t mat)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *(res + (j * mat->c) + i) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
    // TODO use access provided by fmpz matrix
}

void degree_matrix_row_shifted(int64_t *res,
                               const nmod_poly_mat_t mat,
                               const int64_t *shifts)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *(res + (i * mat->r) + j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
    // TODO use access provided by fmpz matrix
    // TODO manage zero entry correctly
}

void degree_matrix_column_shifted(int64_t *res,
                                  const nmod_poly_mat_t mat,
                                  const int64_t *shifts)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *(res + (j * mat->c) + i) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[i];
    // TODO use access provided by fmpz matrix
    // TODO manage zero entry correctly
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
