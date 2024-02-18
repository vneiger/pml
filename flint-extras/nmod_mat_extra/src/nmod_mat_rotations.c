#include "nmod_mat_extra.h"


/**********************************************************************
*                    ROW ROTATION DOWNWARD/UPWARD                    *
**********************************************************************/

void _nmod_mat_rotate_rows_downward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

        mp_limb_t * tmp_mat = mat->rows[j];
        for (slong ii = j; ii > i; ii--)
            mat->rows[ii] = mat->rows[ii-1];
        mat->rows[i] = tmp_mat;
    }
}

void _nmod_mat_rotate_rows_upward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[i];
            for (slong ii = i; ii < j; ii++)
                vec[ii] = vec[ii+1];
            vec[j] = tmp_vec;
        }

        mp_limb_t * tmp_mat = mat->rows[i];
        for (slong ii = i; ii < j; ii++)
            mat->rows[ii] = mat->rows[ii+1];
        mat->rows[j] = tmp_mat;
    }
}


/**********************************************************************
*                COLUMN ROTATION RIGHTWARD/LEFTWARD                  *
**********************************************************************/

void _nmod_mat_rotate_columns_rightward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

        for (slong ii = 0; ii < mat->r; ii++)
        {
            mp_limb_t tmp_mat = nmod_mat_entry(mat, ii, j);
            for (slong jj = j; jj > i; jj--)
                nmod_mat_entry(mat, ii, jj) = nmod_mat_entry(mat, ii, jj-1);
            nmod_mat_entry(mat, ii, i) = tmp_mat;
        }
    }
}

void _nmod_mat_rotate_columns_leftward(nmod_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[i];
            for (slong ii = i; ii < j; ii++)
                vec[ii] = vec[ii+1];
            vec[j] = tmp_vec;
        }

        for (slong ii = 0; ii < mat->r; ii++)
        {
            mp_limb_t tmp_mat = nmod_mat_entry(mat, ii, i);
            for (slong jj = i; jj < j; jj++)
                nmod_mat_entry(mat, ii, jj) = nmod_mat_entry(mat, ii, jj+1);
            nmod_mat_entry(mat, ii, i) = tmp_mat;
        }
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
