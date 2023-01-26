#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* RANDOM - UNIFORM, DENSE                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_rand(nmod_poly_mat_t mat,
                        flint_rand_t state,
                        slong len)
{
    for (slong i = 0; i < mat->r * mat->c; i++)
        nmod_poly_rand(mat->entries + i, state, len);
}

void nmod_poly_mat_rand_row_degree(nmod_poly_mat_t mat,
                                   flint_rand_t state,
                                   const slong * rdeg)
{
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_rand(mat->rows[i] + j, state, rdeg[i]+1);
}

void nmod_poly_mat_rand_column_degree(nmod_poly_mat_t mat,
                                      flint_rand_t state,
                                      const slong * cdeg)
{
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_rand(mat->rows[i] + j, state, cdeg[j]+1);
}

void nmod_poly_mat_rand_degree_matrix(nmod_poly_mat_t mat,
                                      flint_rand_t state,
                                      const fmpz_mat_t dmat)
{
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_rand(mat->rows[i] + j, state, 1 + *fmpz_mat_entry(dmat, i, j));
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* RANDOM - POPOV                                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_rand_popov_rowwise(nmod_poly_mat_t mat,
                                      flint_rand_t state,
                                      const slong * pivind,
                                      const slong * pivdeg,
                                      const slong * shift)
{
    // will store, for each row, the corresponding degree constraints
    slong * dvec = (slong *) flint_malloc(mat->c * sizeof(slong));

    // process each row
    for (slong i = 0; i < mat->r; i++)
    {
        // row-wise constraints for row i
        for (slong j = 0; j < pivind[i]; j++)
            dvec[j] = pivdeg[i]+shift[i]-shift[j];
        for (slong j = pivind[i]+1; j < mat->c; j++)
            dvec[j] = pivdeg[i]+shift[i]-shift[j]-1;
        // column-wise constraints for row i
        for (slong ii = 0; ii < mat->r; ii++)
            if (dvec[pivind[ii]] >= pivdeg[ii])
                dvec[pivind[ii]] = pivdeg[ii]-1;
        // random filling before pivot
        for (slong j = 0; j < pivind[i]; j++)
            nmod_poly_rand(mat->rows[i] + j, state, 1 + dvec[j]);
        // add monic of pivot degree at pivot
        nmod_poly_rand_monic(mat->rows[i] + pivind[i], state, 1+pivdeg[i]);
        // random filling after pivot
        for (slong j = pivind[i]+1; j < mat->c; j++)
            nmod_poly_rand(mat->rows[i] + j, state, 1 + dvec[j]);
    }

    free(dvec);
}

void nmod_poly_mat_rand_popov_columnwise(nmod_poly_mat_t mat,
                                         flint_rand_t state,
                                         const slong * pivind,
                                         const slong * pivdeg,
                                         const slong * shift)
{
    // will store, for each column, the corresponding degree constraints
    slong * dvec = (slong *) flint_malloc(mat->r * sizeof(slong));

    // process each column
    for (slong j = 0; j < mat->c; j++)
    {
        // column-wise constraints for column j
        for (slong i = 0; i < pivind[j]; i++)
            dvec[i] = pivdeg[j]+shift[j]-shift[i];
        for (slong i = pivind[j]+1; i < mat->r; i++)
            dvec[i] = pivdeg[j]+shift[j]-shift[i]-1;
        // row-wise constraints for column j
        for (slong jj = 0; jj < mat->c; jj++)
            if (dvec[pivind[jj]] >= pivdeg[jj])
                dvec[pivind[jj]] = pivdeg[jj]-1;
        // random filling before pivot
        for (slong i = 0; i < pivind[j]; i++)
            nmod_poly_rand(mat->rows[i] + j, state, 1 + dvec[i]);
        // add monic of pivot degree at pivot
        nmod_poly_rand_monic(mat->rows[pivind[j]] + j, state, 1+pivdeg[j]);
        // random filling after pivot
        for (slong i = pivind[j]+1; i < mat->r; i++)
            nmod_poly_rand(mat->rows[i] + j, state, 1 + dvec[i]);
    }

    free(dvec);
}

void nmod_poly_mat_rand_popov(nmod_poly_mat_t mat,
                         flint_rand_t state,
                         const slong * pivind,
                         const slong * pivdeg,
                         const slong * shift,
                         orientation_t row_wise)
{
    const slong * sshift;
    if (shift) // used provided shift
        sshift = shift;
    else // uniform shift
        sshift = flint_calloc(row_wise ? mat->c : mat->r, sizeof(slong));

    // Note: if pivind == NULL, `mat` must be square, this is not checked
    slong * iota;
    if (pivind == NULL)
    {
        iota = flint_malloc((row_wise ? mat->r : mat->c) * sizeof(slong));
        for (slong i = 0; i < (row_wise ? mat->r : mat->c); i++)
            iota[i] = i;
    }
    const slong * ppivind = (pivind) ? pivind : iota;

    if (row_wise)
        nmod_poly_mat_rand_popov_rowwise(mat, state, ppivind, pivdeg, sshift);
    else
        nmod_poly_mat_rand_popov_columnwise(mat, state, ppivind, pivdeg, sshift);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
