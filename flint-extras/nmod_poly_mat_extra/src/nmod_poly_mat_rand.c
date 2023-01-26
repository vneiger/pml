#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* RANDOM - UNIFORM DEGREES                                   */
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

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
