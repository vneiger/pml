#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"

void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, slong degree)
{
    slong rdim = mat->r, cdim = mat->c;
    for(slong i = 0; i < rdim; i++)
        for(slong j = 0; j < cdim; j++)
            nmod_mat_set_entry(res, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j), degree));
}

void nmod_poly_mat_shift(nmod_poly_mat_t res, slong k)
{
    nmod_poly_struct *P;
    if (k > 0)
    {
        for (slong i = 0; i < res->r; i++)
            for (slong j = 0; j < res->c; j++)
            {
                P = nmod_poly_mat_entry(res, i, j);
                if (!nmod_poly_is_zero(P))
                    nmod_poly_shift_left(P, P, k);
            }
        return;
    }


    if (k < 0)
    {
        for (slong i = 0; i < res->r; i++)
            for (slong j = 0; j < res->c; j++)
            {
                P = nmod_poly_mat_entry(res, i, j);
                if (!nmod_poly_is_zero(P))
                    nmod_poly_shift_right(P, P, k);
            }
        return;
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
