#include "nmod_poly_mat_utils.h"

void nmod_mat_permute_rows(nmod_mat_t mat, const slong * perm, slong rdim)
{
    slong *perm_copy = _perm_init(rdim);
    _perm_set(perm_copy, perm, rdim);
    for (slong i = 0; i < rdim; i++)
        for (slong j = 0; j < rdim; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
        }

    _perm_clear(perm_copy);
}

void nmod_poly_mat_permute_columns(nmod_poly_mat_t mat, const slong * perm, slong cdim)
{
    slong *perm_copy = _perm_init(cdim);
    _perm_inv(perm_copy, perm, cdim);
    for (slong i = 0; i < cdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_poly_mat_swap_cols(mat, perm_copy, i, perm_copy[i]);
        }
    }
    _perm_clear(perm_copy);
}

void nmod_poly_mat_permute_rows(nmod_poly_mat_t mat, const slong *perm, slong rdim)
{
    slong *perm_copy = _perm_init(rdim);
    _perm_set(perm_copy, perm, rdim);

    for (slong i = 0; i < rdim; i++)
        for (slong j = 0; j < rdim; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_poly_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
        }
    _perm_clear(perm_copy);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
