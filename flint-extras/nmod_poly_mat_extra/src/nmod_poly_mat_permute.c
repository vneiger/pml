#include "nmod_poly_mat_utils.h"

//void nmod_mat_permute_rows(nmod_mat_t mat, const slong * perm, slong rdim)
//{
//    slong *perm_copy = _perm_init(rdim);
//    _perm_set(perm_copy, perm, rdim);
//    for (slong i = 0; i < rdim; i++)
//        for (slong j = 0; j < rdim; j++)
//        {
//            if (i == perm_copy[i])
//                break;
//            nmod_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
//        }
//
//    _perm_clear(perm_copy);
//}

void nmod_mat_permute_rows(nmod_mat_t mat, slong * perm_store, const slong * perm_act)
{
    // perm_store[i] <- perm_store[perm_act[i]] 
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    // rows[i] <- rows[perm_act[i]] 
    mp_limb_t ** mat_tmp = flint_malloc(mat->r * sizeof(mp_limb_t *));
    for (slong i = 0; i < mat->r; ++i)
        mat_tmp[perm_act[i]] = mat->rows[i];
    for (slong i = 0; i < mat->r; ++i)
        mat->rows[i] = mat_tmp[i];

    flint_free(mat_tmp);
}

//        mp_ptr * Atmp;
//        slong * APtmp;
//        slong i;
//
//        Atmp = flint_malloc(sizeof(mp_ptr) * n);
//        APtmp = flint_malloc(sizeof(slong) * n);
//
//        for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i] + offset];
//        for (i = 0; i < n; i++) A->rows[i + offset] = Atmp[i];
//
//        for (i = 0; i < n; i++) APtmp[i] = AP[P[i] + offset];
//        for (i = 0; i < n; i++) AP[i + offset] = APtmp[i];
//
//        flint_free(Atmp);
//        flint_free(APtmp);



void nmod_poly_mat_permute_columns(nmod_poly_mat_t mat, const slong * perm)
{
    slong *perm_copy = _perm_init(mat->c);
    _perm_inv(perm_copy, perm, mat->c);
    for (slong i = 0; i < mat->c; i++)
    {
        for (slong j = 0; j < mat->c; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_poly_mat_swap_cols(mat, perm_copy, i, perm_copy[i]);
        }
    }
    _perm_clear(perm_copy);
}

void nmod_poly_mat_permute_rows(nmod_poly_mat_t mat, const slong *perm)
{
    slong *perm_copy = _perm_init(mat->r);
    _perm_set(perm_copy, perm, mat->r);

    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->r; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_poly_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
        }
    _perm_clear(perm_copy);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
