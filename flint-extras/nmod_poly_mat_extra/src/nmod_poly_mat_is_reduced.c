#include "nmod_poly_mat_forms.h"

int nmod_poly_mat_is_row_reduced(const nmod_poly_mat_t mat)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_rowwise(lmat, mat);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->r == rank_lead;
}

int nmod_poly_mat_is_row_reduced_shifted(const nmod_poly_mat_t mat, const slong *shift)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_shifted_rowwise(lmat, mat, shift);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->r == rank_lead;
}

int nmod_poly_mat_is_column_reduced(const nmod_poly_mat_t mat)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_columnwise(lmat, mat);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->c == rank_lead;
}

int nmod_poly_mat_is_column_reduced_shifted(const nmod_poly_mat_t mat, const slong *shift)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_shifted_columnwise(lmat, mat, shift);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->c == rank_lead;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
