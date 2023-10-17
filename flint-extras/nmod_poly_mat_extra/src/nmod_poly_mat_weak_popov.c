#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // TODO remove, for debugging

#define MAT(i,j) (mat->rows[i] + j)
#define TSF(i,j) (tsf->rows[i] + j)
#define OTHER(i,j) (other->rows[i] + j)

// normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of other
// (i,j) is position of the pivot entry
// other == NULL or other is another matrix with the same modulus
void _normalize_pivot_general(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j)
{
    if (! nmod_poly_is_monic(MAT(i, j)))
    {
        mp_limb_t inv = n_invmod(MAT(i, j)->coeffs[MAT(i, j)->length - 1], MAT(i, j)->mod.n);
        for (slong jj = 0; jj < mat->c; jj++)
            _nmod_vec_scalar_mul_nmod(MAT(i, jj)->coeffs, MAT(i, jj)->coeffs, MAT(i, jj)->length, inv, MAT(i, jj)->mod);
        if (other)
            for (slong jj = 0; jj < other->c; jj++)
                _nmod_vec_scalar_mul_nmod(OTHER(i, jj)->coeffs, OTHER(i, jj)->coeffs, OTHER(i, jj)->length, inv, OTHER(i, jj)->mod);
    }
}

// reduce entry mat[ii,j] against pivot entry which is at mat[i,j],
// applying the corresponding operation to the whole row ii of mat and the
// corresponding row ii of other
// u, v are used as temporaries and must be already initialized
// other == NULL or other is another matrix with the same modulus
void _reduce_against_pivot_general(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                   slong i, slong j, slong ii,
                                   nmod_poly_t u, nmod_poly_t v)
{
    if (MAT(ii, j)->length >= MAT(i, j)->length)
    {
        // division with remainder: [ii,j] = pivot * u + v
        nmod_poly_divrem(u, v, MAT(ii, j), MAT(i, j));
        // replace [ii,j] by remainder
        nmod_poly_swap(MAT(ii, j), v);
        // apply same transformation on rest of the row: row_ii <- row_ii - u * row_i
        nmod_poly_neg(u, u);
        for (slong jj = 0; jj < j; jj++)
        {
            nmod_poly_mul(v, u, MAT(i, jj));
            nmod_poly_add(MAT(ii, jj), MAT(ii, jj), v);
        }
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            nmod_poly_mul(v, u, MAT(i, jj));
            nmod_poly_add(MAT(ii, jj), MAT(ii, jj), v);
        }
        // apply same transformation on row ii of other
        if (other)
        {
            for (slong jj = 0; jj < other->c; jj++)
            {
                nmod_poly_mul(v, u, OTHER(i, jj));
                nmod_poly_add(OTHER(ii, jj), OTHER(ii, jj), v);
            }
        }
    }
}


slong nmod_poly_mat_weak_popov_mulders_storjohann_upper_rowwise(nmod_poly_mat_t wpf,
                                                                const slong * shift,
                                                                nmod_poly_mat_t tsf);

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
