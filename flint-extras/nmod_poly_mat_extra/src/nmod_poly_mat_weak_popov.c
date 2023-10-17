#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // TODO remove, for debugging

#define MAT(i,j) (mat->rows[i] + j)
#define TSF(i,j) (tsf->rows[i] + j)
#define OTHER(i,j) (other->rows[i] + j)

/**********************************************************************
*                          HELPER FUNCTIONS                          *
**********************************************************************/

// Context: row-wise
// two "shifted weak Popov pivots" have been found in the same column j, i.e.
// j is the shift-pivot index of two distinct rows pi1 and pi2 (thus pi1 != pi2).
// Without loss of generality we assume deg [pi1,j] >= deg [pi2,j].
// -> use pi2,j to kill the leading term of pi1,j, and apply the
// corresponding transformation to the whole row pi1
// --> this is a "simple transformation of the first kind" (terminology from Mulders & Storjohann 2003)
//         mat[pi1,:] = mat[pi1,:] + cst * x**exp * mat[pi2,:]
// where cst = - leading_coeff(mat[pi1,j]) / leading_coeff(mat[pi2,j])
//   and exp = deg(mat[pi1,j]) - deg(mat[pi2,j])
void _atomic_solve_pivot_collision_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                           slong pi1, slong pi2, slong j)
{
    // compute exp = deg(mat[pi1,j]) - deg(mat[pi2,j])
    const slong exp = MAT(pi1, j)->length - MAT(pi2, j)->length;
    // compute cst = - leading_coeff(mat[pi1,j]) / leading_coeff(mat[pi2,j])
    mp_limb_t cst = n_invmod(MAT(pi2, j)->coeffs[MAT(pi2, j)->length -1], mat->modulus);
    cst = nmod_mul(MAT(pi1, j)->coeffs[MAT(pi1, j)->length -1], cst, MAT(pi1, j)->mod);
    cst = nmod_neg(cst, MAT(pi1, j)->mod);
    // update other entries in row pi1
    for (slong jj = 0; jj < mat->c; jj++)
    {
        const slong len = MAT(pi2, jj)->length + exp;
        nmod_poly_fit_length(MAT(pi1, jj), len);
        if (len > MAT(pi1, jj)->length)
        {
            _nmod_vec_zero(MAT(pi1, jj)->coeffs + MAT(pi1, jj)->length, len - MAT(pi1, jj)->length);
            _nmod_poly_set_length(MAT(pi1, jj), len);
        }
        _nmod_vec_scalar_addmul_nmod(MAT(pi1, jj)->coeffs+exp, MAT(pi2, jj)->coeffs, MAT(pi2, jj)->length, cst, MAT(pi1, jj)->mod);
        _nmod_poly_normalise(MAT(pi1, jj));
    }
    // update the row pi1 of other
    if (other) for (slong jj = 0; jj < other->c; jj++)
    {
        const slong len = OTHER(pi2, jj)->length + exp;
        nmod_poly_fit_length(OTHER(pi1, jj), len);
        if (len > OTHER(pi1, jj)->length)
        {
            _nmod_vec_zero(OTHER(pi1, jj)->coeffs + OTHER(pi1, jj)->length, len - OTHER(pi1, jj)->length);
            _nmod_poly_set_length(OTHER(pi1, jj), len);
        }
        _nmod_vec_scalar_addmul_nmod(OTHER(pi1, jj)->coeffs+exp, OTHER(pi2, jj)->coeffs, OTHER(pi2, jj)->length, cst, OTHER(pi1, jj)->mod);
        _nmod_poly_normalise(OTHER(pi1, jj));
    }
}

// Context: row-wise
// reduce entry mat[ii,j] against pivot entry which is at mat[i,j],
// applying the corresponding operation to the whole row ii of mat and the
// corresponding row ii of other
// u, v are used as temporaries and must be already initialized
// other == NULL or other is another matrix with the same modulus (it must have a row ii)
void _reduce_against_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
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


// Context: row-wise
// normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of other
// (i,j) is position of the pivot entry
// other == NULL or other is another matrix with the same modulus
void _normalize_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j)
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

/**********************************************************************
*                     weak Popov form algorithms                     *
**********************************************************************/

// Orientation: lower, row-wise
// iterative weak Popov form algorithm by Mulders and Storjohann, 2003 there
// presented for the uniform shift, here straightforwardly adapted to the
// shifted case.
//
// Introduce rows one after another, only adding next one when current matrix
// is already in shift-ordered weak Popov form. When discovering a zero row, a
// row rotation is used to push it back as the last row.
//
// Strategy for pivot collision selection: target the collision involving the
// smallest possible degree (this means fewer field operations to do for
// transforming `wpf`, but also means greater degrees in the unimodular
// transformation or in `tsf`, which may impact performance if `tsf` is not
// NULL)
//
// pivind must be allocated with at least mat->r entries; it will be populated
// with the shifted pivot index of the output weak Popov form (undefined
// behaviour for entries beyond mat->r)
slong nmod_poly_mat_weak_popov_mulders_storjohann_lower_rowwise(nmod_poly_mat_t mat,
                                                                const slong * shift,
                                                                slong * pivind,
                                                                nmod_poly_mat_t tsf)
{
    const slong m = mat->r;
    const slong n = mat->c;

    if (m == 0 || n == 0)
        return 0;

    // number of found pivots (i.e., rank), and number of found zero rows
    // all along, mat->r == rk + zr + number of rows that remained to be processed
    slong rk = 0;
    slong zr = 0;

    //flint_mpn_store(pivind, mat->r, -1); // fill with -1, presumably all rows are zero

    while (rk + zr < mat->r)
    {
        // consider next row, compute its pivot index
        slong pivdeg;
        _nmod_poly_vec_pivot_profile(pivind+rk, &pivdeg, mat->rows[rk], shift, mat->c);
        if (pivind[rk] == -1)
        {
            // row is zero, rotate to put it last; increment zr
            _nmod_poly_mat_rotate_rows_upward(mat, pivind, rk, mat->r-1);
            zr++;
        }
        else
        {
            // row is nonzero, let's actually do some work
        }
    }
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
