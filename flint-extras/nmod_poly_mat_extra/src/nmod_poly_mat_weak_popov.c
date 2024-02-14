#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"

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
// Complexity (field operations):
//  - for mat: sum_{jj = 0 ... mat->c-1} (deg(mat[pi2,jj]) + 1)
//             <=  mat->c * (rdeg(mat[pi2,:])+1)
//  - for other: sum_{jj = 0 ... other->c-1} (deg(other[pi2,jj]) + 1)
//             <=  other->c * (rdeg(other[pi2,:])+1)
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
// Complexity: 0 if deg(mat[ii,j]) < deg(mat[i,j]), otherwise:
//    . divrem: M(deg(mat[i,j]))
//    . for mat: sum_{jj = 0 .. j-1 and j+1 ... mat->c-1} M(deg(mat[ii,j]) - deg(mat[i,j]) + deg(mat[i,jj]))
//        <=     (mat->c - 1) * M(deg(mat[ii,j]) - deg(mat[i,j]) + rdeg(mat[i,:]))
//    . for other: sum_{jj = 0 ... other->c-1} M(deg(mat[ii,j]) - deg(mat[i,j]) + deg(other[i,jj]))
//        <=     other->c * M(deg(mat[ii,j]) - deg(mat[i,j]) + rdeg(other[i,:]))
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
// . normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of other
// . return the constant used for normalization
// (i,j) is position of the pivot entry
// other == NULL or other is another matrix with the same modulus
// Complexity: if already monic, 0, otherwise:
//    . for mat: sum_{jj = 0 ... mat->c-1} deg(mat[i,jj])
//          <= mat->c * rdeg(mat[i,:])
//    . for other: sum_{jj = 0 ... other->c-1} deg(other[i,jj])
//          <= other->c * rdeg(other[i,:])
mp_limb_t _normalize_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j)
{
    if (nmod_poly_is_monic(MAT(i, j)))
        return 1UL;
    else
    {
        mp_limb_t inv = n_invmod(MAT(i, j)->coeffs[MAT(i, j)->length - 1], MAT(i, j)->mod.n);
        for (slong jj = 0; jj < mat->c; jj++)
            _nmod_vec_scalar_mul_nmod(MAT(i, jj)->coeffs, MAT(i, jj)->coeffs, MAT(i, jj)->length, inv, MAT(i, jj)->mod);
        if (other)
            for (slong jj = 0; jj < other->c; jj++)
                _nmod_vec_scalar_mul_nmod(OTHER(i, jj)->coeffs, OTHER(i, jj)->coeffs, OTHER(i, jj)->length, inv, OTHER(i, jj)->mod);
        return inv;
    }
}


/**********************************************************************
*                     weak Popov form algorithms                     *
**********************************************************************/

// Introduce rows one after another, only adding next one when current matrix
// is already in shift-ordered weak Popov form. Track the number rk of added
// rows, and the number zr of encountered zero rows. Consequence: at any
// moment, we have encountered rk+zr rows, the zr zero ones have been put at
// the bottom of the matrix, and the rk ones are at the top, with the rk-1
// first ones being full rank and in weak Popov form, and the rk-th one being
// the one we are currently working on.
//
// Therefore in the currently considered submatrix there can be at most one
// pivot collision, and up to swapping the two involved rows (among which is
// rk), we ensure that the atomic collision transformation is modifying row rk.
// This transformation will either keep a collision at the same column, or move
// a pivot resulting to another single collision, or move a pivot but not
// resulting in a collision. After sufficiently many steps, we must end up in
// the last situation (see Mulders and Storjohann for proofs and bounds), which
// splits in two cases: either the rk-th row is zero (we put it at the bottom
// and increment zr) or the rk-th row provides a new pivot (we leave it here
// and increment rk).
//
// In this row-by-row framework, the strategy for pivot collision selection is
// imposed, there is actually no choice.
//
// This applies the algorithm to the submatrix mat[rstart:rstart+rdim,cstart:cstart+cdim],
// called submat below. The left unimodular transformations are applied to the whole
// of mat[rstart:rstart+rdim,:], i.e. not restricting the columns to those of submat.
// tsf must be NULL or have at least rstart+rdim rows, the transformation will
// be applied to its rows rstart:rstart+rdim.
slong _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(nmod_poly_mat_t mat,
                                                     const slong * shift,
                                                     nmod_poly_mat_t tsf,
                                                     slong * det,
                                                     slong * pivind,
                                                     slong * rrp,
                                                     slong rstart,
                                                     slong cstart,
                                                     slong rdim,
                                                     slong cdim,
                                                     slong early_exit_zr,
                                                     orientation_t orient)
{
    if (rdim == 0 || cdim == 0)
        return 0;

    const slong zpiv = (orient == ROW_LOWER) ? -1 : cdim;

    // number of found pivots (i.e., rank), and number of found zero rows
    slong rk = 0;
    slong zr = 0;
    // All along, rdim == rk + zr + number of rows that remained to be processed
    // Also, we will ensure that pivind[:rk] remains correct (pivot index of
    // current submat[:rk,:]), and pivind[rdim-zr:rdim] is correct as well
    // (pivot index of the zr zero rows at the bottom of current submatrix)

    // -> pivot_row[j] is either -1 (not among the pivots discovered so far)
    // or is the index (in submatrix) of the row with pivot j
    slong * pivot_row = flint_malloc(cdim * sizeof(slong));
    flint_mpn_store(pivot_row, cdim, -1); // fill with -1

    slong pivdeg; // will be used to store pivot degrees

    while (rk + zr < rdim && zr < early_exit_zr)
    {
        // consider row rk of current matrix (corresponds to row rk+zr of initial matrix)
        // compute its pivot index
        _nmod_poly_vec_pivot_profile(pivind+rk, &pivdeg, mat->rows[rstart+rk]+cstart, shift, cdim, orient);
        if (pivind[rk] == zpiv)
        {
            // row is zero: rotate to put it last, increment zr and go to next row
            _nmod_poly_mat_rotate_rows_upward(mat, NULL, rstart+rk, rstart+rdim-1);
            if (tsf) _nmod_poly_mat_rotate_rows_upward(tsf, NULL, rstart+rk, rstart+rdim-1);
            if (det && (rdim-1-rk) % 2 == 1) *det = - *det;
            zr++;
        }
        else if (pivot_row[pivind[rk]] == -1)
        {
            // row provides new pivot: update pivot_row, increment rk and go to next row
            pivot_row[pivind[rk]] = rk;
            if (rrp) rrp[rk] = rk+zr;
            rk++;
        }
        else
        {
            // pivind[rk] == pivind[pi]: pivot collision, let's actually do some work
            slong pi = pivot_row[pivind[rk]];
            // see who has greatest pivot degree, and swap accordingly
            // (note that this does not disturb the row rank profile)
            if (pivdeg < nmod_poly_degree(MAT(rstart+pi, pivind[rk])))
            {
                nmod_poly_mat_swap_rows(mat, NULL, rstart+pi, rstart+rk);
                if (tsf) nmod_poly_mat_swap_rows(tsf, NULL, rstart+pi, rstart+rk);
                if (det && pi != rk) *det = - *det;
            }

            // perform atomic collision solving:
            //    mat[rstart+rk,:] = mat[rstart+rk,:] + cst * x**exp * mat[rstart+pi,:]
            _atomic_solve_pivot_collision_rowwise(mat, tsf, rstart+rk, rstart+pi, pivind[rk]);
        }
    }

    if (zr >= early_exit_zr)
        return -rk;
    else
        return rk;
}

slong nmod_poly_mat_ordered_weak_popov_iter(nmod_poly_mat_t mat,
                                            const slong * shift,
                                            nmod_poly_mat_t tsf,
                                            slong * pivind,
                                            slong * rrp,
                                            orientation_t orient)
{
    slong rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(mat, shift, tsf, NULL, pivind, rrp, 0, 0, mat->r, mat->c, mat->r, orient);

    slong * perm = flint_malloc(mat->r * sizeof(slong));
    _nmod_poly_mat_permute_rows_by_sorting_vec(mat, rk, pivind, perm);
    if (tsf)
        nmod_poly_mat_permute_rows(tsf, perm, NULL);
    flint_free(perm);
    return rk;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
