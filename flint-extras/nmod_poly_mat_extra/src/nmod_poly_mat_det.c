#include <flint/nmod_vec.h>
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_utils.h" // for permute_rows_by_sorting_vec
#include "nmod_poly_mat_forms.h"

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested
// this modifies mat, if all goes well (input matrix nonsingular) at the end it
// is upper triangular (up to permutation)
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    if (mat->r == 0) { nmod_poly_one(det); return; }

    // determinant (+1 or -1), pivot index, and rank for weak Popov form computations
    int udet = 1; 
    slong rk;
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    // permutation for putting into ordered weak Popov
    slong * perm = flint_malloc(mat->r * sizeof(slong));

    // window of original mat
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, mat, 0, 0, mat->r, mat->r);
    // TODO the following is necessary for Flint <= v3.0.0
    //view->modulus = mat->modulus; 

    for (slong i = mat->r -1; i >= 1; i--)
    {
        //            [ V  * ]                 [ V ]
        // mat is now [ 0  D ] and view is now [ 0 ] of size mat->r x (i+1),
        // with det = det(D) and V of size (i+1) x (i+1)
        //                 [ V'  *  * ]
        // -> transform to [ 0   d  * ] via weak Popov form of V[:,:i]
        //                 [ 0   0  D ]

        // this applies weak Popov form on view[:i+1,:i],
        // with transformations applied to the whole rows view[:i+1,:] = V
        // with early exit if detecting rank < i (in which case rk < 0)
        // with update of the determinant of unimodular transformation (+1 or -1)
        rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(view, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2, ROW_UPPER);

        // early exit if rank-deficient
        if (rk < i || nmod_poly_is_zero(nmod_poly_mat_entry(view, i, i)))
        {
            nmod_poly_zero(det);
            return;
        }

        // permute into ordered weak Popov form
        _nmod_poly_mat_permute_rows_by_sorting_vec(view, rk, pivind, perm);
        _nmod_poly_mat_window_resize_columns(view, -1);
        if (_perm_parity(perm, rk)) // odd permutation, negate udet
            udet = -udet;
    }
    flint_free(pivind);
    flint_free(perm);

    // retrieve determinant as product of diagonal entries
    // use view rather than mat, since mat has not been row-permuted
    // so it is only triangular up to permutation
    _nmod_poly_mat_window_resize_columns(view, mat->r -1); // reset view->c to mat->c
    nmod_poly_set(det, view->rows[0]+0); // recall here mat->r == mat->c > 0
    if (nmod_poly_is_zero(det))
        return; // rank deficient early exit, [0,0] had not been tested yet
    if (udet == -1)
        _nmod_vec_neg(det->coeffs, det->coeffs, det->length, det->mod);
    for (slong i = 1; i < view->r; i++)
        nmod_poly_mul(det, det, view->rows[i]+i);
    nmod_poly_mat_window_clear(view);
}

//void nmod_poly_mat_hnf_ur_(nmod_poly_t det, nmod_poly_mat_t mat)
//{
//    if (mat->r == 0) { nmod_poly_one(det); return; }
//
//    // determinant (+1 or -1), pivot index, and rank for weak Popov form computations
//    int udet = 1; 
//    slong rk;
//    slong * pivind = flint_malloc(mat->r * sizeof(slong));
//    // permutation for putting into ordered weak Popov
//    slong * perm = flint_malloc(mat->r * sizeof(slong));
//
//    for (slong i = mat->r -1; i >= 1; i--)
//    {
//        //            [ V  * ] 
//        // mat is now [ 0  D ] with det = det(D) and V of size (i+1) x (i+1)
//        //                 [ V'  *  * ]
//        // -> transform to [ 0   d  * ] via weak Popov form of V[:i+1,:i]
//        //                 [ 0   0  D ]
//
//        // this applies weak Popov form on mat[:i+1,:i],
//        // with transformations applied to the whole rows mat[:i+1,:]
//        // with early exit if detecting rank < i (in which case rk < 0)
//        // with update of the determinant of unimodular transformation (+1 or -1)
//        rk = _nmod_poly_mat_weak_popov_lr_iter_submat_rowbyrow(mat, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2);
//
//        // early exit if rank-deficient
//        if (rk < i || nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, i)))
//        {
//            nmod_poly_zero(det);
//            return;
//        }
//
//        // permute into ordered weak Popov form
//        _nmod_poly_mat_permute_rows_by_sorting_vec(mat, rk, pivind, perm);
//        if (_perm_parity(perm, rk)) // odd permutation, negate udet
//            udet = -udet;
//    }
//    flint_free(pivind);
//    flint_free(perm);
//
//    // retrieve determinant as product of diagonal entries
//    nmod_poly_set(det, mat->rows[0]+0); // recall here mat->r == mat->c > 0
//    if (nmod_poly_is_zero(det))
//        return; // rank deficient early exit, [0,0] had not been tested yet
//    if (udet == -1)
//        _nmod_vec_neg(det->coeffs, det->coeffs, det->length, det->mod);
//    for (slong i = 1; i < mat->r; i++)
//        nmod_poly_mul(det, det, mat->rows[i]+i);
//}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
