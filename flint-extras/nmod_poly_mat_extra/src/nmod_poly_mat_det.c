#include <flint/nmod_vec.h>
#include "nmod_poly_mat_forms.h"

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested
// this modifies mat, if all goes well (input matrix nonsingular) at the end it
// is upper triangular
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    // determinant (+1 or -1), pivod index, and rank for weak Popov form computations
    int udet; 
    slong rk;
    slong * pivind = flint_malloc(mat->r * sizeof(slong));

    nmod_poly_one(det);
    for (slong i = mat->r -1; i >= 0; i--)
    {
        //            [ S  * ]
        // mat is now [ 0  D ] with det = det(D) and S of size i+1 x i+1
        //                 [ S'  *  * ]
        // -> transform to [ 0   d  * ] by computing a weak Popov form of the first
        //                 [ 0   0  D ]
        //    i columns of S, storing the determinant d' of the transformation
        // -> det becomes d' * d * det

        // apply weak Popov on i+1 x i leading principal submatrix,
        // with transformations applied to the whole rows mat[0:i+1,:]
        // with early exit if detecting rank < i (in which case rk < 0)
        udet = 1;
        rk = _nmod_poly_mat_weak_popov_lr_iter_submat_rowbyrow(mat, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 1);
        if (rk < i)
        {
            nmod_poly_zero(det);
            return;
        }
        // multiply det by entry i+1,i+1  (called d above)
        if (udet == -1)
            _nmod_vec_neg(det->coeffs, det->coeffs, det->length, det->mod);
        nmod_poly_mul(det, det, mat->rows[i]+i);
    }
    flint_free(pivind);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
