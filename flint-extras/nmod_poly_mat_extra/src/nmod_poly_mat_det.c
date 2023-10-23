#include <flint/nmod_vec.h>
#include "nmod_poly_mat_forms.h"

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested
// this modifies mat, if all goes well (input matrix nonsingular) at the end it
// is upper triangular
// TODO early exit when meeting given det bound? does it speed up generic case?
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    // determinant (+1 or -1), pivot index, and rank for weak Popov form computations
    int udet = 1; 
    slong rk;
    slong * pivind = flint_malloc(mat->r * sizeof(slong));

    for (slong i = mat->r -1; i >= 0; i--)
    {
        //            [ S  * ]
        // mat is now [ 0  D ] with det = det(D) and S of size i+1 x i+1
        //                 [ S'  *  * ]
        // -> transform to [ 0   d  * ] via weak Popov form of S[:,:i]
        //                 [ 0   0  D ]

        // apply weak Popov on i+1 x i leading principal submatrix,
        // with transformations applied to the whole rows mat[0:i+1,:]
        // with early exit if detecting rank < i (in which case rk < 0)
        // with update of the determinant of unimodular transformation (+1 or -1)
        rk = _nmod_poly_mat_weak_popov_lr_iter_submat_rowbyrow(mat, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2);
        if (rk < i)
        {
            nmod_poly_zero(det);
            return;
        }
    }
    flint_free(pivind);

    // retrieve determinant as product of diagonal entries
    nmod_poly_one(det);
    if (udet == -1)
        det->coeffs[0] = nmod_neg(det->coeffs[0], det->mod);
    for (slong i = 0; i < mat->r; i++)
        nmod_poly_mul(det, det, mat->rows[i]+i);
    return;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
