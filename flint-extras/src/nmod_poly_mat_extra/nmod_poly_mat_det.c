/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"  /* for permute_rows_by_sorting_vec */
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_extra.h"  /* for prototype of det_iter */

// Mulders-Storjohann triangularization, see the documentation in
// nmod_poly_mat_extra.h
slong _nmod_poly_mat_upper_triangularize_iter(nmod_poly_struct * diag, nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    if (n == 0)
        return 1;

    // determinant (+1 or -1) of the unimodular transformation, 0 if singular
    slong udet = 1;
    slong * pivind = flint_malloc(n * sizeof(slong));
    // permutation for putting into ordered weak Popov
    slong * perm = flint_malloc(n * sizeof(slong));

    // window of original mat
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, mat, 0, 0, n, n);

    for (slong i = n - 1; i >= 1; i--)
    {
        //            [ V  * ]                 [ V ]
        // mat is now [ 0  D ] and view is now [ 0 ] of size n x (i+1),
        // with V of size (i+1) x (i+1)
        //                 [ V'  *  * ]
        // -> transform to [ 0   d  * ] via weak Popov form of V[:,:i]
        //                 [ 0   0  D ]
        // (weak Popov form of view[:i+1,:i], with transformations applied to
        // the whole rows view[:i+1,:] = V, early exit if rank < i)
        const slong rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(view, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2, ROW_UPPER);

        if (rk < i || nmod_poly_is_zero(nmod_poly_mat_entry(view, i, i)))
        {
            udet = 0;
            break;
        }

        // permute into ordered weak Popov form
        _nmod_poly_mat_permute_rows_by_sorting_vec(view, rk, pivind, perm);
        _nmod_poly_mat_window_resize_columns(view, -1);
        if (_perm_parity(perm, rk)) // odd permutation, negate udet
            udet = -udet;
    }

    // read the diagonal through view rather than mat, since with FLINT < 3.3
    // the rows of mat have not been permuted
    _nmod_poly_mat_window_resize_columns(view, n - view->c);
    if (udet != 0 && nmod_poly_is_zero(nmod_poly_mat_entry(view, 0, 0)))
        udet = 0;
    if (udet != 0)
        for (slong i = 0; i < n; i++)
            nmod_poly_swap(diag + i, nmod_poly_mat_entry(view, i, i));

    flint_free(pivind);
    flint_free(perm);
    nmod_poly_mat_window_clear(view);
    return udet;
}

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested; this modifies mat (see
// _nmod_poly_mat_upper_triangularize_iter)
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    nmod_poly_struct * diag = flint_malloc(FLINT_MAX(1, n) * sizeof(nmod_poly_struct));
    for (slong i = 0; i < n; i++)
        nmod_poly_init(diag + i, mat->modulus);

    const slong udet = _nmod_poly_mat_upper_triangularize_iter(diag, mat);
    if (udet == 0)
        nmod_poly_zero(det);
    else if (n == 0)
        nmod_poly_one(det);
    else
    {
        /* determinant = udet * product of the diagonal entries */
        nmod_poly_swap(det, diag + 0);
        for (slong i = 1; i < n; i++)
            nmod_poly_mul(det, det, diag + i);
        if (udet == -1)
            nmod_poly_neg(det, det);
    }

    for (slong i = 0; i < n; i++)
        nmod_poly_clear(diag + i);
    flint_free(diag);
}
