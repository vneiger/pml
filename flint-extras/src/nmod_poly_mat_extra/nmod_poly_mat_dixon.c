/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <stdlib.h>
#include <math.h>

#include <flint/nmod_types.h>
#include <flint/nmod_poly.h>

#include "nmod_poly_mat_dixon.h"


/**
 *  TODO: probably only degree > 0 currently
 *
 *  Truncated inverse of A mod x^order
 *  A is assumed to be invertible for x=0 (checked by PML_ASSERT)
 *
 *  CHECK TODO: uses a method based on a shifted approximant computation
 *              check whether the specification of nmod_poly_mat_pmbasis
 *              is as expected
 *
 */
void nmod_poly_mat_inv_trunc(nmod_poly_mat_t S,
                            const nmod_poly_mat_t A,
                            ulong order)
{
    const slong n = A->r;

    // build H = [A \\ -Id]
    nmod_poly_mat_t H;
    nmod_poly_mat_init(H, 2*n, n, A->modulus);
    /* FIXME sets top n rows to A, but this is not really documented, so might change */
    nmod_poly_mat_set(H, A);

    for (slong i = n; i < 2*n; i++)
    {
        nmod_poly_fit_length(nmod_poly_mat_entry(H, i, i-n), 1);
        nmod_poly_mat_entry(H, i, i-n)->coeffs[0] = A->modulus - 1;
        _nmod_poly_set_length(nmod_poly_mat_entry(H, i, i-n), 1);
    }

    // Will be the approximant basis
    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, 2*n, 2*n, A->modulus);

    // Appropriate shift
    slong * shift = FLINT_ARRAY_ALLOC(2*n, slong);
    for (slong i = 0; i < n; i++)
        shift[i] = 0;
    for (slong i = n; i < 2*n; i++)
        shift[i] = order;

    // Shifted approximant computation
    nmod_poly_mat_pmbasis(tsf, shift, H, order);

    // Check whether the lower right part of the approximant is the identity matrix
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, tsf, n, n, 2*n, 2*n);
    PML_ASSERT(nmod_poly_mat_is_one(view));
    /* fails -> error in nmod_poly_mat_inv_trunc: check of identity, singular at origin? */
    nmod_poly_mat_window_clear(view);

    // Extract the truncated inverse from the approximant basis
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(S, i, j), nmod_poly_mat_entry(tsf, n+i, j));

    flint_free(shift);
    nmod_poly_mat_clear(H);
    nmod_poly_mat_clear(tsf);
}

/**
 *
 * TODO: probably only degree > 0 currently
 *
 * x-adic iterations à la Dixon : AX=B mod x^sigma
 * A is nxn, B is nxm
 * A is assumed to be invertible for x=0 (will produce an error otherwise)
 *
 * Iterations mod x^order to obtain a total approximation mod x^sigma
 *
 */
void nmod_poly_mat_dixon(nmod_poly_mat_t X,
                            const nmod_poly_mat_t A,
                            const nmod_poly_mat_t B,
                            ulong order,
                            ulong sigma)
{
    const slong m = B->c;
    const slong n = A->r;

    // C: inverse of A mod x^order
    nmod_poly_mat_t C;
    nmod_poly_mat_init(C, n, n, A->modulus);

    nmod_poly_mat_inv_trunc(C, A, order);

    // Temporary matrices
    nmod_poly_mat_t S;
    nmod_poly_mat_t T;
    nmod_poly_mat_t BB;
    nmod_poly_mat_init(S, n, m, A->modulus);
    nmod_poly_mat_init(T, n, m, A->modulus);
    nmod_poly_mat_init_set(BB, B);

    // The solution
    nmod_poly_mat_zero(X);

    // Main loop for Dixon's iterations
    for (slong l = 0; l * order < sigma; l++)
    {
        nmod_poly_mat_mul(S, C, BB);
        nmod_poly_mat_truncate(S, order);

        // X is constructed matrix digit by matrix digit
        nmod_poly_mat_shift_left(T, S, l*order);
        nmod_poly_mat_add(X, X, T);

        // New residue for next iteration
        nmod_poly_mat_mul(T, A, S);
        nmod_poly_mat_sub(T, BB, T);
        nmod_poly_mat_shift_right(BB, T, order);
    }

    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(BB);
    nmod_poly_mat_clear(S);
    nmod_poly_mat_clear(T);
}
