/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_types.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_inverse.h"
#include "nmod_poly_mat_linsolve.h"
#include "nmod_poly_mat_multiply.h"


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
slong nmod_poly_mat_dixon(nmod_poly_mat_t X,
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

    slong is_inv = nmod_poly_mat_inv_trunc(C, A, order);
    if (!is_inv)
        return 0;

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
        nmod_poly_mat_multiply(S, C, BB);
        nmod_poly_mat_truncate(S, order);

        // X is constructed matrix digit by matrix digit
        nmod_poly_mat_shift_left(T, S, l*order);
        nmod_poly_mat_add(X, X, T);

        // New residue for next iteration
        nmod_poly_mat_multiply(T, A, S);
        nmod_poly_mat_sub(T, BB, T);
        nmod_poly_mat_shift_right(BB, T, order);
    }

    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(BB);
    nmod_poly_mat_clear(S);
    nmod_poly_mat_clear(T);
    return 1;
}
