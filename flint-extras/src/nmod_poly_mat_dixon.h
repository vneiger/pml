/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_DIXON_H
#define NMOD_POLY_MAT_DIXON_H

#include "nmod_poly_mat_approximant.h"

#ifdef __cplusplus
extern "C" {
#endif


/**
 *
 *  Truncated inverse of A mod x^order
 *  Returns 0 if A is not invertible as a power series (i.e., at x=0),
 *  in which case S is not modified
 *  Returns nonzero otherwise, and S is the truncated inverse
 * 
 *  CHECK TODO: uses a method based on a shifted approximant computation 
 *              check whether the specification of nmod_poly_mat_pmbasis 
 *              is as expected 
 * 
 */

slong nmod_poly_mat_inv_trunc(nmod_poly_mat_t S, 
                              const nmod_poly_mat_t A, 
                              ulong order);


/**
 * x-adic iterations à la Dixon : AX=B mod x^sigma
 * A is nxn, B is nxm 
 *  Returns 0 if A is not invertible as a power series (i.e., at x=0),
 *  in which case X is not modified
 *  Returns nonzero otherwise, and X is the sought solution
 * 
 * Iterations mod x^order to obtain a total approximation mod x^sigma
 * 
 */

slong nmod_poly_mat_dixon(nmod_poly_mat_t X, 
                          const nmod_poly_mat_t A, 
                          const nmod_poly_mat_t B, 
                          ulong order,
                          ulong sigma);

#ifdef __cplusplus
}
#endif

#endif  /* NMOD_POLY_MAT_DIXON_H */
