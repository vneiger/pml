/*
    Copyright (C) 2025 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_LINSOLVE_H
#define NMOD_POLY_MAT_LINSOLVE_H

/** \file nmod_poly_mat_linsolve.h
 *
 * Solving linear systems with univariate polynomial matrices, such as
 * computing the truncated power series solution of a linear system.
 *
 * \todo series solving at low and high precision, high-order lifting, and
 * solving via kernel bases (see ntl-extras, mat_lzz_pX_linsolve.h)
 */

#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * x-adic iterations à la Dixon : AX=B mod x^sigma
 * A is nxn, B is nxm 
 *  Returns 0 if A is not invertible as a power series (i.e., at x=0),
 *  in which case X is not modified
 *  Returns nonzero otherwise, and X is the sought solution
 * 
 * Iterations mod x^order to obtain a total approximation mod x^sigma; the
 * truncated inverse of A at order `order` is computed by
 * ::nmod_poly_mat_inv_trunc (see nmod_poly_mat_inverse.h)
 * 
 * Size of the integers (see @ref IntegerSizes): it is assumed that `order <=
 * WORD_MAX / 2` (see ::nmod_poly_mat_inv_trunc) and `sigma + order <=
 * WORD_MAX` (the iterations use multiples of `order` below `sigma + order`).
 */

slong nmod_poly_mat_dixon(nmod_poly_mat_t X, 
                          const nmod_poly_mat_t A, 
                          const nmod_poly_mat_t B, 
                          ulong order,
                          ulong sigma);

#ifdef __cplusplus
}
#endif

#endif  /* NMOD_POLY_MAT_LINSOLVE_H */
