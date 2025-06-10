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

#include <nmod_poly_mat_approximant.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 *
 *  Truncated inverse of A mod x^order
 *  A is assumed to be invertible for x=0 (will produce an error otherwise) 
 * 
 *  CHECK TODO: uses a method based on a shifted approximant computation 
 *              check whether the specification of nmod_poly_mat_pmbasis 
 *              is as expected 
 * 
 */

void nmod_poly_mat_inv_trunc(nmod_poly_mat_t S, 
                            const nmod_poly_mat_t A, 
                            ulong order);


/**
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
                            ulong sigma);

#ifdef __cplusplus
}
#endif

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


