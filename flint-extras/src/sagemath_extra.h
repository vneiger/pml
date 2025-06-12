/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef SAGEMATH_EXTRA_H
#define SAGEMATH_EXTRA_H

/** \brief Printing towards SageMath format
 *
 * \file sagemath_extras.h
 * \version 0.0
 * \date 2022-06-25
 *
 */

// Kevin Tran: 2022-06: initial version
// Vincent Neiger: 2022-07: cleaned, added nmod_poly_print

#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

// TODO add printing to file?

/** printing a vector of slong */
void slongvec_print_sagemath(const slong * shift, slong length);

/** printing an `nmod_poly` */
static inline void
nmod_poly_print_sagemath(const nmod_poly_t pol, const char * var)
{
    nmod_poly_print_pretty(pol, var);
}

/** printing an `nmod_mat` */
void nmod_mat_print_sagemath(const nmod_mat_t mat);

/** printing an `nmod_poly_mat` **/
void nmod_poly_mat_print_sagemath(const nmod_poly_mat_t mat, const char * var);

#ifdef __cplusplus
}
#endif

#endif // SAGEMATH_EXTRA_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
