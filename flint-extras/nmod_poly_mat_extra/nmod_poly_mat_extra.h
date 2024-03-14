#ifndef NMOD_POLY_MAT_EXTRA_H
#define NMOD_POLY_MAT_EXTRA_H

/**
 * \file nmod_poly_mat_extra.h
 * \brief Main header for matrices with univariate polynomial entries modulo word-size prime
 * \version 0.0
 * \date 2022-06-25
 *
 * This is the main header for functions for matrices over the univariate
 * polynomials with coefficients in a finite field Z/pZ for "small" p
 * (word-size, so that we use a representation with flint's ``nmod``). The
 * purpose of this file is only to include all relevant headers (each one
 * gathers functions for a specific kind of tasks). This file contains general
 * TODOs, and may contain a few declarations that do not require a separate
 * header (for the moment).
 *
 * \todo benchmark performance
 * \todo test for memory leaks
 *
 * \todo Note: all parameters are supposed init
 *
 */

// include flint's matrices
#include <flint/nmod_poly_mat.h>

// include flint-extra's files
#include "nmod_poly_mat_utils.h"

#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_arith.h"

#include "nmod_poly_mat_multiply.h"

#include "nmod_poly_mat_approximant.h"
// #include "nmod_poly_mat_interpolant.h"

// #include "nmod_poly_mat_kernel.h"
// #include "nmod_poly_mat_determinant.h"

// #include "nmod_poly_mat_inverse.h"
// #include "nmod_poly_mat_linsolve.h"

// #include "nmod_poly_mat_linearization.h"


// TODO remove once using flint's comp instead
NMOD_POLY_MAT_INLINE void
apply_perm_to_vector(slong *res, const slong *initial_vect,
                          const slong *perm, slong length)
{
    for (slong i = 0; i < length; i++)
        res[perm[i]] = initial_vect[i];
}

/** \todo doc
 *  \todo move in suitable header
 **/
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat);

// TODO implem + doc
//slong nmod_poly_mat_linsolve_mulders_storjohann(nmod_poly_mat_t mat);

/** TODO: Other todos
 *  \todo algorithm for modular composition / charpoly mod
 *  \todo row reduction as in GJV (and via kernel?)
 *  \todo triangularization / Hermite
 */

#endif // NMOD_POLY_MAT_EXTRA_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
