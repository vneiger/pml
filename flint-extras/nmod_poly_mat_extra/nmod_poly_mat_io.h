#ifndef NMOD_POLY_MAT_IO_H
#define NMOD_POLY_MAT_IO_H

/** \brief Input/output functions routines for univariate polynomial matrices
 * over `nmod`
 *
 * \file nmod_poly_mat_io.h
 * \author Vincent Neiger, Kevin Tran
 * \version 0.0
 * \date 2022-06-25
 *
 */

#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif


/** print for testing result on sage **/
void nmod_poly_mat_print_sage(const nmod_poly_mat_t mat, const char * var);

void slong_print_sage(const slong *shifts, slong length);

void nmod_mat_print_sage(const nmod_mat_t mat);

void slong_mat_print(const slong *mat, slong rdim, slong cdim);



#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_UTILS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
