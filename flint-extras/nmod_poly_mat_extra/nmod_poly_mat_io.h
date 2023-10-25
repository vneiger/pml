#ifndef NMOD_POLY_MAT_IO_H
#define NMOD_POLY_MAT_IO_H

/** \brief Input/output functions routines for univariate polynomial matrices
 * over `nmod`
 *
 * \file nmod_poly_mat_io.h
 * \version 0.0
 * \date 2022-06-25
 *
 */

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_poly_mat_forms.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \todo print to file
 *  \todo pretty print to file
 */

/** Matrix pretty print to standard output */
void nmod_poly_mat_print_pretty(const nmod_poly_mat_t mat, const char * var);


/** Print the degree matrix, see @ref DegreeMatrix */
void nmod_poly_mat_degree_matrix_print_pretty(const nmod_poly_mat_t mat);
void nmod_poly_mat_degree_matrix_shifted_print_pretty(const nmod_poly_mat_t mat,
                                                      const slong *shift,
                                                      orientation_t orient);

/** Print the leading matrix, see @ref LeadingMatrix . For uniform shift, one
 * can input `shift = NULL`. */
void nmod_poly_mat_leading_matrix_print_pretty(const nmod_poly_mat_t mat,
                                               const slong *shift,
                                               orientation_t orient);

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_UTILS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
