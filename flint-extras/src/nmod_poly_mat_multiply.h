/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_MULTIPLY_H
#define NMOD_POLY_MAT_MULTIPLY_H

#include <flint/nmod_types.h>
#include <flint/nmod_poly.h>

#include "pml.h"

// temporarily disabling FFT_SMALL based variants
#define FFT_SMALL_VARIANTS 0

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: 2^(ceiling(log_2(lenA+lenB-1))) divides p-1 (assumption not checked)
 *  uses tft multiplication
 *  \todo temporarily disabled
 */
#if FFT_SMALL_VARIANTS
void nmod_poly_mat_mul_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);
#endif

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: num columns of A < 2^30 and min(deg A, deg B) < 2^30 (assumption not checked)
 *  uses tft multiplication modulo 50 bits fft primes
 *  \todo temporarily disabled
 */
#if FFT_SMALL_VARIANTS
void nmod_poly_mat_mul_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);
#endif

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  v1 requires that the field contains distinct points 0, 1, 2, ..., max_length(A) + max_length(B) - 2
 *     (for Z/pZ, equivalent to cardinality >= max_length(A) + max_length(B) - 1)
 *  v2 requires that the field contains distinct points 1**2, 2**2, ..., (max_length(A) + max_length(B) - 1)**2
 *     (for Z/pZ, implied by cardinality > (max_length(A) + max_length(B) - 1)**2)
 */
void nmod_poly_mat_mul_vandermonde1(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);
void nmod_poly_mat_mul_vandermonde2(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses waksman's algorithm
 *  ASSUME: p != 2
 */
void nmod_poly_mat_mul_waksman(nmod_poly_mat_t C, const nmod_poly_mat_t A,  const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices, using evaluation-interpolation at geometric progression
 *  output can alias input
 *  in _precomp, G is precomputed to support evaluation and interpolation with G->len >= len1+len2-1
 *  the other variant computes G, assuming this is feasible (not checked)
 */
void nmod_poly_mat_mul_geometric(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2);
void _nmod_poly_mat_mul_geometric_precomp(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G);


/** general interface, picks an algorithm depending on parameters
 * TODO thresholds to be tuned
 */
void nmod_poly_mat_multiply(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2);



/** Middle product for polynomial matrices
 * For all functions below:
 * They compute
 *       res  <-  (pmat1 * pmat2 mod x^nhi) div x^nlo
 * i.e., they set res to the first nhi - nlo middle coefficients [nlo, nhi) of
 * the product pmat1 * pmat2
 * Input len1 and len2 is such that pmat1 has length <= len1 and pmat2 has
 * length <= len2
 * Output can alias input
 *
 * TODO computation flow currently not ideal: reductions (e.g., if pmat1 has length << nlo)
 * are done in the general interface, and so will not be done in case one wants to directly
 * call the geometric variant
 *
 * TODO for evaluation-interpolation (currently, geometric), would be good to
 * have separate functions for matrix evaluation and matrix interpolation, and just call them?
 * (but then we also need it for "matrix evaluation of reverse"...)
 *
 * TODO faster matrix eval/interp by better handling the memory aspects

 * TODO for naive, add and use function for combined shift+truncate
 *
 * TODO improve the case of constant pmat1 or pmat2, currently no specific function
 */

/** general interface
 * tries a few reductions and calls _nmod_poly_mat_mulmid */
void nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                          const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                          slong nlo, slong nhi);

/** main interface
 * picks the fastest variant depending on parameters */
void _nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                           const nmod_poly_mat_t pmat1, slong len1,
                           const nmod_poly_mat_t pmat2, slong len2,
                           slong nlo, slong nhi);

/** actual worker: naive method, multiply + shift+truncate */
void _nmod_poly_mat_mulmid_naive(nmod_poly_mat_t res,
                                 const nmod_poly_mat_t pmat1, slong len1,
                                 const nmod_poly_mat_t pmat2, slong len2,
                                 slong nlo, slong nhi);

#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
/** actual worker: transposed multiplication, via evaluation-interpolation at geometric progression
 * - requires len1 <= nlo+1 or len2 <= nlo+1
 * - requires the existence of a geometric progression of order (at least) nhi
 * - finds such a progression, builds temporary precomputation data,
 *   and calls `mulmid_geometric{1,2}_precomp`
 */
void _nmod_poly_mat_mulmid_geometric(nmod_poly_mat_t res,
                                     const nmod_poly_mat_t pmat1, slong len1,
                                     const nmod_poly_mat_t pmat2, slong len2,
                                     slong nlo, slong nhi);

/** same as above, with additional G, precomputed data for evaluation and
 * interpolation at a geometric progression of order at least nhi
 * requires len1 <= nlo+1 or len2 <= nlo+1
 */
void _nmod_poly_mat_mulmid_geometric_precomp(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi, nmod_geometric_progression_t G);
#endif


/*------------------------------------------------------------*/
/* TODO currently disabled variants                           */
/*------------------------------------------------------------*/

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  ASSUME: 2^(ceiling(log_2(dA + dB + 1))) divides p-1
 *  uses tft multiplication
 */
#if FFT_SMALL_VARIANTS
void nmod_poly_mat_middle_product_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                      const ulong dA, const ulong dB);
#endif

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  uses tft middle product modulo 50 bits fft primes
 */
#if FFT_SMALL_VARIANTS
void nmod_poly_mat_middle_product_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                      const ulong dA, const ulong dB);
#endif


/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  ASSUME: large enough field 
 */
void nmod_poly_mat_middle_product_vdm1(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                       const ulong dA, const ulong dB);

#endif
