#ifndef NMOD_POLY_MAT_MULTIPLY_H
#define NMOD_POLY_MAT_MULTIPLY_H

#include <flint/nmod_poly_mat.h>
#include "nmod_poly_extra.h"

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses tft multiplication
 *  ASSUME: 2^(ceiling(log_2(lenA+lenB-1))) divides p-1
 */
void nmod_poly_mat_mul_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses tft multiplication modulo 50 bits fft primes
 */
void nmod_poly_mat_mul_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses geometric evaluation and interpolation
 *  ASSUME: existence of primitive root
 */
void nmod_poly_mat_mul_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  ASSUME: large enough field 
 */
void nmod_poly_mat_mul_vdm1(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);
void nmod_poly_mat_mul_vdm2(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses waksman's algorithm
 *  ASSUME: p != 2
 */
void nmod_poly_mat_mul_waksman(nmod_poly_mat_t C, const nmod_poly_mat_t A,  const nmod_poly_mat_t B);

#endif
