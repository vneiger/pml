#ifndef NMOD_POLY_MAT_MULTIPLY_H
#define NMOD_POLY_MAT_MULTIPLY_H

#include <flint/nmod_types.h>

// several functions allocate arrays of matrices
// setting this flag allocates all memory at once
// slightly faster than using nmod_mat_init
#define DIRTY_ALLOC_MATRIX


/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: 2^(ceiling(log_2(lenA+lenB-1))) divides p-1 (assumption not checked)
 *  uses tft multiplication
 */
void nmod_poly_mat_mul_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: num columns of A < 2^30 and min(deg A, deg B) < 2^30 (assumption not checked)
 *  uses tft multiplication modulo 50 bits fft primes
 */
void nmod_poly_mat_mul_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: existence of primitive root (assumption not checked)
 *  uses geometric evaluation and interpolation
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



/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1), assuming deg(A) <= dA and deg(B) <= dA + dB
 *  output can alias input
 *  naive implementation (multiply, shift, truncate)
 */
void nmod_poly_mat_middle_product_naive(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                        const ulong dA, const ulong dB);

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  ASSUME: 2^(ceiling(log_2(dA + dB + 1))) divides p-1
 *  uses tft multiplication
 */
void nmod_poly_mat_middle_product_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                      const ulong dA, const ulong dB);

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  uses tft middle product modulo 50 bits fft primes
 */
void nmod_poly_mat_middle_product_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                      const ulong dA, const ulong dB);

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  ASSUME: existence of primitive root
 *  uses geometric evaluation and interpolation
 */
void nmod_poly_mat_middle_product_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                            const ulong dA, const ulong dB);


/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  ASSUME: large enough field 
 */
void nmod_poly_mat_middle_product_vdm1(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                       const ulong dA, const ulong dB);

#endif
