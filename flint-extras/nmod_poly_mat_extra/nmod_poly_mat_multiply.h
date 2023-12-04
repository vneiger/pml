#ifndef NMOD_POLY_MAT_MULTIPLY_H
#define NMOD_POLY_MAT_MULTIPLY_H

#include <flint/nmod_poly_mat.h>

#include "nmod_poly_extra.h"

// assume 2^(ceiling(log_2(lenA+lenB-1))) divides p-1
void nmod_poly_mat_mul_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

void nmod_poly_mat_mul_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

#endif
