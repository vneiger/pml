#ifndef __FMPZ_MAT_EXTRA_H
#define __FMPZ_MAT_EXTRA_H

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "pml.h"

/** ------------------------------------------------------------ */
/** a multimodular algorithm for matrix multiplication           */
/** based on flint's fmpz_mat_mul_multi_mod implementation       */
/** uses our fmpz_multimod_CRT nmod_mat_mul_small_modulus        */ 
/** ------------------------------------------------------------ */
void fmpz_mat_mul_multimod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

#endif
