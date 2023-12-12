#ifndef __FMPZ_POLY_MAT_EXTRA_H
#define __FMPZ_POLY_MAT_EXTRA_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly_mat.h>

/** ------------------------------------------------------------ */
/** Waksman's algorithm for matrix multiplication                */
/** does n^3/2+O(n^2) products, but many additions               */
/** good for small matrices with large entries                   */
/** ------------------------------------------------------------ */
void fmpz_poly_mat_mul_waksman(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);


/** ------------------------------------------------------------ */
/** a multimodular algorithm for matrix multiplication           */
/** ------------------------------------------------------------ */
void fmpz_poly_mat_mul_multimod(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

#endif
