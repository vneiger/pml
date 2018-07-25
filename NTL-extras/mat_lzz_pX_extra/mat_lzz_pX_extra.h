#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* random matrix of a given degree                            */
/*------------------------------------------------------------*/
void random_mat_zz_pX(Mat<zz_pX>& a, long n, long m, long d);

/*------------------------------------------------------------*/
/* maximum degree of the entries of a                         */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & a);

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

#endif
