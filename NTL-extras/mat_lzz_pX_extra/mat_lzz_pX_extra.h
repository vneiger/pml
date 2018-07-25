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
/* matrix of integers; deg(0) = -1                            */
/*------------------------------------------------------------*/
void degree_matrix(Mat<long> &a, const Mat<zz_pX> &b, const Vec<long>& shift=Vec<long>(), const bool row_wise=true);

/*------------------------------------------------------------*/
/* max degree of row entries                                  */
/*------------------------------------------------------------*/
void row_degree(Vec<long> &a, const Mat<zz_pX> &b, const Vec<long>& shift=Vec<long>()); 

/*------------------------------------------------------------*/
/* max degree of col entries                                  */
/*------------------------------------------------------------*/
void col_degree(Vec<long> &a, const Mat<zz_pX> &b,const Vec<long>& shift=Vec<long>()); 

/*------------------------------------------------------------*/
/* leading matrix of b                                        */
/*------------------------------------------------------------*/
void leading_matrix(Mat<zz_p> &a, const Mat<zz_pX> &b,const Vec<long> & shift = Vec<long>(), const bool row_wise = true);

/*------------------------------------------------------------*/
/* returns true if b is reduced                               */
/*------------------------------------------------------------*/
bool is_reduced (const Mat<zz_pX> &b,const Vec<long> & shift = Vec<long>(), const bool row_wise = true);

/*------------------------------------------------------------*/
/* finds the pivot indices; returns the row/col degs          */
/*------------------------------------------------------------*/
Vec<long> pivot_index (Vec<long> &index, const Mat<zz_pX> &b,const Vec<long> & shift = Vec<long>(), const bool row_wise = true);

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

#endif






































