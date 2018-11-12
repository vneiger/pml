#ifndef MAT_LZZ_PX_INVERSE__H
#define MAT_LZZ_PX_INVERSE__H

#include <memory> // for std::unique_ptr
#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                      INVERSE EXPANSION                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, quadratic algorithm               */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void plain_inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m);

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, Newton iteration                  */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/* crossover point with naive algo is m = 2^thresh            */
/* thresh = -1 means we use predetermined values              */
/*------------------------------------------------------------*/
void newton_inv_trunc_FFT(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh = -1);
void newton_inv_trunc_middle_product(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh = -1);
void newton_inv_trunc_geometric(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh = -1);

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m                                    */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m);

inline Mat<zz_pX> inv_trunc(const Mat<zz_pX>& a, long m)
{ Mat<zz_pX> y; inv_trunc(y, a, m); return y; }

/*------------------------------------------------------------*/
/* for i >= 0, define Si = coefficients of A^{-1} of degrees  */
/*             i-(2d-1) .. i-1, with d=deg(A)                 */
/* given src = Si, this computes S_{2i-d}                     */
/* invA = A^{-1} mod x^d                                      */
/* note: deg(Si) < 2d-1                                       */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void high_order_lift_inverse_odd(
                                 Mat<zz_pX> & next,
                                 const Mat<zz_pX>& src, 
                                 std::unique_ptr<mat_lzz_pX_lmultiplier> & A, 
                                 std::unique_ptr<mat_lzz_pX_lmultiplier> & invA,
                                 long d
                                );


#endif /* end of include guard: MAT_LZZ_PX_INVERSE__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
