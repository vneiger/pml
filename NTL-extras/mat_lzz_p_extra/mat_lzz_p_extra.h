#ifndef MAT_LZZ_P_EXTRA__H
#define MAT_LZZ_P_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector>

NTL_CLIENT

Mat<zz_p> multiply_left_upper_toeplitz(const Mat<zz_p>& A, const Vec<zz_p> u);
Mat<zz_p> multiply_right_upper_toeplitz(const Mat<zz_p>& A, const Vec<zz_p> u);
Mat<zz_p> multiply_left_lower_toeplitz(const Mat<zz_p>& A, const Vec<zz_p> ell);
Mat<zz_p> multiply_right_lower_toeplitz(const Mat<zz_p>& A, const Vec<zz_p> ell);

#endif


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
