#ifndef MAT_LZZ_P_EXTRA__H
#define MAT_LZZ_P_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector>

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0 0 c]                                                  */
/* [1 0 0 0]                                                  */
/* [0 1 0 0]                                                  */
/* [0 0 1 0]                                                  */
/* in size n                                                  */
/*------------------------------------------------------------*/
Mat<zz_p> Z_lzz_p(long n, const zz_p& c);

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0 0 1]                                                  */
/* [0 0 1 0]                                                  */
/* [0 1 0 0]                                                  */
/* [1 0 0 0]                                                  */
/* in size n                                                  */
/*------------------------------------------------------------*/
Mat<zz_p> J_lzz_p(long n);

#endif


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
