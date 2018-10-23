#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "mat_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0 0 c]                                                  */
/* [1 0 0 0]                                                  */
/* [0 1 0 0]                                                  */
/* [0 0 1 0]                                                  */
/* in size n                                                  */
/*------------------------------------------------------------*/
Mat<zz_p> Z_lzz_p(long n, const zz_p& c)
{
    Mat<zz_p> res;
    res.SetDims(n, n);
    for (long i = 1; i < n; i++)
        res[i][i-1] = 1;
    res[0][n-1] = c;
    return res;
}

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0 0 1]                                                  */
/* [0 0 1 0]                                                  */
/* [0 1 0 0]                                                  */
/* [1 0 0 0]                                                  */
/* in size n                                                  */
/*------------------------------------------------------------*/
Mat<zz_p> J_lzz_p(long n)
{
    Mat<zz_p> res;
    res.SetDims(n, n);
    for (long i = 0; i < n; i++)
        res[i][n - 1 - i] = 1;
    return res;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

