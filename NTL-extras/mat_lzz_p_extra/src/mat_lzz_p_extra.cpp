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
    Mat<zz_p> res(INIT_SIZE, n, n);
    for (long i = 1; i < n; ++i)
        set(res[i][i-1]);
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
    Mat<zz_p> res(INIT_SIZE, n, n);
    res.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        set(res[i][n-1-i]);
    return res;
}

/*------------------------------------------------------------*/
/* diagonal matrix with diagonal d                            */
/*------------------------------------------------------------*/
Mat<zz_p> diagonal_matrix(const Vec<zz_p> & d)
{
    const long n = d.length();
    Mat<zz_p> res(INIT_SIZE, n, n);
    for (long i = 0; i < n; ++i)
        res[i][i] = d[i];
    return res;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
