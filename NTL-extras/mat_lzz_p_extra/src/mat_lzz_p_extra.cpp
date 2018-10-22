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
/* [a0 a2 a1 a0 a2]                                           */
/* [a1 a0 a2 a1 a0]                                           */
/* [a2 a1 a0 a2 a1]                                           */
/* in size nrows, ncols. A must have length nrows             */
/*------------------------------------------------------------*/
Mat<zz_p> circulant_lzz_p(const Vec<zz_p>& A, long nrows, long ncols)
{
    if (A.length() != nrows)
        Error("Bad vector size for circulant_lzz_p");

    Mat<zz_p> res;
    res.SetDims(nrows, ncols);
    for (long j = 0; j < ncols; j++)
    {
        long idx = j % nrows;
        for (long i = 0; i < nrows; i++)
        {
            res[idx][j] = A[i];
            idx++;
            if (idx == nrows)
                idx = 0;
        }
    }
    return res;
}

/*------------------------------------------------------------*/
/* returns                                                    */
/* [a0  0  0 0]                                               */
/* [a1 a0  0 0]                                               */
/* [a2 a1 a0 0]                                               */
/* in size nrows, ncols. A must have length nrows             */
/*------------------------------------------------------------*/
Mat<zz_p> lower_toeplitz_lzz_p(const Vec<zz_p>& A, long nrows, long ncols)
{
    if (A.length() != nrows)
        Error("Bad vector size for circulant_lzz_p");

    Mat<zz_p> res;
    res.SetDims(nrows, ncols);
    for (long j = 0; j < min(nrows, ncols); j++)
    {
        for (long i = j; i < nrows; i++)
            res[i][j] = A[i - j];
    }
    return res;
}

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0  0  a0]                                               */
/* [0 0  a0 a1]                                               */
/* [0 a0 a1 a2]                                               */
/* in size nrows, ncols. A must have length nrows             */
/*------------------------------------------------------------*/
Mat<zz_p> lower_hankel_lzz_p(const Vec<zz_p>& A, long nrows, long ncols)
{
    if (A.length() != nrows)
        Error("Bad vector size for circulant_lzz_p");

    Mat<zz_p> res;
    res.SetDims(nrows, ncols);
    for (long j = 0; j < min(nrows, ncols); j++)
    {
        for (long i = j; i < nrows; i++)
            res[i][ncols -1 -j] = A[i - j];
    }
    return res;
}

/*------------------------------------------------------------*/
/* returns                                                    */
/* [a2 a1 a0 0]                                               */
/* [a1 a0  0 0]                                               */
/* [a0  0  0 0]                                               */
/* in size nrows, ncols. A must have length nrows             */
/*------------------------------------------------------------*/
Mat<zz_p> upper_toeplitz_lzz_p(const Vec<zz_p>& A, long nrows, long ncols)
{
    if (A.length() != nrows)
        Error("Bad vector size for circulant_lzz_p");

    Mat<zz_p> res;
    res.SetDims(nrows, ncols);
    for (long j = 0; j < min(nrows, ncols); j++)
    {
        for (long i = j; i < nrows; i++)
            res[nrows -1 -i][j] = A[i - j];
    }
    return res;
}

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0  0  a0]                                               */
/* [0 0  a0 a1]                                               */
/* [0 a0 a1 a2]                                               */
/* in size nrows, ncols. A must have length nrows             */
/*------------------------------------------------------------*/
Mat<zz_p> upper_hankel_lzz_p(const Vec<zz_p>& A, long nrows, long ncols)
{
    if (A.length() != nrows)
        Error("Bad vector size for circulant_lzz_p");

    Mat<zz_p> res;
    res.SetDims(nrows, ncols);
    for (long j = 0; j < min(nrows, ncols); j++)
    {
        for (long i = j; i < nrows; i++)
            res[nrows -1 -i][ncols -1 -j] = A[i - j];
    }
    return res;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

