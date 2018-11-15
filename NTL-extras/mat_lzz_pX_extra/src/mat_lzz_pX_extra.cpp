#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota
#include <NTL/BasicThreadPool.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* UTILS                                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* horizontal join                                            */
/* requires a.NumRows() == b.NumRows()                        */
/*------------------------------------------------------------*/
void horizontal_join(Mat<zz_pX>& c, const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{
    long r = a.NumRows();
    if (r != b.NumRows())
        LogicError("Dimension mismatch for horizontal join");
    long ca = a.NumCols();
    long cb = b.NumCols();

    Mat<zz_pX> c2;
    c2.SetDims(r, ca + cb);
    for (long i = 0; i < r; i++)
        for (long j = 0; j < ca; j++)
            c2[i][j] = a[i][j];
    for (long i = 0; i < r; i++)
        for (long j = 0; j < cb; j++)
            c2[i][j+ca] = b[i][j];
    c = c2;
}

/*------------------------------------------------------------*/
/* collapses s consecutive columns of a into one column of c  */
/* let t=a.NumCols(). For i=0..t/s-1, the i-th column of c is */
/* a[i*s] + x^d a[i*s+1] + ... + x^{(s-1)*d} a[i*s+s-1)]      */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_consecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s)
{
    if (&c == &a)
    {
        c = collapse_consecutive_columns(a, d, s);
        return;
    }
    long rho = a.NumRows();
    long t = a.NumCols();
    long r = t / s;
    if (t != r * s)
        LogicError("Bad parameters for collapse_consecutive_columns");

    c.SetDims(rho, r);
    for (long j = 0; j < rho; j++)
        for (long i = 0; i < r; i++)
        {
            c[j][i] = 0;
            c[j][i].SetLength(d * s);
            for (long k = 0; k < s; k++)
                for (long ell = 0; ell < d; ell++)
                    SetCoeff(c[j][i], k*d + ell, coeff(a[j][i*s + k], ell));
        }
}

/*------------------------------------------------------------*/
/* collapses columns with stepsize s of a into a column of c  */
/* let t=a.NumCols(). For i=0..s-1, the i-th column of c is   */
/* a[i] + x^d a[i+s] + ... + x^{(t/s-1)*d} a[i+(t/s-1)*s)]    */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_nonconsecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s)
{
    if (&c == &a)
    {
        c = collapse_nonconsecutive_columns(a, d, s);
        return;
    }
    long rho = a.NumRows();
    long t = a.NumCols();
    long r = t / s;
    if (t != r * s)
        LogicError("Bad parameters for collapse_nonconsecutive_columns");

    c.SetDims(rho, s);
    for (long j = 0; j < rho; j++)
        for (long i = 0; i < s; i++)
        {
            c[j][i] = 0;
            c[j][i].SetLength(d * r);
            for (long k = 0; k < r; k++)
                for (long ell = 0; ell < d; ell++)
                    SetCoeff(c[j][i], k*d + ell, coeff(a[j][i + k*s], ell));
        }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
