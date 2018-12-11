#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* Waksman's algorithm                                        */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> &c, const Mat<zz_pX> &a, const Mat<zz_pX> &b)
{
    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_waksman(c2, a, b);
        c.swap(c2);
        return;
    }

    Vec<zz_pX> d, e;
    zz_pX val0, val1;

    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    c.SetDims(m, p);
    clear(c);

    d.SetLength(p);
    e.SetLength(m);

    for (long j=1; j<=n/2; ++j)
    {
        const long j2=2*j-1;

        for (long k=0; k<p; k++)
        {
            c[0][k] += (a[0][j2-1]+b[j2][k]) * (a[0][j2]+b[j2-1][k]);
            d[k] += (a[0][j2-1]-b[j2][k]) * (a[0][j2]-b[j2-1][k]);
        }

        for (long l=1; l<m; l++)
        {
            c[l][0] += (a[l][j2-1]+b[j2][0]) * (a[l][j2]+b[j2-1][0]);
            e[l] += (a[l][j2-1]-b[j2][0]) * (a[l][j2]-b[j2-1][0]);
        }

        for (long k=1; k<p; k++)
        {
            for (long l=1; l<m; l++)
            {
                c[l][k] += (a[l][j2-1]+b[j2][k]) * (a[l][j2]+b[j2-1][k]);
            }
        }
    }

    for (long l=1; l<m; l++)
    {
        e[l] = (e[l]+c[l][0])/2;
        c[l][0] -= e[l];
    }

    val0 = (d[0]+c[0][0])/2;
    c[0][0] -= val0;
    for (long k=1; k<p; k++)
    {
        val1 = (d[k]+c[0][k])/2;
        c[0][k] -= val1;
        val1 -= val0;
        for (long l=1; l<m; l++)
        {
            c[l][k] -= val1 + e[l];
        }
    }

    if ( (n&1) == 1)
        for (long l=0; l<m; l++)
            for (long k=0; k<p; k++)
                c[l][k] += a[l][n-1]*b[n-1][k];
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* naive algorithm                                            */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply_naive(c2, a, b);
        c.swap(c2);
        return;
    }

    const long u = a.NumRows();
    const long v = a.NumCols();
    const long w = b.NumCols();

    c.SetDims(u, w);
    zz_pX buf;
    for (long i = 0; i < u; ++i)
        for (long j = 0; j < w; ++j)
        {
            mul(c[i][j], a[i][0], b[0][j]);
            for (long k = 1; k < v; ++k)
            {
                mul(buf, a[i][k], b[k][j]);
                add(c[i][j], c[i][j], buf);
            }
        }
}


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* naive algorithm                                            */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void middle_product_naive(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    if (&b == &a || &b == &c)
    {
        Mat<zz_pX> b2;
        middle_product_naive(b2, a, c, dA, dB);
        b.swap(b2);
        return;
    }

    long u = a.NumRows();
    long v = a.NumCols();
    long w = c.NumCols();

    b.SetDims(u, w);
    zz_pX buf;
    for (long i = 0; i < u; ++i)
        for (long j = 0; j < w; ++j)
        {
            middle_product(b[i][j], a[i][0], c[0][j], dA, dB);
            for (long k = 1; k < v; ++k)
            {
                middle_product(buf, a[i][k], c[k][j], dA, dB);
                add(b[i][j], b[i][j], buf);
            }
        }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
