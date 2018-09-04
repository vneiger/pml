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
void multiply_waksman(Mat<zz_pX> &Cout, const Mat<zz_pX> &A, const Mat<zz_pX> &B)
{
    Mat<zz_pX> C;
    Vec<zz_pX> Arow, Acol, D, E;
    zz_pX val0, val1, val2, crow;

    long m = A.NumRows();
    long n = A.NumCols();
    long p = B.NumCols();

    C.SetDims(m, p);

    D.SetLength(p);
    E.SetLength(m);

    long np = n>>1;

    for (long j=1; j<=np; j++)
    {
        const long j2=(j<<1)-1;

        for (long k=0; k<p; k++)
        {
            C[0][k] += (A[0][j2-1]+B[j2][k]) * (A[0][j2]+B[j2-1][k]);
            D[k] += (A[0][j2-1]-B[j2][k]) * (A[0][j2]-B[j2-1][k]);
        }

        for (long l=1; l<m; l++)
        {
            C[l][0] += (A[l][j2-1]+B[j2][0]) * (A[l][j2]+B[j2-1][0]);
            E[l] += (A[l][j2-1]-B[j2][0]) * (A[l][j2]-B[j2-1][0]);
        }

        for (long k=1; k<p; k++)
        {
            for (long l=1; l<m; l++)
            {
                C[l][k] += (A[l][j2-1]+B[j2][k]) * (A[l][j2]+B[j2-1][k]);
            }
        }
    }

    for (long l=1; l<m; l++)
    {
        E[l] = (E[l]+C[l][0])/2;
        C[l][0] -= E[l];
    }

    val0 = (D[0]+C[0][0])/2;
    C[0][0] -= val0;
    for (long k=1; k<p; k++)
    {
        val1 = (D[k]+C[0][k])/2;
        C[0][k] -= val1;
        val1 -= val0;
        for (long l=1; l<m; l++)
        {
            C[l][k] -= val1 + E[l];
        }
    }

    if ( (n&1) == 1)
    {
        for (long l=0; l<m; l++)
        {
            for (long k=0; k<p; k++)
            {
                C[l][k] += A[l][n-1]*B[n-1][k];
            }
        }
    }
    Cout = C;
}



/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* naive algorithm                                            */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_naive(Mat<zz_pX> & c_out, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{

    Mat<zz_pX> c;
    long u = a.NumRows();
    long v = a.NumCols();
    long w = b.NumCols();

    c.SetDims(u, w);

    for (long i = 0; i < u; i++)
    {
        for (long j = 0; j < w; j++)
        {
            c[i][j] = 0;
            for (long k = 0; k < v; k++)
            {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    c_out = c;
}


/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* naive algorithm                                            */
/* output may alias input; b does not have to be zero matrix  */
/*------------------------------------------------------------*/
void middle_product_naive(Mat<zz_pX> & b_out, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    Mat<zz_pX> b;
    long u = a.NumRows();
    long v = a.NumCols();
    long w = c.NumCols();

    b.SetDims(u, w);

    for (long i = 0; i < u; i++)
    {
        for (long j = 0; j < w; j++)
        {
            b[i][j] = 0;
            for (long k = 0; k < v; k++)
            {
                zz_pX tmp;
                middle_product(tmp, a[i][k], c[k][j], dA, dB);
                b[i][j] += tmp;
            }
        }
    }
    b_out = b;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
