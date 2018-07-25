#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* random matrix of a given degree                            */
/*------------------------------------------------------------*/
void random_mat_zz_pX(Mat<zz_pX>& a, long n, long m, long d)
{
    a.SetDims(n, m);
    for (long i = 0; i < n; i++)
    {
	for (long j = 0; j < m; j++)
	{
	    a[i][j] = random_zz_pX(d);
	}
    }
}

/*------------------------------------------------------------*/
/* maximum degree of the entries of a                         */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & a)
{
    long d = -1;
    for (long i = 0; i < a.NumRows(); i++)
    {
	for (long j = 0; j < a.NumCols(); j++)
	{
	    d = max(d, deg(a[i][j]));
	}
    }
    return d;
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* Waksman's algorithm                                        */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> &C, const Mat<zz_pX> &A, const Mat<zz_pX> &B)
{
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
}



/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* naive algorithm                                            */
/*------------------------------------------------------------*/
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{

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
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points                                           */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
  
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    // double t;
    zz_pX_Multipoint_Geometric ev_geom = get_geometric_points(sz);

    // t = GetTime();
    Vec<Mat<zz_p>> valA, valB, valC;
    ev_geom.evaluate_matrix(valA, a);
    ev_geom.evaluate_matrix(valB, b);
    // cout << "geom eval: " << GetTime()-t << endl;

    // t = GetTime();
    long len = ev_geom.length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
    	mul(valC[i], valA[i], valB[i]);
    }
    // cout << "muls eval: " << GetTime()-t << endl;

    // t = GetTime();
    ev_geom.interpolate_matrix(c, valC);
    // cout << "geom interp: " << GetTime()-t << endl;
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* FFT points                                                 */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
  
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint_FFT ev_FFT = get_FFT_points(sz);

    Vec<Mat<zz_p>> valA, valB, valC;
    ev_FFT.evaluate_matrix(valA, a);
    ev_FFT.evaluate_matrix(valB, b);

    long len = ev_FFT.length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
    	mul(valC[i], valA[i], valB[i]);
    }

    ev_FFT.interpolate_matrix(c, valC);
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* chooses the kind of points                                 */
/*------------------------------------------------------------*/
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    if (is_FFT_ready())
    {
	multiply_evaluate_FFT(c, a, b);
    }
    else
    {
	multiply_evaluate_geometric(c, a, b);
    }
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);

    zz_pX_Transform_naive trs_naive(max(dA, dB) + 1);

    Vec<Mat<zz_p>> valA, valB, valC;
    trs_naive.forward_left_matrix(valA, a);
    trs_naive.forward_right_matrix(valB, b);

    long len = trs_naive.transform_length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
    	mul(valC[i], valA[i], valB[i]);
    }

    trs_naive.backward_matrix(c, valC);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    zz_pX_Transform_karatsuba trs_karatsuba = zz_pX_Transform_karatsuba();

    Vec<Mat<zz_p>> valA, valB, valC;
    trs_karatsuba.forward_left_matrix(valA, a);
    trs_karatsuba.forward_right_matrix(valB, b);

    long len = trs_karatsuba.transform_length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
    	mul(valC[i], valA[i], valB[i]);
    }

    trs_karatsuba.backward_matrix(c, valC);
}

