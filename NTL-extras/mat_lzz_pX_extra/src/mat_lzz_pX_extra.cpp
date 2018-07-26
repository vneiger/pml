#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* zz_pX addition                                             */
/*------------------------------------------------------------*/
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
    {
	LogicError("dimension mismatch in matrix addition");
    }

    c.SetDims(m, n);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < n; j++)
	{
	    c[i][j] = a[i][j] + b[i][j];
	}
    }
}

/*------------------------------------------------------------*/
/* addition, rhs is constant                                  */
/*------------------------------------------------------------*/
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
    {
	LogicError("dimension mismatch in matrix addition");
    }

    c.SetDims(m, n);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < n; j++)
	{
	    c[i][j] = a[i][j] + b[i][j];
	}
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* zz_pX subtraction                                          */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
    {
	LogicError("dimension mismatch in matrix subtraction");
    }

    c.SetDims(m, n);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < n; j++)
	{
	    c[i][j] = a[i][j] - b[i][j];
	}
    }
}

/*------------------------------------------------------------*/
/* subtraction, rhs is constant                               */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
    {
	LogicError("dimension mismatch in matrix subtraction");
    }

    c.SetDims(m, n);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < n; j++)
	{
	    c[i][j] = a[i][j] - b[i][j];
	}
    }
}

/*------------------------------------------------------------*/
/* subtraction, lhs is constant                               */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
    {
	LogicError("dimension mismatch in matrix subtraction");
    }

    c.SetDims(m, n);
    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j < n; j++)
	{
	    c[i][j] = a[i][j] - b[i][j];
	}
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* constant matrix multiplication                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* multiplication, rhs is a constant matrix                   */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b)
{
    long d = deg(a);

    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    if (n != b.NumRows())
    {
	LogicError("dimension mismatch in constant matrix multplication");
    }

    c.SetDims(m, p);
    Mat<zz_p> tmp, res;
    tmp.SetDims(m, p);
    for (long i = 0; i <= d; i++)
    {
	for (long u = 0; u < m; u++)
	{
	    for (long v = 0; v < n; v++)
	    {
		tmp[u][v] = coeff(a[u][v], i);
	    }
	}
	res = tmp * b;
	for (long u = 0; u < m; u++)
	{
	    for (long v = 0; v < p; v++)
	    {
		SetCoeff(c[u][v], i, res[u][v]);
	    }
	}
    }
}

/*------------------------------------------------------------*/
/* multiplication, lhs is a constant matrix                   */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{
    long d = deg(b);

    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    if (n != b.NumRows())
    {
	LogicError("dimension mismatch in constant matrix multplication");
    }

    c.SetDims(m, p);
    Mat<zz_p> tmp, res;
    tmp.SetDims(m, p);
    for (long i = 0; i <= d; i++)
    {
	for (long u = 0; u < n; u++)
	{
	    for (long v = 0; v < p; v++)
	    {
		tmp[u][v] = coeff(b[u][v], i);
	    }
	}
	res = a * tmp;
	for (long u = 0; u < m; u++)
	{
	    for (long v = 0; v < p; v++)
	    {
		SetCoeff(c[u][v], i, res[u][v]);
	    }
	}
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* scalar multiplication                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b){

    long m = a.NumRows();
    long n = a.NumCols();

    c.SetDims(m, n);

    for (long u = 0; u < m; u++)
    {
	for (long v = 0; v < n; v++)
	{
	    c[u][v] = b * a[u][v];
	}
    }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* negate                                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void neg(Mat<zz_pX> & x, const Mat<zz_pX> & a)
{
   long m = a.NumRows();
   long n = a.NumCols();
   
   x.SetDims(m, n);
   
    for (long u = 0; u < m; u++)
    {
	for (long v = 0; v < n; v++)
	{
	    NTL::negate(x[u][v], a[u][v]);
	}
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* get / set coefficients                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets x = ith coefficient of a                              */
/*------------------------------------------------------------*/
void GetCoeff(Mat<zz_p>& x, const Mat<zz_pX>& a, long i)
{
    long m = a.NumRows();
    long n = a.NumCols();
   
    x.SetDims(m, n);
    
    for (long u = 0; u < m; u++)
    {
	for (long v = 0; v < n; v++)
	{
	    x[u][v] = coeff(a[u][v], i);
	}
    }
}


/*------------------------------------------------------------*/
/* returns the matrix of leading coefficients                 */
/*------------------------------------------------------------*/
Mat<zz_p> matrix_of_leading_coefficients(const Mat<zz_pX>& a)
{
    long m = a.NumRows();
    long n = a.NumCols();
    
    Mat<zz_p> x;
    x.SetDims(m, n);
    
    for (long u = 0; u < m; u++)
    {
	for (long v = 0; v < n; v++)
	{
	    x[u][v] = LeadCoeff(a[u][v]);
	}
    }

    return x;
}


/*------------------------------------------------------------*/
/* sets ith coefficient of x to a                             */
/*------------------------------------------------------------*/
void SetCoeff(Mat<zz_pX>& x, long i, Mat<zz_p> &a)
{
    long m = x.NumRows();
    long n = x.NumCols();
    
    if (m != a.NumRows() || n != a.NumCols())
    {
	LogicError("dimension mismatch in matrix SetCoeff");
    }

    if (i < 0)
    {
	LogicError("negative index in matrix SetCoeff");
    }

    for (long u = 0; u < m; u++)
    {
	for (long v = 0; v < n; v++)
	{
	    SetCoeff(x[u][v], i, a[u][v]);
	}
    }
}


/*------------------------------------------------------------*/
/* random (n, m) matrix of degree < d                         */
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

    zz_pX_Multipoint_Geometric ev_geom = get_geometric_points(sz);

    Vec<Mat<zz_p>> valA, valB, valC;
    ev_geom.evaluate_matrix(valA, a);
    ev_geom.evaluate_matrix(valB, b);

    long len = ev_geom.length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
    	mul(valC[i], valA[i], valB[i]);
    }

    ev_geom.interpolate_matrix(c, valC);
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

