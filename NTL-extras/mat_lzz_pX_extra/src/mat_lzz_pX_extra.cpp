#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <algorithm> // for manipulating std::vector (min, max, ..)

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


/*------------------------------------------------------------*/
/* some comment                                               */
/*------------------------------------------------------------*/
void check_shift(bool &shifted, const std::vector<long> &shift, const Mat<zz_pX> &pmat, const bool row_wise = true)
{
	shifted = false;
	if (!shift.empty()) 
	{
		shifted = true;
		if (row_wise && (long)shift.size() != pmat.NumCols())
		{
			throw "==check_shift== Provided shift does not have the right dimension (working row-wise)";
		}
		if (!row_wise && (long)shift.size() != pmat.NumRows())
		{
			throw "==check_shift== Provided shift does not have the right dimension (working column-wise)";
		}
	}
}

/*------------------------------------------------------------*/
/* some comment                                               */
/* a supposed to be empty, just initialized */
/*------------------------------------------------------------*/
void degree_matrix(Mat<long> &degmat, const Mat<zz_pX> &pmat, 
                   const std::vector<long> & shift,
                   const bool row_wise)
{
	// check if shifted + shift dimension
	bool shifted;
	check_shift(shifted, shift, pmat, row_wise);

	// compute minimum shift entry (will be used for degree of zero entries of b)
	long min_shift = *std::min_element(shift.begin(),shift.end());

	// set the dimensions of degmat and populate it with degrees
	degmat.SetDims(pmat.NumRows(), pmat.NumCols());
	for (long i = 0; i < pmat.NumRows(); i++)
	{
		for (long j = 0; j < pmat.NumCols(); j++)
		{
			degmat[i][j] = deg(pmat[i][j]);
			if (shifted)
			{
				if (pmat[i][j] == -1)
				{
					degmat[i][j] += min_shift;
				}
				else if (row_wise)
				{
					degmat[i][j] += shift[j];
				}
				else
				{
					degmat[i][j] += shift[i];
				}
			}
		}
	}
}

/*------------------------------------------------------------*/
/* some comment                                               */
/* requirement: length of rdeg is correct */
/*------------------------------------------------------------*/
void row_degree(std::vector<long> &rdeg, const Mat<zz_pX> &pmat,
                const std::vector<long> &shift)
{ 
	// check if shifted + shift dimension
	bool shifted;
	check_shift(shifted,shift,pmat);

	// check rdeg has the right length
	if ((long)rdeg.size() != pmat.NumRows()) {
		throw "==row_degree== Provided vector does not have size = NumRows";
	}

	// retrieve the shifted degree matrix
	Mat<long> degmat;
	degree_matrix(degmat,pmat,shift,true);

	// take the max of each row of degmat
	for (long r = 0; r < pmat.NumRows(); r++)
	{
		auto max_deg = degmat[r][0];
		for (long c = 1; c < pmat.NumCols(); c++)
		{
			if (max_deg < degmat[r][c]) 
			{
				max_deg = degmat[r][c];
			}
		}
		rdeg[r] = max_deg;
	}
}

/*------------------------------------------------------------*/
/* some comment                                               */
/* requirement: length of cdeg is correct */
/*------------------------------------------------------------*/
void col_degree(std::vector<long> &cdeg, const Mat<zz_pX> &pmat,
                const std::vector<long> &shift)
{
	// check if shifted + shift dimension
	bool shifted;
	check_shift(shifted,shift,pmat,false);

	// check cdeg has the right length
	if ((long)cdeg.size() != pmat.NumCols()) {
		throw "==col_degree== Provided vector does not have size = NumCols";
	}

	// retrieve the shifted degree matrix
	Mat<long> degmat;
	degree_matrix(degmat,pmat,shift,false);

	// take the max of each column of degmat
	for (long c = 0; c < pmat.NumCols(); c++)
	{
		auto max_deg = degmat[0][c];
		for (long r = 1; r < pmat.NumRows(); r++)
		{
			if (max_deg < degmat[r][c]) 
			{
				max_deg = degmat[r][c];
			}
		}
		cdeg[c] = max_deg;
	}
} 

/*------------------------------------------------------------*/
/* some comment                                               */
/*------------------------------------------------------------*/
void leading_matrix(Mat<zz_p> &lmat,
		const Mat<zz_pX> &pmat,
		const std::vector<long> & shift,
		const bool row_wise)
{
	// check if shifted + shift dimension
	bool shifted;
	check_shift(shifted, shift, pmat, row_wise);

	// retrieve the row degree (or column degree)
	std::vector<long> degree;
	if (row_wise)
	{
		degree.resize(pmat.NumRows());
		row_degree(degree,pmat,shift);
	}
	else
	{
		degree.resize(pmat.NumCols());
		col_degree(degree,pmat,shift);
	}

	// initialize space for lmat
	lmat.SetDims(pmat.NumRows(), pmat.NumCols());
	// retrieve the leading coefficients
	for (long r = 0; r < pmat.NumRows(); r++)
	{
		for (long c = 0; c < pmat.NumCols(); c++)
		{
			long d; // actual shifted degree of the r,c entry
			// FIXME issue with zero entries? why not using degree matrix?
			d = deg(pmat[r][c]);
			if (shifted)
			{
				d += (row_wise ? shift[c] : shift[r]);
			}
			if (d == (row_wise ? degree[r] : degree[c]))
			{
				lmat[r][c] = LeadCoeff(pmat[r][c]);
			}
		}
	}
}

/*------------------------------------------------------------*/
/* some comment                                               */
/*------------------------------------------------------------*/
bool is_reduced (const Mat<zz_pX> & pmat,const std::vector<long> & shift, const bool row_wise)
{
	Mat<zz_p> lead_mat;
	leading_matrix(lead_mat,pmat,shift,row_wise);
	auto rank = gauss(lead_mat);
	return rank == (row_wise ? pmat.NumRows() : pmat.NumCols());
}

/*------------------------------------------------------------*/
/* some comment                                               */
/*------------------------------------------------------------*/
// FIXME return vector?
std::vector<long> pivot_index (
		std::vector<long> & index,
		const Mat<zz_pX> & b,
		const std::vector<long> & shift,
		const bool row_wise)
{
	bool shifted;
	check_shift(shifted, shift, b, row_wise);

	Mat<long> deg_mat;
	degree_matrix(deg_mat,b,shift,row_wise);

	std::vector<long> degree;
	if (row_wise)
	{
		row_degree(degree,b,shift);
	}
	else
	{
		col_degree(degree,b,shift);
	}

	long zero_degree = -1;
	if (shifted)
	{
		zero_degree = *std::min_element(shift.begin(),shift.end()) -1;
	}

	if (row_wise)
	{
		if ((long)index.size() != b.NumRows()) {
			throw "==pivot_index== Provided vector does not have size = NumRows";
		}
		for (long r = 0; r < b.NumRows(); r++)
		{
			if (degree[r] == zero_degree) 
			{
				index[r] = -1;
			}
			else
			{
				for (long c = 0; c <b.NumCols(); c++)
				{
					if (deg_mat[r][c] == degree[r]) 
					{
						index[r] = c;
					}
				}
			}
		}
	}
	else
	{
		if ((long)index.size() != b.NumCols()) {
			throw "==pivot_index== Provided vector does not have size = NumCols";
		}
		for(long c = 0; c < b.NumCols(); c++)
		{
			if (degree[c] == zero_degree) 
			{
				index[c] = -1;
			}
			else
			{
				for (long r = 0; r < b.NumRows(); r++)
				{
					if (deg_mat[r][c] == degree[c]) 
					{
						index[c] = r;
					}
				}
			}
		}
	}
	return degree;
}

/*
  bool is_weak_popov (const Mat<zz_pX> &m, const Vec<long> &shift = Vec<long>(), const bool row_wise = true){
  Vec<long> pivots;
  pivot_index(pivots, b, shift, row_wise);
	
  }

  bool is_popov (const Mat<zz_pX> &m, const Vec<long> &shift = Vec<long>(), const bool row_wise = true){
  Vec<long> pivots;
  Vec<long degree = pivot_index(pivots, b, shift, row_wise);
  }
*/
