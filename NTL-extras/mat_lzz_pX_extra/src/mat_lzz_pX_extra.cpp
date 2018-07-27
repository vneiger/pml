#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* UTILS                                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s){
	out << "[ ";
	for (auto &i: s)
		out << i << " ";
	return out << "]";
}

/*------------------------------------------------------------*/
/* operators                                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
Mat<zz_pX>& operator<<=(Mat<zz_pX>& x, long n){
	for (long r = 0; r < x.NumRows(); r++)
		for (long c = 0; c < x.NumCols(); c++)
			x[r][c] <<= n;
	return x;
}

Mat<zz_pX>& operator>>=(Mat<zz_pX>& x, long n){
	for (long r = 0; r < x.NumRows(); r++)
		for (long c = 0; c < x.NumCols(); c++)
			x[r][c] >>= n;
	return x;
}

Mat<zz_pX> operator<< (const Mat<zz_pX> &a, long n){
	Mat<zz_pX> res{a};
	res <<= n;
	return res;
}

Mat<zz_pX> operator>> (const Mat<zz_pX> &a, long n){
	Mat<zz_pX> res{a};
	res >>= n;
	return res;
}

void LeftShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n){
	x = a << n;
}

Mat<zz_pX> LeftShift(const Mat<zz_pX>& a, long n){
	return a << n;
}

void LeftShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n){
	x = a;
	for (long c = 0; c < x.NumCols(); c++)
		x[r][c] <<= n;
}

Mat<zz_pX> LeftShiftRow(const Mat<zz_pX>& a, const long r, long n){
	auto x = a;
	for (long c = 0; c < x.NumCols(); c++)
		x[r][c] <<= n;
	return x;
}

void LeftShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n){
	x = a;
	for (long r = 0; r < x.NumRows(); r++)
		x[r][c] <<= n;
}

Mat<zz_pX> LeftShiftCow(const Mat<zz_pX>& a, const long c, long n){
	auto x = a;
	for (long r = 0; r < x.NumCols(); r++)
		x[r][c] <<= n;
	return x;
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

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
			throw std::invalid_argument("==check_shift== Provided shift does not have the right dimension (working row-wise)");
		}
		if (!row_wise && (long)shift.size() != pmat.NumRows())
		{
			throw std::invalid_argument("==check_shift== Provided shift does not have the right dimension (working column-wise)");
		}
	}
}

/*------------------------------------------------------------*/
/* some comment                                               */
/* degmat supposed to be empty, just initialized */
/*------------------------------------------------------------*/
void degree_matrix(Mat<long> &degmat, const Mat<zz_pX> &pmat, 
                   const std::vector<long> & shift,
                   const bool row_wise)
{
	// check if shifted + shift dimension
	bool shifted;
	check_shift(shifted, shift, pmat, row_wise);

	// compute minimum shift entry (will be used for degree of zero entries of b)
	long min_shift = (shifted ? *std::min_element(shift.begin(),shift.end()) : 0);
	// set the dimensions of degmat and populate it with degrees
	degmat.SetDims(pmat.NumRows(), pmat.NumCols());
	for (long i = 0; i < pmat.NumRows(); i++)
	{
		for (long j = 0; j < pmat.NumCols(); j++)
		{
			degmat[i][j] = deg(pmat[i][j]);
			if (shifted)
			{
				if (degmat[i][j] == -1)
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
		throw std::invalid_argument("==row_degree== Provided vector does not have size = NumRows");
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
		throw std::invalid_argument("==col_degree== Provided vector does not have size = NumCols");
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
			long d = deg(pmat[r][c]); // actual shifted degree of the r,c entry

			if (shifted)
				d += (row_wise ? shift[c] : shift[r]);

			if (d == (row_wise ? degree[r] : degree[c]))
				lmat[r][c] = LeadCoeff(pmat[r][c]);
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
void pivot_index (
		std::vector<long> & pivind,
		std::vector<long> & pivdeg,
		const Mat<zz_pX> & pmat,
		const std::vector<long> & shift,
		const bool row_wise)
{
	// check if shifted + shift dimension
	bool shifted;
	check_shift(shifted, shift, pmat, row_wise);

	// retrieve (shifted) degree matrix
	Mat<long> deg_mat;
	degree_matrix(deg_mat,pmat,shift,row_wise);

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

	long zero_degree = -1;
	if (shifted)
		zero_degree = *std::min_element(shift.begin(),shift.end()) -1;

	if (row_wise)
	{
		if ((long)pivind.size() != pmat.NumRows())
			throw std::invalid_argument("==pivot_index== Provided vector does not have size = NumRows");

		for (long r = 0; r < pmat.NumRows(); ++r)
		{
			if (degree[r] == zero_degree) 
			{
				pivdeg[r] = -1;
				pivind[r] = -1;
			}
			else
			{
				for (long c = 0; c <pmat.NumCols(); ++c)
				{
					if (deg_mat[r][c] == degree[r]) 
					{
						pivdeg[r] = deg(pmat[r][c]);
						pivind[r] = c;
					}
				}
			}
		}
	}
	else
	{
		if ((long)pivind.size() != pmat.NumCols())
			throw std::invalid_argument("==pivot_index== Provided vector does not have size = NumCols");

		for(long c = 0; c < pmat.NumCols(); c++)
		{
			if (degree[c] == zero_degree) 
			{
				pivdeg[c] = -1;
				pivind[c] = -1;
			}
			else
			{
				for (long r = 0; r < pmat.NumRows(); r++)
				{
					if (deg_mat[r][c] == degree[c]) 
					{
						pivdeg[c] = deg(pmat[r][c]);
						pivind[c] = r;
					}
				}
			}
		}
	}
}

bool is_weak_popov (
		const Mat<zz_pX> &pmat,
		const std::vector<long> &shift,
		const bool row_wise,
		const bool ordered)
{
	//retrieve pivot index
	std::vector<long> pivots;
	std::vector<long> pivdeg;
	pivots.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
	pivdeg.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
	pivot_index(pivots, pivdeg, pmat, shift, row_wise);

	// forbide zero vectors
	if( std::find(pivots.begin(),pivots.end(),-1) != pivots.end() )
		return false;

	if (!ordered)
	{ // check for pairwise distinct: sort and check no adjacent equal elements
		std::sort(pivots.begin(),pivots.end());
		if (std::adjacent_find(pivots.begin(),pivots.end()) != pivots.end())
			return false;
	}
	else
	{ // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
		if (std::adjacent_find(pivots.begin(),pivots.end(),std::greater_equal<long>()) != pivots.end())
			return false;
	}
	return true;
}

bool is_monic(const zz_pX &p){
	return IsOne(LeadCoeff(p));
}

bool is_popov (const Mat<zz_pX> &pmat,
		const std::vector<long> &shift,
		const bool row_wise,
		const bool up_to_permutation){
	if (!is_weak_popov(pmat,shift,row_wise,!up_to_permutation))
		return false;
	
	
	std::vector<long> pivots;
	std::vector<long> degrees;
	pivots.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
	degrees.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
	pivot_index(pivots,degrees,pmat,shift,row_wise);
	for (unsigned long i = 0; i < pivots.size(); i++){
		auto index = pivots[i];
		if (index >= 0){
			if(row_wise){
				if (!is_monic(pmat[i][index]))
					return false;
				for (long k=0; k < pmat.NumRows(); k++){
					if (deg(pmat[k][index]) >= degrees[i] && ((unsigned long)k != i))
						return false;
				}
			}else{ // col-wise
				if (!is_monic(pmat[index][i]))
					return false;
				for (long k = 0; k < pmat.NumCols(); k++)
					if (deg(pmat[index][k]) >= degrees[i] && ((unsigned long)k != i))
						return false;
			}
		}
	}
	return true;
}

Mat<zz_pX> identity_matrix(const long n){
	Mat<zz_pX> res;
	res.SetDims(n,n);
	for (long i = 0; i < n; i++)
		res[i][i] = zz_pX(1);
	return res;
}

/*
void weak_popov_form(Mat<zz_pX> &wpf, const Mat<zz_pX> &pmat, const std::vector<long> &shift){
	wpf = pmat;
	m = wpf.NumRows();
	n = wpf.NumCols();
	trans = identity_matrix(m);
	
	// populate shift with zeros if empty
	if (shift.length() == 0){ 
		for (long i = 0; i < n; i++)
			shift.emplace_back(0);
	}
	
	// pivots[i] = shift-pivot index of the row i of wpf
	// degrees = shift-row degree of wpf
	std::vector<long> pivots;
	std::vector<long> degrees;
	pivot_index(pivots,degress,wpf,shift,true);
	
	// rows_with_pivot[p] = indices of the rows that have pivot index p
	std::vector<std::vector<long>> rows_with_pivot(n, std::vector<long>());
	for (long i = 0; i < n; i++)
		rows_with_pivot[pivots[i]].emplace_back(i);
}
*/




PolMatForm get_polmatform(
		const Mat<zz_pX> &pmat,
		const std::vector<long> &shift,
		const bool row_wise
		)
{
	if (is_popov(pmat,shift,row_wise)) {   // TODO waiting for is_popov
		return POPOV;
	}
	if (is_weak_popov(pmat,shift,row_wise,true))
		return ORD_WEAK_POPOV;
	else if (is_weak_popov(pmat,shift,row_wise))
		return WEAK_POPOV;
	else if (is_reduced(pmat,shift,row_wise))
		return REDUCED;
	else
		return NONE;
}

bool is_polmatform(
		const Mat<zz_pX> &pmat,
		const PolMatForm form,
		const std::vector<long> &shift,
		const bool row_wise
		)
{
	switch (form)
	{
		case NONE: return true;
		case REDUCED: return is_reduced(pmat,shift,row_wise);
		case WEAK_POPOV: return is_weak_popov(pmat,shift,row_wise,false);
		case ORD_WEAK_POPOV: return is_weak_popov(pmat,shift,row_wise,true);
		case POPOV: return is_popov(pmat,shift,row_wise);
		default: throw std::invalid_argument("==is_polmatform== Unknown required polynomial matrix form.");
	}
}


