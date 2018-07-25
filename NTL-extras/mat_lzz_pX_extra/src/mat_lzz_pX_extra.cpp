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

<<<<<<< HEAD
long min (const Vec<long> &v){
	if (v.length() == 0) return 0;
	long m = v[0];
	for (long i = 1; i < v.length(); i++)
		if (v[i] < m) m = v[i];
	return m;
}

long max(const Vec<long> &v){
	if (v.length() == 0) return 0;
	long m = v[0];
	for (long i = 1; i < v.length(); i++)
		if (v[i] > m) m = v[i];
	return m;
}

void check_shift(bool &shifted, const Vec<long> &shift, const Mat<zz_pX> &b, const bool row_wise = true){
	shifted = false;
	if (shift.length() != 0) {
		shifted = true;
		if (row_wise && shift.length() != b.NumCols())
			throw "bad row shift dimensions";
		if (!row_wise && shift.length() != b.NumRows())
			throw "bad col shift  dimensions";
	}
}

void degree_matrix(Mat<long> &a, const Mat<zz_pX> &b, 
                   const Vec<long> & shift = Vec<long>(),
                   const bool row_wise = true){
	// check for shifts
	bool shifted;
	check_shift(shifted, shift, b, row_wise);
	long min_shift = min(shift);
	
	a.SetDims(b.NumRows(), b.NumCols());
	for (long i = 0; i < b.NumRows(); i++){
		for (long j = 0; j < b.NumCols(); j++){
			a[i][j] = deg(b[i][j]);
			if (shifted){
				if (a[i][j] == -1) a[i][j] += min_shift;
				else if (row_wise) a[i][j] += shift[j];
				else a[i][j] += shift[i];
			}
		}
	}
}

void row_degree(Vec<long> &a, const Mat<zz_pX> &b,
                const Vec<long> &shift = Vec<long>()){ 
  bool shifted;
  check_shift(shifted,shift,b);
  
	a.SetLength(b.NumRows());
	Mat<long> deg_mat;
	degree_matrix(deg_mat,b,shift,true);
	
	for (long r = 0; r < b.NumRows(); r++){
		auto max_deg = deg_mat[r][0];
		for (long c = 1; c < b.NumCols(); c++){
			if (max_deg < deg_mat[r][c]) max_deg = deg_mat[r][c];
		}
		a[r] = max_deg;
	}
}

void col_degree(Vec<long> &a, const Mat<zz_pX> &b,
                const Vec<long> &shift = Vec<long>()){
	bool shifted;
	check_shift(shifted,shift,b,false);
	
	a.SetLength(b.NumCols());
	Mat<long> deg_mat;
	degree_matrix(deg_mat,b,shift,false);
	
	a.SetLength(b.NumCols());
	for (long c = 0; c < b.NumCols(); c++){
		auto max_deg = deg_mat[0][c];
		for (long r = 1; r < b.NumRows(); r++){
			if (max_deg < deg_mat[r][c]) max_deg = deg_mat[r][c];
		}
		a[c] = max_deg;
	}
} 

void leading_matrix(Mat<zz_p> &a, const Mat<zz_pX> &b, 
                   const Vec<long> & shift = Vec<long>(),
                   const bool row_wise = true){
	// check for shifts
	bool shifted;
	check_shift(shifted, shift, b, row_wise);
	
	Vec<long> degree;
	if (row_wise)
		row_degree(degree,b,shift);
	else
		col_degree(degree,b,shift);
	
	a.SetDims(b.NumRows(), b.NumCols());
	for (long r = 0; r < b.NumRows(); r++)
		for (long c = 0; c < b.NumCols(); c++){
			long d,d2;
			if (row_wise) d = degree[r];
			else d = degree[c];
			d2 = deg(b[r][c]);
			if (shifted){
				if (row_wise) d2 += shift[c];
				else d2 += shift[r];
			}
			if (d2 >= d)
				a[r][c] = LeadCoeff(b[r][c]);
		}
}

bool is_reduced (const Mat<zz_pX> &b,const Vec<long> & shift = Vec<long>(), const bool row_wise = true){
	Mat<zz_p> lead_mat;
	leading_matrix(lead_mat,b,shift,row_wise);
	auto rank = gauss(lead_mat);
	if (row_wise)
		return b.NumRows() == rank;
	return b.NumCols() == rank;
}

Vec<long> pivot_index (Vec<long> &index, const Mat<zz_pX> &b,const Vec<long> & shift = Vec<long>(), const bool row_wise = true){
	bool shifted;
	check_shift(shifted, shift, b, row_wise);

	Mat<long> deg_mat;
	degree_matrix(deg_mat,b,shift,row_wise);
	
	Vec<long> degree;
	if (row_wise) row_degree(degree,b,shift);
	else col_degree(degree,b,shift);
	long zero_degree = -1;
	if (shifted)
		zero_degree = min(shift) -1;
	
	if (row_wise){
		index.SetLength(b.NumRows());
		for (long r = 0; r < b.NumRows(); r++){
			if (degree[r] == zero_degree) index[r] = -1;
			else
				for (long c = 0; c <b.NumCols(); c++)
					if (deg_mat[r][c] == degree[r]) index[r] = c;
		}
	}else{
		index.SetLength(b.NumCols());
		for(long c = 0; c < b.NumCols(); c++){
			if (degree[c] == zero_degree) index[c] = -1;
			else
				for (long r = 0; r < b.NumRows(); r++)
					if (deg_mat[r][c] == degree[c]) index[c] = r;
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



















































=======
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

>>>>>>> da85d4067df80b2e237028f257b724480c1ebbd6
