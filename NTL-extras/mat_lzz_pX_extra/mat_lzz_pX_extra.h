#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector> // std vector, for shifts, degrees, pivot indices

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* UTILS                                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO Identity matrix, zero matrix


/*------------------------------------------------------------*/
/* maximum degree of the entries of pmat                      */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* print vector -- move elsewhere ???                         */
/*------------------------------------------------------------*/
std::ostream &operator<<(std::ostream &out, const std::vector<long> &s);

//TODO need equivalents of those to make storage transparent???
//const zz_p coeff(const zz_pX& a, long i);
// returns the coefficient of X^i, or zero if i not in range
//void SetCoeff(zz_pX& x, long i, zz_p a);
//void SetCoeff(zz_pX& x, long i, long a);
// makes coefficient of X^i equal to a; error is raised if i < 0
//void SetCoeff(zz_pX& x, long i);
// makes coefficient of X^i equal to 1;  error is raised if i < 0

//truncate mod X^... , for all the matrix or some columns/rows of it
//TODO different truncation orders on the different columns/rows
// full matrix versions
void trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
Mat<zz_pX> trunc(const Mat<zz_pX>& a, long n);

// row versions
void truncRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, long r, long n);
Mat<zz_pX> truncRow(const Mat<zz_pX>& a, long r, long n);

// col versions
void truncCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, long c, long n);
Mat<zz_pX> truncCol(const Mat<zz_pX>& a, long c, long n);

//TODO: multiply row or column of matrix (vec_lzz_pX) by constant

/**************************************************************************\

                               Shift Operations

LeftShift by n means multiplication by X^n
RightShift by n means division by X^n

A negative shift amount reverses the direction of the shift.

\**************************************************************************/

// operator notation:
// TODO versions with different shifting orders on different rows/columns

// full matrix, 1 row, 1 col

// full matrix shift
Mat<zz_pX> operator<< (const Mat<zz_pX> &a, long n);
Mat<zz_pX> operator>> (const Mat<zz_pX> &a, long n);

Mat<zz_pX>& operator<<=(Mat<zz_pX>& x, long n);
Mat<zz_pX>& operator>>=(Mat<zz_pX>& x, long n);

//// procedural versions:

// full matrix left shifts
void LeftShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
Mat<zz_pX> LeftShift(const Mat<zz_pX>& a, long n);

// single row shifts
void LeftShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n);
Mat<zz_pX> LeftShiftRow(const Mat<zz_pX>& a, const long r, long n);

// single col shifts
void LeftShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n);
Mat<zz_pX> LeftShiftCol(const Mat<zz_pX>& a, const long c, long n);

// full matrix right shifts
void RightShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
Mat<zz_pX> RightShift(const Mat<zz_pX>& a, long n);

// single row shifts
void RightShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n);
Mat<zz_pX> RightShiftRow(const Mat<zz_pX>& a, const long r, long n);

// single col shifts
void RightShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n);
Mat<zz_pX> RightShiftCol(const Mat<zz_pX>& a, const long c, long n);

// TODO reverse operations
//void reverse(zz_pX& x, const zz_pX& a, long hi);
//zz_pX reverse(const zz_pX& a, long hi);
//
//void reverse(zz_pX& x, const zz_pX& a);
//zz_pX reverse(const zz_pX& a);
// x = reverse of a[0]..a[hi] (hi >= -1);
// hi defaults to deg(a) in second version




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d                         */
/*------------------------------------------------------------*/
void random_mat_zz_pX(Mat<zz_pX>& a, long m, long n, long d);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* BASIC ARITHMETIC                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);
inline void add(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{
    add(c, b, a);
}

inline Mat<zz_pX> operator+(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    add(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator+(const Mat<zz_pX>& a, const Mat<zz_p>& b)
{ 
    Mat<zz_pX> x; 
    add(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator+(const Mat<zz_p>& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    add(x, a, b); 
    return x; 
}

inline Mat<zz_pX> & operator+=(Mat<zz_pX> & x, const Mat<zz_pX>& b)
{
    add(x, x, b); 
    return x; 
}

inline Mat<zz_pX> & operator+=(Mat<zz_pX> & x, const Mat<zz_p>& b)
{
    add(x, x, b); 
    return x; 
}


/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);
void sub(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

inline Mat<zz_pX> operator-(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    sub(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator-(const Mat<zz_pX>& a, const Mat<zz_p>& b)
{ 
    Mat<zz_pX> x; 
    sub(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator-(const Mat<zz_p>& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    sub(x, a, b); 
    return x; 
}

inline Mat<zz_pX> & operator-=(Mat<zz_pX> & x, const Mat<zz_pX>& b)
{
    sub(x, x, b); 
    return x; 
}

inline Mat<zz_pX> & operator-=(Mat<zz_pX> & x, const Mat<zz_p>& b)
{
    sub(x, x, b); 
    return x; 
}

/*------------------------------------------------------------*/
/* constant matrix multiplication                             */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);
void mul(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const Mat<zz_p>& b)
{ 
    Mat<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator*(const Mat<zz_p>& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}

inline Mat<zz_pX> & operator*=(Mat<zz_pX> & x, const Mat<zz_p>& b)
{
    mul(x, x, b); 
    return x; 
}

/*------------------------------------------------------------*/
/* scalar multiplication                                      */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b);

inline void mul(Mat<zz_pX> & c, const zz_p & a, const Mat<zz_pX> & b)
{
    mul(c, b, a);
}

inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const zz_p& b)
{ 
    Mat<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator*(const zz_p& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* negate                                                     */
/*------------------------------------------------------------*/
void neg(Mat<zz_pX> & x, const Mat<zz_pX> & a);

inline Mat<zz_pX> operator-(const Mat<zz_pX> & a)
{
    Mat<zz_pX> x; 
    neg(x, a); 
    return x;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SETTING AND GETTING COEFFICIENTS                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* sets x = ith coefficient of a                              */
/*------------------------------------------------------------*/
void GetCoeff(Mat<zz_p>& x, const Mat<zz_pX>& a, long i);

/*------------------------------------------------------------*/
/* returns ith coefficient matrix of a                        */
/*------------------------------------------------------------*/
inline Mat<zz_p> coeff(const Mat<zz_pX>& a, long i)
{
    Mat<zz_p> x;
    GetCoeff(x, a, i);
    return x;
}

/*------------------------------------------------------------*/
/* returns the matrix of leading coefficients                 */
/*------------------------------------------------------------*/
Mat<zz_p> matrix_of_leading_coefficients(const Mat<zz_pX>& a);

/*------------------------------------------------------------*/
/* returns constant coefficient matrix of a                   */
/*------------------------------------------------------------*/
inline Mat<zz_p> constant_coefficient(const Mat<zz_pX>& a)
{
    return coeff(a, 0);
}

/*------------------------------------------------------------*/
/* sets ith coefficient of x to a                             */
/*------------------------------------------------------------*/
void SetCoeff(Mat<zz_pX>& x, long i, Mat<zz_p> &a);


/*------------------------------------------------------------*/
/* convert from Mat<zz_p>                                     */
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_montgomery3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba4(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long len);

inline void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b){
    multiply_transform(c, a, b, max(deg(a), deg(b)) + 1);
}

void multiply_3_primes(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);


void multiply(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

// TODO: multiply given upper degree bound --> use transforms


/* inline Mat<zz_pX> & operator*=(Mat<zz_pX> & x, const Mat<zz_pX>& b) */
/* { */
/*     mul(x, x, b);  */
/*     return x;  */
/* } */


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* DEGREES, PIVOTS, LEADING MATRIX                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Shifted reduced forms of polynomials matrices. Recall that */
/* Popov => ordered weak Popov => weak Popov => Reduced       */
/*------------------------------------------------------------*/

enum PolMatForm {
	NONE = 0,
	REDUCED = 1, 
	WEAK_POPOV = 2,
	ORD_WEAK_POPOV = 3,
	POPOV = 4,
};

typedef std::vector<long> Shift;
typedef std::vector<long> DegVec;

/*------------------------------------------------------------*/
/* matrix of integers; deg(0) = -1                            */
/*------------------------------------------------------------*/
void degree_matrix(
		Mat<long> &degmat,
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift(),
		const bool row_wise=true
		);

/*------------------------------------------------------------*/
/* max degree of row entries                                  */
/*------------------------------------------------------------*/
void row_degree(
		DegVec & rdeg,
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift()
		); 

/*------------------------------------------------------------*/
/* max degree of col entries                                  */
/*------------------------------------------------------------*/
void col_degree(
		DegVec & cdeg,
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift()
		); 

/*------------------------------------------------------------*/
/* finds the pivot indices; returns the row/col degs          */
/*------------------------------------------------------------*/
void pivot_index(
		std::vector<long> & pivind,
		DegVec & pivdeg,
		const Mat<zz_pX> & pmat,
		const Shift & shift = Shift(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/* leading matrix of b                                        */
/*------------------------------------------------------------*/
void leading_matrix(
		Mat<zz_p> &lmat,
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING REDUCED/NORMAL FORMS                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* returns true if b is reduced                               */
/*------------------------------------------------------------*/
bool is_reduced(
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/* returns true if b is in weak popov form (forbide 0-row/col */
/*------------------------------------------------------------*/
bool is_weak_popov(
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift(),
		const bool row_wise = true,
		const bool ordered= false
		);

bool is_popov(
		const Mat<zz_pX> &pmat,
		const Shift & shift = Shift(),
		const bool row_wise = true,
		const bool up_to_permutation = false
		);

/*------------------------------------------------------------*/
/* if b is in some shifted reduced form,                      */
/* return the strongest detected form                         */
/* otherwise return NONE (= 0)                                */
/*------------------------------------------------------------*/
PolMatForm get_polmatform(
		const Mat<zz_pX> & pmat,
		const Shift & shift = Shift(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/* check that pmat is in the shifted reduced form             */
/* indicated by argument 'form'                               */
/*------------------------------------------------------------*/
bool is_polmatform(
		const Mat<zz_pX> & pmat,
		const PolMatForm form,
		const Shift & shift = Shift(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MINIMAL APPROXIMANT BASES                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

typedef std::vector<long> Order;

// TODO write proper docstrings

////Definition (approximant basis)
// Given:
//   * m x n matrix of univariate polynomials 'pmat',
//   * approximation order 'order' (list of n positive integers),
// An approximant basis for (pmat,order) is a matrix over K[X]
// whose rows form a basis for the K[X]-module
// { 'app' in K[X]^{1 x m}  |  the column j of 'app' 'pmat' is 0 modulo X^{order[j]} }

//// Minimal and Popov approximant bases
// Given in addition:
//   * a degree shift 'shifts' (list of m integers)
// then an approximant basis for (pmat,order) is said to be
// "a shift-minimal" (resp. "the shift-Popov") approximant basis
// if it shift-reduced (resp. in shift-Popov form)
// Idem for shift-ordered weak Popov
// Cf. literature for definitions

/*------------------------------------------------------------*/
/* TODO: documentation to explain what this computes and what */
/* are the options                                            */
/*------------------------------------------------------------*/

// TODO generic: from shift, deduce the pivot degree expected generically and
// use this as a shift instead of 'shift', then obtain directly Popov approx
// basis; skipping the first call to find the pivot degree (in addition, will
// this be more efficient because shift is nicer?)

// TODO generic: try base case at nd = m instead of d = 1 (linearization
// approach), see what improvement this brings (will improve, but for what kind
// of n?)

// TODO mbasis: threshold res_update
// TODO pmbasis: threshold mbasis

// TODO mbasis (non-res-update): at the beginning, in a single round, gather
// matrix coefficients a[i], i<order, such that pmat = sum_i a[i] X^i,
// and work with them

// TODO mbasis parallelize computation of residual (for big matrices)

// TODO mbasis-resupdate: what representation of pmat to use?

// Guarantee: output is at least ordered weak Popov
// return value is pivot degree

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

// TODO: when algorithms are ready, find threshold and write the definition
DegVec approximant_basis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Order & order,
		const Shift & shift = Shift(),
		const PolMatForm form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool generic = false
		);

// TODO: when the above function is ready, uncomment the definition below
DegVec approximant_basis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift = Shift(),
		const PolMatForm form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool generic = false
		);
//{
//	Order orders(mat.NumRows(),order);
//	return approximant_basis(appbas,mat,orders,shift,canonical,row_wise,generic);
//}


/*------------------------------------------------------------*/
/* Verifying that appbas is a shift-minimal approximant       */
/* basis for input matrix 'pmat' and order 'order'            */
/* 'form' gives the minimal requirement to check (matrix must */
/* be at least in the form 'form')                            */
/* 'randomized' says whether using a Monte Carlo or Las Vegas */
/* verification algorithm is acceptable                       */
/* Note: currently, deterministic verification is, for most   */
/* instances, as long as re-computing the basis               */
/*------------------------------------------------------------*/

bool is_approximant_basis(
		const Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Order & order,
		const Shift & shift = Shift(),
		const PolMatForm & form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool randomized = false
		);

bool is_approximant_basis(
		const Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift = Shift(),
		const PolMatForm & form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool randomized = false
		);

/*------------------------------------------------------------*/
/* Iterative algorithm for general order and shift            */
/* References:                                                */
/*   - Beckermann 1992                                        */
/*   - Van Barel-Bultheel 1991+1992                           */
/*   - Beckermann-Labahn 2000 (ensuring s-Popov)              */
/*------------------------------------------------------------*/
DegVec appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Order & order,
		const Shift & shift,
		bool order_wise=true
		);

DegVec popov_appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Order & order,
		const Shift & shift,
		bool order_wise=true
		);

/*------------------------------------------------------------*/
/* M-Basis algorithm for approximant order = 1                */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018 (ensuring s-Popov)       */
/*------------------------------------------------------------*/
// input: kerbas is constant, will contain the left kernel of pmat in reduced REF
// output: pivot degrees of the approx basis (also indicates where the rows of
// kernel should appear in approx basis)
DegVec popov_mbasis1(
		Mat<zz_p> & kerbas,
		const Mat<zz_p> & pmat,
		const Shift & shift
		);

/*------------------------------------------------------------*/
/* M-Basis algorithm for uniform approximant order            */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/
DegVec mbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift
		);

// TODO some thresholding to be done, so that mbasis does the
// resupdate strategy when it is faster
// --> organize code so that mbasis is always ~the fastest
DegVec mbasis_resupdate(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift
		);

DegVec popov_mbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift
		);

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform approximant order           */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shifts) */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/
DegVec pmbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift
		);

DegVec popov_pmbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const Shift & shift
		);


/**********************************************************************
*                     MINIMAL INTERPOLANT BASES                      *
**********************************************************************/

// list of points: for each column, we have a list of points (elt from zz_p and
// multiplicity)
// FIXME not thought thorougly yet, subject to change
typedef std::vector<std::vector<std::pair<zz_p,long>>> Points;

////Definition (interpolant basis)
// Given:
//   * m x n matrix of univariate polynomials 'pmat',
//   * list 'points' of n pairs [root,multiplicity], which define n products of linear factors M_0,...,M_{n-1},
// An interpolant basis for (pmat,points) is a matrix over K[X]
// whose rows form a basis for the K[X]-module
// { 'int' in K[X]^{1 x m}  |  the column j of 'app' 'pmat' is 0 modulo M_j }

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

// TODO
DegVec interpolant_basis(
		Mat<zz_pX> & intbas,
		const Mat<zz_pX> & pmat,
		const Points & pts,
		const Shift & shift = Shift(),
		const PolMatForm form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool generic = false
		);

// TODO uncomment below when above is ready
// below, the following is called "uniform interpolation case" (or case with uniform points)
// this means that we have the same points on all columns, all with multiplicity one)
// (FIXME could be easily generalized to any constant multiplicity for all...?)
DegVec interpolant_basis(
		Mat<zz_pX> & intbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift = Shift(),
		const PolMatForm form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool generic = false
		);
//{
//	std::vector<std::pair<zz_p,long>> list_pts(pts.size());
//	for ( long i=0; i<list_pts.size(); ++i )
//	{
//		list_pts[i] = std::pair<zz_p,long>(pts[i],1);
//	}
//	Points points(mat.NumCols(),list_pts);
//	return interpolant_basis(appbas,mat,points,shift,canonical,row_wise,generic);
//}


/*------------------------------------------------------------*/
/* Verifying that intbas is a shift-minimal interpolant       */
/* basis for input matrix 'pmat' and points 'points'          */
/* 'form' gives the minimal requirement to check (matrix must */
/* be at least in the form 'form')                            */
/* 'randomized' says whether using a Monte Carlo or Las Vegas */
/* verification algorithm is acceptable                       */
/* Note: currently, deterministic verification is, for most   */
/* instances, as long as re-computing the basis               */
/*------------------------------------------------------------*/

// TODO not implemented yet
bool is_interpolant_basis(
		const Mat<zz_pX> & intbas,
		const Mat<zz_pX> & pmat,
		const Points & pts,
		const Shift & shift = Shift(),
		const PolMatForm & form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool randomized = false
		);

// TODO (uniform interpolation variant)
bool is_interpolant_basis(
		const Mat<zz_pX> & intbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift = Shift(),
		const PolMatForm & form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool randomized = false
		);

/*------------------------------------------------------------*/
/* Iterative algorithm for arbitrary points and shift         */
/* References:                                                */
/*   - Beckermann 1992                                        */
/*   - Van Barel-Bultheel 1991+1992                           */
/*   - Beckermann-Labahn 2000 (ensuring s-Popov)              */
/*------------------------------------------------------------*/
DegVec intbas_iterative(
		Mat<zz_pX> & intbas,
		const Mat<zz_pX> & pmat,
		const Points & pts,
		const Shift & shift,
		bool order_wise=true
		);

DegVec popov_intbas_iterative(
		Mat<zz_pX> & intbas,
		const Mat<zz_pX> & pmat,
		const Points & pts,
		const Shift & shift,
		bool order_wise=true
		);


/*------------------------------------------------------------*/
/* Adaptation of M-Basis for uniform interpolation points     */
/*------------------------------------------------------------*/

// --> popov_mbasis1 can be called as such (with, as input, pmat evaluated at a
// point)

DegVec mbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift
		);

// TODO some thresholding to be done, so that mbasis does the
// resupdate strategy when it is faster
// --> organize code so that mbasis is always ~the fastest
DegVec mbasis_resupdate(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift
		);

DegVec popov_mbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift
		);

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform interpolation points        */
/*------------------------------------------------------------*/
DegVec pmbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift
		);

DegVec popov_pmbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & pmat,
		const Vec<zz_p> & pts,
		const Shift & shift
		);

/**********************************************************************
*                            KERNEL BASIS                            *
**********************************************************************/

/****************************
*  Kernel via approximant  *
****************************/

// TODO generic case

// TODO general case via fast s-Popov appbas (fastest known approach for bad shifts)

/*************************************************
*  Zhou-Labahn-Storjohann minimal kernel basis  *
*************************************************/

// TODO

/******************************************
*  Application to linear system solving  *
******************************************/

// TODO

/*************************************************
*  Fast left kernel for small column dimension  *
*************************************************/

// TODO-long shot via relation basis (itself todo), cf work with Xuan

/**********************************************************************
*                             ROW BASIS                              *
**********************************************************************/

/*************************************
*  Li and Storjohann's compression  *
*************************************/

// TODO, cf Chao Li's 2007 Master's thesis
// does not really give a row basis, but allows to reduce to almost
// square case by an efficient Las Vegas algorithm

/***************************
*  Zhou-Labahn algorithm  *
***************************/

// TODO

/**********************************************************************
*                         INVERSE EXPANSION                          *
**********************************************************************/

// TODO Newton iteration
// --> polynomial matrix division with remainder (cf. e.g. Neiger-Vu 2017)

// TODO see Storjohann 2003 for high-order lifting
// TODO application to linear system solving
// --> see which other consequences of high-order lifting may be worth implementing


/**********************************************************************
*                             INVERSION                              *
**********************************************************************/

// TODO evaluation-interpolation approach

// TODO Jeannerod-Villard: worth implementing?

// TODO Storjohann 2015 fast Las Vegas

// TODO Zhou-Labahn-Storjohann 2014 fast deterministic

/**********************************************************************
*                          BASIS REDUCTION                           *
*            (shifted reduced form and shifted normal forms)         *
**********************************************************************/

// TODO general reduction to uniform shift via pre-multiplication
// worthwile at least when shift close to uniform

// TODO naive algorithms (see Mulders-Storjohann for good reference)

// TODO general shifted Popov form via kernel (itself via approximant basis)

// TODO understand if there is any chance Alekhnovich improves over the
// kernel approach

// TODO nonsingular: Giorgi-Jeannerod-Villard's Las Vegas reduction
// (worth implementing for shifts other than uniform?)


/**********************************************************************
*                    TRIANGULARIZATION ALGORITHMS                    *
**********************************************************************/

/**************************************************
*  Labahn-Neiger-Zhou partial triangularization  *
**************************************************/

// TODO one step (uses kernel + row basis)
// TODO diagonal entries of Hermite form

/***********************************************
*  Labahn-Neiger-Zhou Hermite form algorithm  *
***********************************************/

// TODO (requires partial linearization + basis reduction)


/**********************************************************************
*                       DETERMINANT ALGORITHMS                       *
**********************************************************************/

/*******************************************************************
*  Labahn-Neiger-Zhou: via diagonal entries of triangularization  *
*******************************************************************/

// TODO general version
// TODO generic version: much simpler (only requires generic kernel + products)

// TODO other determinant algorithms??
// --> could rely on x-Smith decomposition of Gupta et al (worth
// implementing??), cf Appendix of LaNeZh17 
// --> no randomized faster approach?


#endif // MAT_LZZ_PX_EXTRA__H
