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

/*------------------------------------------------------------*/
/* print vector -- move elsewhere ???                         */
/*------------------------------------------------------------*/
std::ostream &operator<<(std::ostream &out, const std::vector<long> &s);

//TODO equivalents of those:
//const zz_p coeff(const zz_pX& a, long i);
// returns the coefficient of X^i, or zero if i not in range
//void SetCoeff(zz_pX& x, long i, zz_p a);
//void SetCoeff(zz_pX& x, long i, long a);
// makes coefficient of X^i equal to a; error is raised if i < 0
//void SetCoeff(zz_pX& x, long i);
// makes coefficient of X^i equal to 1;  error is raised if i < 0

//TODO : left and right shifts (multiplication/division by powers of X),
//for all the matrix or some columns/rows of it

//TODO : truncate mod X^... , for all the matrix or some columns/rows of it
//void trunc(zz_pX& x, const zz_pX& a, long n); // x = a % X^n
//zz_pX trunc(const zz_pX& a, long n);

//TODO: multiply row or column of matrix (vec_lzz_pX) by constant

/**************************************************************************\

                               Shift Operations

LeftShift by n means multiplication by X^n
RightShift by n means division by X^n

A negative shift amount reverses the direction of the shift.

\**************************************************************************/

// operator notation:

// full matrix, 1 row, 1 col

//zz_pX operator<<(const zz_pX& a, long n);
// full matrix shift
Mat<zz_pX> operator<< (const Mat<zz_pX> &a, long n);

//zz_pX operator>>(const zz_pX& a, long n);
Mat<zz_pX> operator>> (const Mat<zz_pX> &a, long n);

//
//zz_pX& operator<<=(zz_pX& x, long n);
Mat<zz_pX>& operator<<=(Mat<zz_pX>& x, long n);
//zz_pX& operator>>=(zz_pX& x, long n);
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


//
//void RightShift(zz_pX& x, const zz_pX& a, long n);
//zz_pX RightShift(const zz_pX& a, long n);
// full matrix left shifts
void RightShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
Mat<zz_pX> RightShift(const Mat<zz_pX>& a, long n);

// single row shifts
void RightShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n);
Mat<zz_pX> RightShiftRow(const Mat<zz_pX>& a, const long r, long n);

// single col shifts
void RightShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n);
Mat<zz_pX> RightShiftCol(const Mat<zz_pX>& a, const long c, long n);

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
/* maximum degree of the entries of pmat                      */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* matrix of integers; deg(0) = -1                            */
/*------------------------------------------------------------*/
void degree_matrix(
		Mat<long> &degmat,
		const Mat<zz_pX> &pmat,
		const std::vector<long>& shift=std::vector<long>(),
		const bool row_wise=true
		);

/*------------------------------------------------------------*/
/* max degree of row entries                                  */
/*------------------------------------------------------------*/
void row_degree(
		std::vector<long> &rdeg,
		const Mat<zz_pX> &pmat,
		const std::vector<long>& shift=std::vector<long>()
		); 

/*------------------------------------------------------------*/
/* max degree of col entries                                  */
/*------------------------------------------------------------*/
void col_degree(
		std::vector<long> &cdeg,
		const Mat<zz_pX> &pmat,
		const std::vector<long>& shift=std::vector<long>()
		); 

/*------------------------------------------------------------*/
/* finds the pivot indices; returns the row/col degs          */
/*------------------------------------------------------------*/
void pivot_index(
		std::vector<long> &pivind,
		std::vector<long> &pivdeg,
		const Mat<zz_pX> &pmat,
		const std::vector<long> & shift = std::vector<long>(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/* leading matrix of b                                        */
/*------------------------------------------------------------*/
void leading_matrix(
		Mat<zz_p> &lmat,
		const Mat<zz_pX> &pmat,
		const std::vector<long> & shift = std::vector<long>(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING REDUCED/NORMAL FORMS                               */
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

/*------------------------------------------------------------*/
/* returns true if b is reduced                               */
/*------------------------------------------------------------*/
bool is_reduced(
		const Mat<zz_pX> &pmat,
		const std::vector<long> & shift = std::vector<long>(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/* returns true if b is in weak popov form (forbide 0-row/col */
/*------------------------------------------------------------*/
bool is_weak_popov(
		const Mat<zz_pX> &pmat,
		const std::vector<long> &shift = std::vector<long>(),
		const bool row_wise = true,
		const bool ordered= false
		);

bool is_popov(
		const Mat<zz_pX> &pmat,
		const std::vector<long> &shift = std::vector<long>(),
		const bool row_wise = true,
		const bool up_to_permutation = false
		);

/*------------------------------------------------------------*/
/* if b is in some shifted reduced form,                      */
/* return the strongest detected form                         */
/* otherwise return NONE (= 0)                                */
/*------------------------------------------------------------*/
PolMatForm get_polmatform(
		const Mat<zz_pX> &pmat,
		const std::vector<long> &shift = std::vector<long>(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/* check that pmat is in the shifted reduced form             */
/* indicated by argument 'form'                               */
/*------------------------------------------------------------*/
bool is_polmatform(
		const Mat<zz_pX> &pmat,
		const PolMatForm form,
		const std::vector<long> &shift = std::vector<long>(),
		const bool row_wise = true
		);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MINIMAL APPROXIMANT BASES                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

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

// Guarantee: output is at least ordered weak Popov
// return value is pivot degree

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

std::vector<long> approximant_basis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const std::vector<long> & order,
		const std::vector<long> & shift = std::vector<long>(),
		const PolMatForm form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool generic = false
		);

std::vector<long> approximant_basis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift = std::vector<long>(),
		const PolMatForm form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool generic = false
		);
//{
//	std::vector<long> orders(mat.NumRows(),order);
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
		const std::vector<long> & order,
		const std::vector<long> & shift = std::vector<long>(),
		const PolMatForm & form = ORD_WEAK_POPOV,
		const bool row_wise = true,
		const bool randomized = false
		);

bool is_approximant_basis(
		const Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift = std::vector<long>(),
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
std::vector<long> appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const std::vector<long> order,
		const std::vector<long> & shift,
		bool order_wise=true
		);

std::vector<long> popov_appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const std::vector<long> order,
		const std::vector<long> & shift,
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
std::vector<long> popov_mbasis1(
		Mat<zz_p> & kerbas,
		const Mat<zz_p> & pmat,
		const std::vector<long> & shift
		);

/*------------------------------------------------------------*/
/* M-Basis algorithm for uniform approximant order            */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/
std::vector<long> mbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift
		);

std::vector<long> popov_mbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift
		);

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform approximant order           */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shifts) */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/
std::vector<long> pmbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift
		);

std::vector<long> popov_pmbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift
		);


#endif // MAT_LZZ_PX_EXTRA__H
