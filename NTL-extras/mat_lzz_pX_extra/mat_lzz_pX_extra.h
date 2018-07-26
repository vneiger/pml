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
/* print vector                                               */
/*------------------------------------------------------------*/
std::ostream &operator<<(std::ostream &out, const std::vector<long> &s);


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


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MINIMAL APPROXIMANT BASES                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* TODO: documentation to explain what this computes and what */
/* are the options                                            */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/
// choice: output s-pivot degree
std::vector<long> approximant_basis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & mat,
		const std::vector<unsigned long> & order,
		const std::vector<long> & shift = std::vector<long>(),
		const bool canonical = true,
		const bool row_wise = true,
		const bool generic = false
		);

std::vector<long> approximant_basis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & mat,
		const unsigned long order,
		const std::vector<long> & shift = std::vector<long>(),
		const bool canonical = true,
		const bool row_wise = true,
		const bool generic = false
		);
//{
//	std::vector<unsigned long> orders(mat.NumRows(),order);
//	return approximant_basis(appbas,mat,orders,shift,canonical,row_wise,generic);
//}

/*------------------------------------------------------------*/
/* Iterative algorithm for general order and shift            */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo for any shift)        */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov)                   */
/*------------------------------------------------------------*/
std::vector<long> appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & mat,
		const long ord,
		const std::vector<long> & shift
		);

std::vector<long> popov_iter_appbas(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & mat,
		const long ord,
		const std::vector<long> & shift
		);

/*------------------------------------------------------------*/
/* M-Basis algorithm for approximant order = 1                */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo for any shift)        */
/*   - Jeannerod-Neiger-Villard 2018 (ensuring s-Popov)       */
/*------------------------------------------------------------*/
std::vector<long> mbasis1(
		Mat<zz_pX> & appbas,
		const Mat<zz_p> & mat,
		const std::vector<long> & shift
		);

/*------------------------------------------------------------*/
/* M-Basis algorithm for uniform approximant order            */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo for any shift)        */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov)                   */
/*------------------------------------------------------------*/
std::vector<long> mbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & mat,
		const long order,
		const std::vector<long> & shift
		);

std::vector<long> popov_mbasis(
		Mat<zz_pX> &appbas,
		const Mat<zz_pX> & mat,
		const long order,
		const std::vector<long> & shift
		);


#endif // MAT_LZZ_PX_EXTRA__H
