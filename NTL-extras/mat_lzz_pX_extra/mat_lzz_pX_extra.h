#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector> // std vector, for shifts, degrees, pivot indices

#include "lzz_p_extra.h"
#include "thresholds_matrix_multiply.h"
#include "thresholds_matrix_middle_product.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* UTILS                                                      */
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
/* Types for integer tuples: degrees and shifts               */
/*------------------------------------------------------------*/

typedef std::vector<long> Shift;
typedef std::vector<long> DegVec;

// TODO : kill? swap? vec_lzz_pX ?

/*------------------------------------------------------------*/
/* clears the matrix  (pmat = 0 with same dimensions)         */
/*------------------------------------------------------------*/
void clear(Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* tests whether pmat is the zero matrix (whatever its dims)  */
/*------------------------------------------------------------*/
long IsZero(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* set pmat to be the identity                                */
/* (same size, assuming square / size dim)                    */
/*------------------------------------------------------------*/
void set(Mat<zz_pX> & pmat);
void set(Mat<zz_pX> & pmat, long dim);

/*------------------------------------------------------------*/
/* return the identity matrix of size dim                     */
/*------------------------------------------------------------*/
Mat<zz_pX> identity(long dim);

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix                  */
/*------------------------------------------------------------*/
long IsIdentity(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* maximum degree of the entries of pmat                      */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* evaluate at a given point                                  */
/*------------------------------------------------------------*/
void eval(Mat<zz_p> & evmat, const Mat<zz_pX> & pmat, zz_p pt);

inline Mat<zz_p> eval(const Mat<zz_pX> & pmat, zz_p pt)
{
    Mat<zz_p> evmat;
    eval(evmat, pmat, pt);
    return evmat;
}


/*------------------------------------------------------------*/
/* transpose                                                  */
/*------------------------------------------------------------*/
void transpose(Mat<zz_pX>& x, const Mat<zz_pX>& a);

inline Mat<zz_pX> transpose(const Mat<zz_pX> & a)
{ 
    Mat<zz_pX> x; 
    transpose(x, a); 
    return x; 
}



/*------------------------------------------------------------*/
/* truncate mod X^..., for all the matrix / some columns/rows */
/* TODO: different truncation orders on columns/rows          */
/*------------------------------------------------------------*/
// full matrix versions
void trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
Mat<zz_pX> trunc(const Mat<zz_pX>& a, long n);

// row versions
void truncRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, long r, long n);
Mat<zz_pX> truncRow(const Mat<zz_pX>& a, long r, long n);

// col versions
void truncCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, long c, long n);
Mat<zz_pX> truncCol(const Mat<zz_pX>& a, long c, long n);


/*------------------------------------------------------------*/
/*                        Shift Operations                    */
/* LeftShift by n means multiplication by X^n                 */
/* RightShift by n means division by X^n                      */
/* A negative shift reverses the direction of the shift.      */
/* TODO                                                       */
/* versions with different shifting orders on rows/columns    */
/*------------------------------------------------------------*/
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


/*------------------------------------------------------------*/
/* reverse operations                                         */
/* x = reverse of a[0]..a[hi] (hi >= -1);                     */
/* hi defaults to deg(a) in second version                    */
/* TODO versions with different shifts on different cols/rows */
/*------------------------------------------------------------*/
void reverse(Mat<zz_pX>& x, const Mat<zz_pX>& a, long hi);

static inline Mat<zz_pX> reverse(const Mat<zz_pX>& a, long hi)
{
    Mat<zz_pX> x;
    reverse(x, a, hi);
    return x;
}

static inline void reverse(Mat<zz_pX>& x, const Mat<zz_pX>& a)
{
    reverse(x, a, deg(a));
}

static inline Mat<zz_pX> reverse(const Mat<zz_pX>& a)
{
    Mat<zz_pX> x;
    reverse(x, a);
    return x;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d                         */
/*------------------------------------------------------------*/
void random_mat_zz_pX(Mat<zz_pX>& pmat, long m, long n, long d);

/*------------------------------------------------------------*/
/* random (m, n) matrix of row degree < rdeg                  */
/*------------------------------------------------------------*/
void random_mat_zz_pX_rdeg(Mat<zz_pX>& pmat, long m, long n, DegVec rdeg);

/*------------------------------------------------------------*/
/* random (m, n) matrix of column degree < cdeg               */
/*------------------------------------------------------------*/
void random_mat_zz_pX_cdeg(Mat<zz_pX>& pmat, long m, long n, DegVec cdeg);

// TODO random matrix with given PolMatForm

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
/* TODO                                                       */
/* multiply row or column of matrix (vec_lzz_pX) by constant  */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
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
void SetCoeff(Mat<zz_pX>& x, long i, const Mat<zz_p> &a);

/*------------------------------------------------------------*/
/* convert from Mat<zz_p>                                     */
/*------------------------------------------------------------*/
void conv(Mat<zz_pX>& mat, const Mat<zz_p>& coeff);

inline Mat<zz_pX> conv(const Mat<zz_p>& coeff)
{
    Mat<zz_pX> mat;
    conv(mat, coeff);
    return mat;
}

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (degree deduced from input)                                */
/*------------------------------------------------------------*/
void conv(Vec<Mat<zz_p>>& coeffs, const Mat<zz_pX>& mat);

inline Vec<Mat<zz_p>> conv(const Mat<zz_pX>& mat)
{
    Vec<Mat<zz_p>> coeffs;
    conv(coeffs, mat);
    return coeffs;
}

void conv(Mat<zz_pX>& mat, const Vec<Mat<zz_p>>& coeffs);

inline Mat<zz_pX> conv(const Vec<Mat<zz_p>>& coeffs)
{
    Mat<zz_pX> mat;
    conv(mat, coeffs);
    return mat;
}

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (user provided truncation order)                           */
/*------------------------------------------------------------*/
void conv(Vec<Mat<zz_p>>& coeffs, const Mat<zz_pX>& mat, const long order);

inline Vec<Mat<zz_p>> conv(const Mat<zz_pX>& mat, const long order)
{
    Vec<Mat<zz_p>> coeffs;
    conv(coeffs, mat, order);
    return coeffs;
}

void conv(Mat<zz_pX>& mat, const Vec<Mat<zz_p>>& coeffs, const long order);

inline Mat<zz_pX> conv(const Vec<Mat<zz_p>>& coeffs, const long order)
{
    Mat<zz_pX> mat;
    conv(mat, coeffs, order);
    return mat;
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_geometric_using_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_geometric_no_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_montgomery3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba4(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long len);

inline void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    multiply_transform(c, a, b, max(deg(a), deg(b)) + 1);
}

void multiply_3_primes(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

void multiply(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long is_prime = 1);

inline void mul_trunc(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long n, long is_prime = 1)
{
    multiply(c, a, b, is_prime);
    trunc(c, c, n);
}

// inline Mat<zz_pX> & operator*=(Mat<zz_pX> & x, const Mat<zz_pX>& b)
// {
//     mul(x, x, b);
//     return x;
// }

/*------------------------------------------------------------*/
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* a must have degree at most dA                              */
/* c must have degree at most dA + dB                         */
/*------------------------------------------------------------*/
void t_multiply_evaluate_geometric(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/*------------------------------------------------------------*/
void middle_product_naive(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_FFT(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_3_primes(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_evaluate(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, long is_prime = 1);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* DEGREES, PIVOTS, LEADING MATRIX                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Degree matrix: matrix of the degree of each entry          */
/* Convention: deg(0) = -1, more generally the shifted degree */
/* of a zero entry is min(shift)-1                            */
/*------------------------------------------------------------*/
void degree_matrix(
                   Mat<long> &degmat,
                   const Mat<zz_pX> &pmat,
                   const Shift & shift = Shift(),
                   const bool row_wise=true
                  );

inline Mat<long> degree_matrix(
                               const Mat<zz_pX> &pmat,
                               const Shift & shift = Shift(),
                               const bool row_wise=true
                              )
{
    Mat<long> degmat;
    degree_matrix(degmat,pmat,shift,row_wise);
    return degmat;
}

/*------------------------------------------------------------*/
/* tuple (shifted deg row 1, shifted deg row 2, ...)          */
/* where shifted deg row k is the maximum of the shifted      */
/* degrees of the entries of row k                            */
/*------------------------------------------------------------*/
void row_degree(
                DegVec & rdeg,
                const Mat<zz_pX> &pmat,
                const Shift & shift = Shift()
               ); 

/*------------------------------------------------------------*/
/* similar function for column degrees (see row_degree)       */
/*------------------------------------------------------------*/
void col_degree(
                DegVec & cdeg,
                const Mat<zz_pX> &pmat,
                const Shift & shift = Shift()
               ); 

/*------------------------------------------------------------*/
/* similar function with row-wise option and returning degree */
/*------------------------------------------------------------*/
inline DegVec vector_degree(
                            const Mat<zz_pX> &pmat,
                            const Shift & shift = Shift(),
                            const bool row_wise = true
                           )
{
    DegVec degs;
    if (row_wise)
    {
        degs.resize(pmat.NumRows());
        row_degree(degs,pmat,shift);
    }
    else
    {
        degs.resize(pmat.NumCols());
        col_degree(degs,pmat,shift);
    }
    return degs;
}


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
/* leading matrix of pmat                                     */
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
//  Order orders(mat.NumRows(),order);
//  return approximant_basis(appbas,mat,orders,shift,canonical,row_wise,generic);
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

// TODO check if serious difference of time if not returning Popov but just
// minimal, like done in GJV03 and GL14 (implies slightly less permutation
// work: the final permutation of the rows is not necessary)

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

// variant which first converts to vector of constant matrices
// TODO see if this is ever slower than the above
DegVec mbasis_vector(
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
//  std::vector<std::pair<zz_p,long>> list_pts(pts.size());
//  for ( long i=0; i<list_pts.size(); ++i )
//  {
//    list_pts[i] = std::pair<zz_p,long>(pts[i],1);
//  }
//  Points points(mat.NumCols(),list_pts);
//  return interpolant_basis(appbas,mat,points,shift,canonical,row_wise,generic);
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

// TODO (naive version written)
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & pmat, // vector of evaluations
                          const Vec<zz_p> & pts, // "uniform" case
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
                        bool point_wise=true // TODO to be thought about
                       );

DegVec popov_intbas_iterative(
                              Mat<zz_pX> & intbas,
                              const Mat<zz_pX> & pmat,
                              const Points & pts,
                              const Shift & shift,
                              bool point_wise=true // TODO to be thought about
                             );


/*------------------------------------------------------------*/
/* Adaptation of M-Basis for uniform interpolation points     */
/*------------------------------------------------------------*/

// --> popov_mbasis1 can be called as such (with, as input, pmat evaluated at a
// point)

// TODO input pmat = polynomial matrix, not implemented yet
DegVec mbasis(
              Mat<zz_pX> & intbas,
              const Mat<zz_pX> & pmat,
              const Vec<zz_p> & pts,
              const Shift & shift
             );

// input pmat = list of evaluations, implemented
DegVec mbasis(
              Mat<zz_pX> & intbas,
              const Vec<Mat<zz_p>> & evals,
              const Vec<zz_p> & pts,
              const Shift & shift
             );

DegVec popov_mbasis(
                    Mat<zz_pX> &intbas,
                    const Mat<zz_pX> & pmat,
                    const Vec<zz_p> & pts,
                    const Shift & shift
                   );

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform interpolation points        */
/*------------------------------------------------------------*/

// TODO there are two variants, test them to be sure if they are similar / which is faster
//   either compute more in the evaluated world and interpolate intbas at the end,
//   or compute in the polynomial world and evaluate to obtain the residuals
// (in any case, there will still be interpolation/evaluation in the middle)
DegVec pmbasis(
               Mat<zz_pX> & intbas,
               const Mat<zz_pX> & pmat,
               const Vec<zz_p> & pts,
               const Shift & shift
              );

DegVec popov_pmbasis(
                     Mat<zz_pX> &intbas,
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

/**********************************************************************
 * plain truncated inverse                                            *
 * Requires: a(0) is invertible                                       *
 **********************************************************************/
void plain_inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m);
void newton_inv_trunc_FFT(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m);

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

// TODO diagonal entries of Hermite form, not implemented yet
void diagonal_of_hermite(Vec<zz_pX> & diag, const Mat<zz_pX> & pmat);

/***********************************************
 *  Labahn-Neiger-Zhou Hermite form algorithm  *
 ***********************************************/

// TODO (requires partial linearization + basis reduction)


/**********************************************************************
 *                       DETERMINANT ALGORITHMS                       *
 **********************************************************************/

// general user interface
// TODO (not implemented yet)
void determinant(zz_pX & det, const Mat<zz_pX> & pmat);

inline zz_pX determinant(const Mat<zz_pX> & pmat)
{
    zz_pX det;
    determinant(det, pmat);
    return det;
}

// verifies that det = c det(pmat),
// for some c a nonzero field element if up_to_constant==false; and c=1 otherwise
// if randomized==true, it is allowed to use a Monte Carlo randomized approach
// TODO: only randomized implemented for now
bool verify_determinant(const zz_pX & det, const Mat<zz_pX> & pmat, bool up_to_constant, bool randomized);


/*******************************************************************
 *  Labahn-Neiger-Zhou: via diagonal entries of triangularization  *
 *******************************************************************/

void determinant_via_diagonal_of_hermite(zz_pX & det, const Mat<zz_pX> & pmat);

// TODO version assuming diagonal of (lower triangular) Hermite form is
// [det 1 .. 1] and that we have generic degree profiles in the partial
// triangularization of Labahn-Neiger-Zhou

// Version 1 (deterministic algorithm): degree of determinant is known (e.g. if
// matrix is reduced or if the determinant corresponds to some invariant of an
// object with known degree). Runs the partial triangularization and returns
// false if the computed diagonal entry of Hermite form does not have the right
// degree. True is returned iff determinant is correct. 
// TODO currently returns determinant up to constant factor!!
bool determinant_generic_knowing_degree(zz_pX & det, const Mat<zz_pX> & pmat, long degree);

// Version 2 (Las Vegas randomized algorithm): runs the partial
// triangularization and checks determinant by Zippel-Schwartz; returns false
// if determinant is wrong or if field size is too small for allowing the
// Zippel-Schwartz check. If true is returned, then determinant is correct.
// TODO: not implemented yet
bool determinant_generic_las_vegas(zz_pX & det, const Mat<zz_pX> & pmat);

// TODO other determinant algorithms??
// --> could rely on x-Smith decomposition of Gupta et al (worth
// implementing??), cf Appendix of LaNeZh17 
// --> no randomized faster approach?


#endif // MAT_LZZ_PX_EXTRA__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
