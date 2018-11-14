#ifndef MAT_LZZ_PX_UTILS__H
#define MAT_LZZ_PX_UTILS__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <vector>

NTL_CLIENT

/*------------------------------------------------------------*/
/* Basic routines for handling Mat<zz_pX>                     */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ZERO MATRIX; IDENTITY MATRIX                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* clears the matrix  (pmat = 0 with same dimensions)         */
/*------------------------------------------------------------*/
void clear(Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* set pmat to be the identity                                */
/* (same size, assuming square / size dim)                    */
/*------------------------------------------------------------*/
void ident(Mat<zz_pX> & pmat, long dim);

/*------------------------------------------------------------*/
/* return the identity matrix of size dim                     */
/*------------------------------------------------------------*/
Mat<zz_pX> ident_mat_zz_pX(long dim);

/*------------------------------------------------------------*/
/* tests whether m is zero (whatever its dims)                */
/*------------------------------------------------------------*/
long IsZero(const Vec<zz_pX> & pvec);
long IsZero(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix                  */
/*------------------------------------------------------------*/
long IsIdent(const Mat<zz_pX> & pmat);

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix of size 'dim'    */
/*------------------------------------------------------------*/
long IsIdent(const Mat<zz_pX> & pmat, long dim);




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MATRIX DEGREE                                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* maximum degree of the entries of pmat                      */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & pmat);
long deg(const Vec<zz_pX> & pvec); // TODO in its own file


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
{ Mat<zz_p> x; GetCoeff(x, a, i); return x; }

/*------------------------------------------------------------*/
/* sets ith coefficient of x to a                             */
/*------------------------------------------------------------*/
void SetCoeff(Mat<zz_pX>& x, long i, const Mat<zz_p> &a);





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TRANSPOSE, TRUNC, SHIFT, REVERSE, EVAL                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Transpose                                                  */
/*------------------------------------------------------------*/
void transpose(Mat<zz_pX>& x, const Mat<zz_pX>& a);

inline Mat<zz_pX> transpose(const Mat<zz_pX> & a)
{ Mat<zz_pX> x; transpose(x, a); return x; }

/*------------------------------------------------------------*/
/* Truncate mod X^..., for all the matrix / some columns/rows */
/* output can alias input                                     */
/*------------------------------------------------------------*/

// vector versions
void trunc(Vec<zz_pX>& x, const Vec<zz_pX>& a, long n);
inline Vec<zz_pX> trunc(const Vec<zz_pX>& a, long n)
{ Vec<zz_pX> x; trunc(x, a, n); return x; }

// full matrix versions
void trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
inline Mat<zz_pX> trunc(const Mat<zz_pX>& a, long n)
{ Mat<zz_pX> x; trunc(x, a, n); return x; }

// row versions
void truncRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, long r, long n);
inline Mat<zz_pX> truncRow(const Mat<zz_pX>& a, long r, long n)
{ Mat<zz_pX> x; truncRow(x, a, r, n); return x; }

// col versions
void truncCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, long c, long n);
inline Mat<zz_pX> truncCol(const Mat<zz_pX>& a, long c, long n)
{ Mat<zz_pX> x; truncCol(x, a, c, n); return x; }

/* TODO: different truncation orders on columns/rows          */

/*------------------------------------------------------------*/
/* Shift operations:                                          */
/*  - LeftShift by n means multiplication by X^n              */
/*  - RightShift by n means division by X^n                   */
/*  - a negative shift reverses the direction of the shift.   */
/*------------------------------------------------------------*/

/* TODO                                                       */
/* versions with different shifting orders on rows/columns    */
/* shiftAdd, shiftSub                                         */

// left shift, full matrix 
void LeftShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
inline Mat<zz_pX> LeftShift(const Mat<zz_pX>& a, long n)
{ Mat<zz_pX> x; LeftShift(x, a, n); return x; }

// right shift, full matrix 
void RightShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n);
inline Mat<zz_pX> RightShift(const Mat<zz_pX>& a, long n)
{ Mat<zz_pX> x; RightShift(x, a, n); return x; }

// left shift, single row
void LeftShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n);
inline Mat<zz_pX> LeftShiftRow(const Mat<zz_pX>& a, const long r, long n)
{ Mat<zz_pX> x; LeftShiftRow(x, a, r, n); return x; }

// right shifts, single row
void RightShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n);
inline Mat<zz_pX> RightShiftRow(const Mat<zz_pX>& a, const long r, long n)
{ Mat<zz_pX> x; RightShiftRow(x, a, r, n); return x; }

// left shift, single column
void LeftShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n);
inline Mat<zz_pX> LeftShiftCol(const Mat<zz_pX>& a, const long c, long n)
{ Mat<zz_pX> x; LeftShiftCol(x, a, c, n); return x; }

// right shifts, single column
void RightShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n);
inline Mat<zz_pX> RightShiftCol(const Mat<zz_pX>& a, const long c, long n)
{ Mat<zz_pX> x; RightShiftCol(x, a, c, n); return x; }

// operators for full matrix shift
inline Mat<zz_pX> operator<<(const Mat<zz_pX> &a, long n)
{ Mat<zz_pX> x; LeftShift(x, a, n); return x; }
inline Mat<zz_pX> operator>>(const Mat<zz_pX> &a, long n)
{ Mat<zz_pX> x; RightShift(x, a, n); return x; }

inline Mat<zz_pX>& operator<<=(Mat<zz_pX>& x, long n)
{ LeftShift(x, x, n); return x; }
inline Mat<zz_pX>& operator>>=(Mat<zz_pX>& x, long n)
{ RightShift(x, x, n); return x; }


/*------------------------------------------------------------*/
/* reverse the order of the entries in a vector               */
/* x = a[n - 1 -i], i=0..n-1, with n=length(a)                */
/*------------------------------------------------------------*/
void reverse_vector(Vec<zz_pX>& x, const Vec<zz_pX>& a);
inline Vec<zz_pX> reverse_vector(const Vec<zz_pX>& a)
{ Vec<zz_pX> x; reverse_vector(x, a); return x; }

/*------------------------------------------------------------*/
/* Reverse operations:                                        */
/* x = reverse of a[0]..a[hi] (hi >= -1);                     */
/* user provided 'hi'                                         */
/*------------------------------------------------------------*/
void reverse(Mat<zz_pX>& x, const Mat<zz_pX>& a, long hi);

inline Mat<zz_pX> reverse(const Mat<zz_pX>& a, long hi)
{ Mat<zz_pX> x; reverse(x, a, hi); return x; }

/* TODO versions with different degree on different cols/rows */
void reverse(
             Mat<zz_pX> &x, 
             const Mat<zz_pX> &a, 
             const VecLong &hi,
             const bool row_wise = true
            );

/*------------------------------------------------------------*/
/* Reverse operations:                                        */
/* x = reverse of a[0]..a[deg(a)]                             */
/*------------------------------------------------------------*/

inline void reverse(Mat<zz_pX>& x, const Mat<zz_pX>& a)
{ reverse(x, a, deg(a)); }

inline Mat<zz_pX> reverse(const Mat<zz_pX>& a)
{ Mat<zz_pX> x; reverse(x, a, deg(a)); return x; }


/*------------------------------------------------------------*/
/* Evaluate at one point                                      */
/*------------------------------------------------------------*/
void eval(Mat<zz_p> & evmat, const Mat<zz_pX> & pmat, zz_p pt);

inline Mat<zz_p> eval(const Mat<zz_pX> & pmat, zz_p pt)
{ Mat<zz_p> evmat; eval(evmat, pmat, pt); return evmat; }





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CREATE RANDOM MATRICES                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* random vector of length n, degree < d                      */
/*------------------------------------------------------------*/
void random(Vec<zz_pX> & pvec, long n, long d);
inline Vec<zz_pX> random_vec_zz_pX(long n, long d)
{ Vec<zz_pX> pvec; random(pvec, n, d); return pvec; }

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d                         */
/*------------------------------------------------------------*/
void random(Mat<zz_pX> & pmat, long m, long n, long d);
inline Mat<zz_pX> random_mat_zz_pX(long n, long m, long d)
{ Mat<zz_pX> pmat; random(pmat, n, m, d); return pmat; }

/*------------------------------------------------------------*/
/* random (m, n) matrix of row degree < rdeg                  */
/*------------------------------------------------------------*/
void random_mat_zz_pX_rdeg(Mat<zz_pX> & pmat, long m, long n, VecLong rdeg);

/*------------------------------------------------------------*/
/* random (m, n) matrix of column degree < cdeg               */
/*------------------------------------------------------------*/
void random_mat_zz_pX_cdeg(Mat<zz_pX> & pmat, long m, long n, VecLong cdeg);
// TODO replace with Vec<long> ??





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CONVERT                                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* convert from Mat<zz_p>                                     */
/*------------------------------------------------------------*/
void conv(Mat<zz_pX>& mat, const Mat<zz_p>& coeff);

inline Mat<zz_pX> conv(const Mat<zz_p>& coeff)
{ Mat<zz_pX> mat; conv(mat, coeff); return mat; }

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (degree deduced from input)                                */
/*------------------------------------------------------------*/
void conv(Vec<Mat<zz_p>>& coeffs, const Mat<zz_pX>& mat);

inline Vec<Mat<zz_p>> conv(const Mat<zz_pX>& mat)
{ Vec<Mat<zz_p>> coeffs; conv(coeffs, mat); return coeffs; }

void conv(Mat<zz_pX>& mat, const Vec<Mat<zz_p>>& coeffs);

inline Mat<zz_pX> conv(const Vec<Mat<zz_p>>& coeffs)
{ Mat<zz_pX> mat; conv(mat, coeffs); return mat; }

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (user provided truncation order)                           */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* coeffs will have length order independently of deg(mat)    */
/*------------------------------------------------------------*/
void conv(Vec<Mat<zz_p>>& coeffs, const Mat<zz_pX>& mat, const long order);

inline Vec<Mat<zz_p>> conv(const Mat<zz_pX>& mat, const long order)
{ Vec<Mat<zz_p>> coeffs; conv(coeffs, mat, order); return coeffs; }

void conv(Mat<zz_pX>& mat, const Vec<Mat<zz_p>>& coeffs, const long order);

inline Mat<zz_pX> conv(const Vec<Mat<zz_p>>& coeffs, const long order)
{ Mat<zz_pX> mat; conv(mat, coeffs, order); return mat; }

#endif // MAT_LZZ_PX_UTILS__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
