#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* random (n, m) matrix of degree < d                         */
/*------------------------------------------------------------*/
void random_mat_zz_pX(Mat<zz_pX>& a, long n, long m, long d);

/*------------------------------------------------------------*/
/* maximum degree of the entries of a                         */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & a);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic arithmetic                                           */
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
/* negate                                                     */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* get / set coefficients                                     */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* convert from Mat<zz_p>                                     */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/*------------------------------------------------------------*/



/*------------------------------------------------------------*/
/* leading and trailing coefficient matrices                  */
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
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

#endif
