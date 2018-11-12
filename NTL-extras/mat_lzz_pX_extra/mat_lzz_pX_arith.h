#ifndef MAT_LZZ_PX_ARITH__H
#define MAT_LZZ_PX_ARITH__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                  BASIC ARITHMETIC                          */
/* addition, subtraction, negation,                           */
/* multiplication by constant / by single polynomial, ...     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* matrix addition                                            */
/*------------------------------------------------------------*/
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);

inline void add(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{ add(c, b, a); }

inline Mat<zz_pX> operator+(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> x; add(x, a, b); return x; }

inline Mat<zz_pX> operator+(const Mat<zz_pX>& a, const Mat<zz_p>& b)
{ Mat<zz_pX> x; add(x, a, b); return x; }

inline Mat<zz_pX> operator+(const Mat<zz_p>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> x; add(x, a, b); return x; }

inline Mat<zz_pX> & operator+=(Mat<zz_pX> & x, const Mat<zz_pX>& b)
{ add(x, x, b); return x; }

inline Mat<zz_pX> & operator+=(Mat<zz_pX> & x, const Mat<zz_p>& b)
{ add(x, x, b); return x; }


/*------------------------------------------------------------*/
/* vector addition                                            */
/*------------------------------------------------------------*/
// TODO in its own file for Vec<zz_pX>?
void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b);
void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b);

inline void add(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b)
{ add(c, b, a); }

inline Vec<zz_pX> operator+(const Vec<zz_pX>& a, const Vec<zz_pX>& b)
{ Vec<zz_pX> x; add(x, a, b); return x; }

inline Vec<zz_pX> operator+(const Vec<zz_pX>& a, const Vec<zz_p>& b)
{ Vec<zz_pX> x; add(x, a, b); return x; }

inline Vec<zz_pX> operator+(const Vec<zz_p>& a, const Vec<zz_pX>& b)
{ Vec<zz_pX> x; add(x, a, b); return x; }

inline Vec<zz_pX> & operator+=(Vec<zz_pX> & x, const Vec<zz_pX>& b)
{ add(x, x, b); return x; }

inline Vec<zz_pX> & operator+=(Vec<zz_pX> & x, const Vec<zz_p>& b)
{ add(x, x, b); return x; }



/*------------------------------------------------------------*/
/* matrix subtraction                                         */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);
void sub(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

inline Mat<zz_pX> operator-(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> x; sub(x, a, b); return x; }

inline Mat<zz_pX> operator-(const Mat<zz_pX>& a, const Mat<zz_p>& b)
{ Mat<zz_pX> x; sub(x, a, b); return x; }

inline Mat<zz_pX> operator-(const Mat<zz_p>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> x; sub(x, a, b); return x; }

inline Mat<zz_pX> & operator-=(Mat<zz_pX> & x, const Mat<zz_pX>& b)
{ sub(x, x, b); return x; }

inline Mat<zz_pX> & operator-=(Mat<zz_pX> & x, const Mat<zz_p>& b)
{ sub(x, x, b); return x; }


/*------------------------------------------------------------*/
/* vector subtraction                                         */
/*------------------------------------------------------------*/
void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b);
void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b);
void sub(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b);

inline Vec<zz_pX> operator-(const Vec<zz_pX>& a, const Vec<zz_pX>& b)
{ Vec<zz_pX> x; sub(x, a, b); return x; }

inline Vec<zz_pX> operator-(const Vec<zz_pX>& a, const Vec<zz_p>& b)
{ Vec<zz_pX> x; sub(x, a, b); return x; }

inline Vec<zz_pX> operator-(const Vec<zz_p>& a, const Vec<zz_pX>& b)
{ Vec<zz_pX> x; sub(x, a, b); return x; }

inline Vec<zz_pX> & operator-=(Vec<zz_pX> & x, const Vec<zz_pX>& b)
{ sub(x, x, b); return x; }

inline Vec<zz_pX> & operator-=(Vec<zz_pX> & x, const Vec<zz_p>& b)
{ sub(x, x, b); return x; }


/*------------------------------------------------------------*/
/* multiplication by a constant matrix                        */
/*------------------------------------------------------------*/
// TODO: for mul with rhs constant
// -- it seems more efficient to rather expand 'a' as single big constant
// matrix 'cmat', and compute b*cmat, and retrieve back the entries in 'c'
// -- this would require first computing cdeg(a); it may be given by the user
// as an optional parameter since it is sometimes known from pivdeg and such
// ---->>> make this depend on conv?

// TODO: for mul with lhs constant
// -- it seems more efficient to rather expand 'a' as single big constant
// matrix 'cmat', and compute b*cmat, and retrieve back the entries in 'c'
// -- this would require first computing cdeg(a); it may be given by the user
// as an optional parameter since it is sometimes known from pivdeg and such


// TODO mul/multiply? unify names
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b);
void mul(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b);

inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const Mat<zz_p>& b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

inline Mat<zz_pX> operator*(const Mat<zz_p>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

inline Mat<zz_pX> & operator*=(Mat<zz_pX> & x, const Mat<zz_p>& b)
{ mul(x, x, b); return x; }

/*------------------------------------------------------------*/
/* scalar multiplication for vectors                          */
/*------------------------------------------------------------*/
void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_p & b);

inline void mul(Vec<zz_pX> & c, const zz_p & a, const Vec<zz_pX> & b)
{ mul(c, b, a); }

inline Vec<zz_pX> operator*(const Vec<zz_pX>& a, const zz_p& b)
{ Vec<zz_pX> x; mul(x, a, b); return x; }

inline Vec<zz_pX> operator*(const zz_p& a, const Vec<zz_pX>& b)
{ Vec<zz_pX> x; mul(x, a, b); return x; }

/*------------------------------------------------------------*/
/* scalar multiplication                                      */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b);

inline void mul(Mat<zz_pX> & c, const zz_p & a, const Mat<zz_pX> & b)
{ mul(c, b, a); }

inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const zz_p& b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

inline Mat<zz_pX> operator*(const zz_p& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> x; mul(x, a, b); return x; }

/*------------------------------------------------------------*/
/* polynomial multiplication                                  */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_pX & b);

inline void mul(Mat<zz_pX> & c, const zz_pX & a, const Mat<zz_pX> & b)
{
    mul(c, b, a);
}

inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const zz_pX& b)
{ 
    Mat<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}

inline Mat<zz_pX> operator*(const zz_pX& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}

/*------------------------------------------------------------*/
/* polynomial multiplication for vectors                      */
/*------------------------------------------------------------*/
void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_pX & b);

inline void mul(Vec<zz_pX> & c, const zz_pX & a, const Vec<zz_pX> & b)
{
    mul(c, b, a);
}

inline Vec<zz_pX> operator*(const Vec<zz_pX>& a, const zz_pX& b)
{ 
    Vec<zz_pX> x; 
    mul(x, a, b); 
    return x; 
}

inline Vec<zz_pX> operator*(const zz_pX& a, const Vec<zz_pX>& b)
{ 
    Vec<zz_pX> x; 
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

#endif /* ifndef MAT_LZZ_PX_ARITH__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
