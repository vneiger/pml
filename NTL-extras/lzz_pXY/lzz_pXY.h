#ifndef LZZ_PXY__H
#define LZZ_PXY__H

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over zz_p         */
/* a zz_pXY is simply a vector of zz_pX                       */
/* with the convention f = sum_i rep[i](X) Y^i                */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pXY
{
public:
    /*------------------------------------------------------------*/
    /* array of coefficients                                      */
    /*------------------------------------------------------------*/
    Vec<zz_pX> rep;

    /*------------------------------------------------------------*/
    /* zero test                                                  */
    /*------------------------------------------------------------*/
    long is_zero() const;

    /*------------------------------------------------------------*/
    /* set to zero                                                */
    /*------------------------------------------------------------*/
    void zero();

    /*------------------------------------------------------------*/
    /* total degree                                               */
    /*------------------------------------------------------------*/
    long tdeg() const;

    /*------------------------------------------------------------*/
    /* degree in X                                                */
    /*------------------------------------------------------------*/
    long degX() const;
    
    /*------------------------------------------------------------*/
    /* degree in Y                                                */
    /*------------------------------------------------------------*/
    long degY() const;
    
    /*------------------------------------------------------------*/
    /* resizes the array of coefficients                          */
    /* to remove the trailing entries that are zero, if any       */
    /*------------------------------------------------------------*/
    void normalize();
    
    /*------------------------------------------------------------*/
    /* plain constructor                                          */
    /*------------------------------------------------------------*/
    zz_pXY(){}
    
    /*------------------------------------------------------------*/
    /* builds from a vector of zz_pX                              */
    /*------------------------------------------------------------*/
    zz_pXY(const Vec<zz_pX> & rep);

    /*------------------------------------------------------------*/
    /* builds from a zz_pX                                        */
    /*------------------------------------------------------------*/
    zz_pXY(const zz_pX & rep);
    
    /*------------------------------------------------------------*/
    /* plain destructor                                           */
    /*------------------------------------------------------------*/
    ~zz_pXY(){}
};

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
ostream& operator<<(ostream& s, const zz_pXY& a);

/*------------------------------------------------------------*/
/* random element f with deg(f,x) < dx and deg(f,y) < dy      */
/*------------------------------------------------------------*/
void random_zz_pXY(zz_pXY & f, long dx, long dy);

inline zz_pXY random_zz_pXY(long dx, long dy)
{
    zz_pXY c;
    random_zz_pXY(c, dx, dy);
    return c;
}

/*------------------------------------------------------------*/
/* evaluation with respect to X at a point                    */
/*------------------------------------------------------------*/
void evaluate(zz_pX & value, const zz_pXY & f, const zz_p & point);

inline zz_pX evaluate(const zz_pXY & f, const zz_p & point)
{
    zz_pX c;
    evaluate(c, f, point);
    return c;
}

/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
void add(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

inline zz_pXY add(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    add(c, a, b);
    return c;
}

inline zz_pXY operator+(const zz_pXY& a, const zz_pXY& b)
{
    return add(a, b);
}

/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
void sub(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

inline zz_pXY sub(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    sub(c, a, b);
    return c;
}

inline zz_pXY operator-(const zz_pXY& a, const zz_pXY& b)
{
    return sub(a, b);
}


/*------------------------------------------------------------*/
/* naive multiplication                                       */
/*------------------------------------------------------------*/
void mul_naive(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

inline zz_pXY mul_naive(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    mul_naive(c, a, b);
    return c;
}


/*------------------------------------------------------------*/
/* kronecker substitution using vector arguments              */
/*------------------------------------------------------------*/
void to_kronecker(zz_pX& out, const Vec<zz_pX>& a, long stepX);
void from_kronecker(Vec<zz_pX>& out, const zz_pX& a, long stepX);

/*------------------------------------------------------------*/
/* kronecker substitution                                     */
/*------------------------------------------------------------*/
void to_kronecker(zz_pX& out, const zz_pXY& a, long stepX);
void from_kronecker(zz_pXY& out, const zz_pX& a, long stepX);

/*------------------------------------------------------------*/
/* kronecker multiplication                                   */
/*------------------------------------------------------------*/
void mul_kronecker(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

inline zz_pXY mul_kronecker(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    mul_kronecker(c, a, b);
    return c;
}

/*------------------------------------------------------------*/
/* multiplication = kronecker for the time being              */
/*------------------------------------------------------------*/
inline void mul(zz_pXY& c, const zz_pXY& a, const zz_pXY& b)
{
    mul_kronecker(c, a, b);
}

inline zz_pXY mul(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    mul(c, a, b);
    return c;
}

inline zz_pXY operator*(const zz_pXY& a, const zz_pXY& b)
{
    return mul(a, b);
}

/*------------------------------------------------------------*/
/* equality test                                              */
/*------------------------------------------------------------*/
inline long operator==(const zz_pXY& a, const zz_pXY& b)
{
   return a.rep == b.rep;
}

inline long operator!=(const zz_pXY& a, const zz_pXY& b)
{ 
    return !(a == b); 
}

/*------------------------------------------------------------*/
/* shift in x.                                                */ 
/* i >=0 means multiply by x^i                                */
/* i <=0 means div by x^(-i)                                  */
/*------------------------------------------------------------*/
void shift_x(zz_pXY& b, const zz_pXY &a, long i);

inline zz_pXY shift_x(const zz_pXY &a, long i)
{
    zz_pXY b;
    shift_x(b, a, i);
    return b;
}

/*------------------------------------------------------------*/
/* shift in y.                                                */ 
/* i >=0 means multiply by y^i                                */
/* i <=0 means div by y^(-i)                                  */
/*------------------------------------------------------------*/
void shift_y(zz_pXY& b, const zz_pXY &a, long i);

inline zz_pXY shift_y(const zz_pXY &a, long i)
{
    zz_pXY b;
    shift_y(b, a, i);
    return b;
}

/*------------------------------------------------------------*/
/* trunc in x (modulo x^i)                                    */ 
/*------------------------------------------------------------*/
void trunc_x(zz_pXY& b, const zz_pXY &a, long i);

inline zz_pXY trunc_x(const zz_pXY &a, long i)
{
    zz_pXY b;
    trunc_x(b, a, i);
    return b;
}

/*------------------------------------------------------------*/
/* trunc in y (modulo y^i)                                    */ 
/*------------------------------------------------------------*/
void trunc_y(zz_pXY& b, const zz_pXY &a, long i);

inline zz_pXY trunc_y(const zz_pXY &a, long i)
{
    zz_pXY b;
    trunc_y(b, a, i);
    return b;
}

/*------------------------------------------------------------*/
/* setting coefficients                                       */
/*------------------------------------------------------------*/
void SetCoeff(zz_pXY& a, long i, const zz_pX& c);


/*------------------------------------------------------------*/
/* computes the resultant of f and g                          */
/* returns 0 iff the computation failed                       */
/*------------------------------------------------------------*/
long resultant(zz_pX& res, const zz_pXY& f, const zz_pXY& g);


#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
