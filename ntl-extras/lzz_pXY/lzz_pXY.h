#ifndef LZZ_PXY__H
#define LZZ_PXY__H

/** \brief Basic arithmetic for bivariate polynomials over `zz_p`
 *
 * \file lzz_pXY.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-02-01
 *
 * This file contains the declarations a helper class performing basic
 * arithmetic with bivariate polynomials over `zz_p`. A `zz_pXY` is simply 
 * a vector of zz_pX, with the convention f = sum_i rep[i](X) Y^i.
 * Minimal functionalities are provided (add, sub, multiplication),
 * together with (experimental) implementations of bivariate resultants
 *
 */

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

NTL_CLIENT

/**  Main class for bivariate polynomials.
 *   
 */
class zz_pXY
{
public:
    /*------------------------------------------------------------*/
    /* array of coefficients                                      */
    /*------------------------------------------------------------*/
    /** Array of coefficients. The last entry is maintained to be non-zero,
     *  unless `this` is zero.
     */
    Vec<zz_pX> rep;

    /*------------------------------------------------------------*/
    /* zero test                                                  */
    /*------------------------------------------------------------*/
    /** Returns 1 if `this` is zero, 0 otherwise */
    long is_zero() const;

    /*------------------------------------------------------------*/
    /* set to zero                                                */
    /*------------------------------------------------------------*/
    /** Sets `this` to zero by emptying `rep` */
    void zero();

    /*------------------------------------------------------------*/
    /* total degree                                               */
    /*------------------------------------------------------------*/
    /** Returns the total degree of `this` */
    long tdeg() const;

    /*------------------------------------------------------------*/
    /* degree in X                                                */
    /*------------------------------------------------------------*/
    /** Returns the X-degree of `this` */
    long degX() const;
    
    /*------------------------------------------------------------*/
    /* degree in Y                                                */
    /*------------------------------------------------------------*/
    /** Returns the Y-degree of `this`, which is the length of
     *  `rep` minus one.
     */
    long degY() const;
    
    /*------------------------------------------------------------*/
    /* resizes the array of coefficients                          */
    /* to remove the trailing entries that are zero, if any       */
    /*------------------------------------------------------------*/
    /** Resizes the array of coefficients, by removing the trailing   
     *  zero entries, if any.
     */
    void normalize();
    
    /*------------------------------------------------------------*/
    /* plain constructor                                          */
    /*------------------------------------------------------------*/
    /** Creates the zero bivariate polynomial */
    zz_pXY(){}
    
    /*------------------------------------------------------------*/
    /* builds from a vector of zz_pX                              */
    /*------------------------------------------------------------*/
    /** Builds from a vector of `zz_pX` */
    zz_pXY(const Vec<zz_pX> & rep);

    /*------------------------------------------------------------*/
    /* builds from a zz_pX                                        */
    /*------------------------------------------------------------*/
    /** Builds from a single `zz_pX` */
    zz_pXY(const zz_pX & rep);
    
    /*------------------------------------------------------------*/
    /* plain destructor                                           */
    /*------------------------------------------------------------*/
    /** Destructor */
    ~zz_pXY(){}
};


/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
/** @name I/O
 *  @{
 */ 

/** Prints `a` to `s` using NTL's representation */
ostream& operator<<(ostream& s, const zz_pXY& a);

/** @}  */       // doxygen I/O


/*------------------------------------------------------------*/
/* random element f with deg(f,x) < dx and deg(f,y) < dy      */
/*------------------------------------------------------------*/
/** @name Random generation
 *  @{
 */ 

/** Sets `f` as a random `zz_pXY` with deg(`f`,X) < `dx` and deg(`f`,Y) < `dy` */ 
void random_zz_pXY(zz_pXY & f, long dx, long dy);

/** Returns a random `zz_pXY` with deg(`f`,X) < `dx` and deg(`f`,Y) < `dy` */ 
inline zz_pXY random_zz_pXY(long dx, long dy)
{
    zz_pXY c;
    random_zz_pXY(c, dx, dy);
    return c;
}
/** @}  */       // doxygen random


/*------------------------------------------------------------*/
/* evaluation with respect to X at a point                    */
/*------------------------------------------------------------*/
/** @name Evaluation 
 *  @{
 */ 

/** 
 * Sets `value` as `f`(`point`, `x`).
 * `value` is a `zz_pX` 
 */ 
void evaluate(zz_pX & value, const zz_pXY & f, const zz_p & point);

/** Returns  `f`(`point`, `x`) as a `zz_px` */ 
inline zz_pX evaluate(const zz_pXY & f, const zz_p & point)
{
    zz_pX c;
    evaluate(c, f, point);
    return c;
}
/** @}  */       // doxygen evaluation


/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
/** @name Addition
 *  @{
 */ 

/** Sets `c` = `a` + `b` */
void add(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

/** Returns `a` + `b` */
inline zz_pXY add(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    add(c, a, b);
    return c;
}

/** Returns `a` + `b` */
inline zz_pXY operator+(const zz_pXY& a, const zz_pXY& b)
{
    return add(a, b);
}
/** @}  */       // doxygen addition


/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
/** @name Subtraction
 *  @{
 */ 

/** Sets `c` = `a` - `b` */
void sub(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

/** Returns `a` - `b` */
inline zz_pXY sub(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    sub(c, a, b);
    return c;
}

/** Returns `a` - `b` */
inline zz_pXY operator-(const zz_pXY& a, const zz_pXY& b)
{
    return sub(a, b);
}
/** @}  */       // doxygen subtraction


/*------------------------------------------------------------*/
/* multiplication                                             */
/*------------------------------------------------------------*/
/** @name Multiplication and related functions
 *  @{
 */ 


/** Sets `c` = `a` * `b`, `b` is a constant. */
void mul(zz_pXY& c, const zz_pXY& a, const zz_p b);

/** Returns `a` * `b`, `b` is a constant. */
inline zz_pXY mul(const zz_pXY& a, const zz_p b)
{
    zz_pXY c;
    mul(c, a, b);
    return c;
}

/** Returns `a` * `b`, `b` is a constant. */
inline zz_pXY operator*(const zz_pXY& a, const zz_p b)
{
    return mul(a, b);
}

/** Sets `c` = `a` * `b`, computed using the naive algorithm */
void mul_naive(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

/** Returns `a` * `b`, computed using the naive algorithm */
inline zz_pXY mul_naive(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    mul_naive(c, a, b);
    return c;
}

/*------------------------------------------------------------*/
/* kronecker substitution                                     */
/*------------------------------------------------------------*/
/**  Sets `out` to be A(X, X^`stepX`+1), if `a` represents a
 *   bivariate polynomial A. 
 */ 
void to_kronecker(zz_pX& out, const Vec<zz_pX>& a, long stepX);

/** Recovers an array of coefficients `out` by applying the inverse of the
 *  function mapping A to A(X, X^`stepX`+1) to `a`.
 */
void from_kronecker(Vec<zz_pX>& out, const zz_pX& a, long stepX);

/** Sets `out` to be `a`(X, X^`stepX`+1) */ 
void to_kronecker(zz_pX& out, const zz_pXY& a, long stepX);

/** Recovers `out` by applying the inverse of the
 *  function mapping A to A(X, X^`stepX`+1) to `a`.
 */
void from_kronecker(zz_pXY& out, const zz_pX& a, long stepX);

/*------------------------------------------------------------*/
/* kronecker multiplication                                   */
/*------------------------------------------------------------*/
/** Sets `c` = `a` * `b`, computed using Kronecker substitution */
void mul_kronecker(zz_pXY& c, const zz_pXY& a, const zz_pXY& b);

/** Returns `a` * `b`, computed using Kronecker substitution */
inline zz_pXY mul_kronecker(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    mul_kronecker(c, a, b);
    return c;
}

/*------------------------------------------------------------*/
/* multiplication = kronecker for the time being              */
/*------------------------------------------------------------*/
/** Sets `c` = `a` * `b`. For the time being, we only use
 *  using Kronecker substitution 
 */
inline void mul(zz_pXY& c, const zz_pXY& a, const zz_pXY& b)
{
    mul_kronecker(c, a, b);
}

/** Returns `a` * `b`. For the time being, we only use
 *  using Kronecker substitution 
 */
inline zz_pXY mul(const zz_pXY& a, const zz_pXY& b)
{
    zz_pXY c;
    mul(c, a, b);
    return c;
}

/** Returns `a` * `b`. For the time being, we only use
 *  using Kronecker substitution 
 */
inline zz_pXY operator*(const zz_pXY& a, const zz_pXY& b)
{
    return mul(a, b);
}


/** @}  */       // doxygen multiplication

/*------------------------------------------------------------*/
/* equality test                                              */
/*------------------------------------------------------------*/
/** @name Equality test
 *  @{
 */ 

/**  Returns non-zero if `a` == `b`, 0 otherwise  */ 
inline long operator==(const zz_pXY& a, const zz_pXY& b)
{
   return a.rep == b.rep;
}

/**  Returns non-zero if `a` != `b`, 0 otherwise  */ 
inline long operator!=(const zz_pXY& a, const zz_pXY& b)
{ 
    return !(a == b); 
}
/** @}  */       // doxygen equality



/*------------------------------------------------------------*/
/* shift                                                      */ 
/*------------------------------------------------------------*/
/** @name Shifts 
 *  Index `i` >= 0 means multiplying by the `i`th power of the variable.
 *  Index `i` <= 0 means dividing by the -`i`th power of the variable.
 *  @{
 */ 

/** Sets `b` to be `a` shifted `i` places in X. */
void shift_x(zz_pXY& b, const zz_pXY &a, long i);

/** Returns `a` shifted `i` places in X. */
inline zz_pXY shift_x(const zz_pXY &a, long i)
{
    zz_pXY b;
    shift_x(b, a, i);
    return b;
}

/** Sets `b` to be `a` shifted `i` places in Y. */
void shift_y(zz_pXY& b, const zz_pXY &a, long i);

/** Returns `a` shifted `i` places in Y. */
inline zz_pXY shift_y(const zz_pXY &a, long i)
{
    zz_pXY b;
    shift_y(b, a, i);
    return b;
}

/** @}  */       // doxygen shift


/*------------------------------------------------------------*/
/* truncation                                                 */ 
/*------------------------------------------------------------*/
/** @name Truncation 
 *  @{
 */ 

/**  Sets `b` to be `a` mod X^`i`.   */ 
void trunc_x(zz_pXY& b, const zz_pXY &a, long i);

/**  Returns `a` mod X^`i`.   */ 
inline zz_pXY trunc_x(const zz_pXY &a, long i)
{
    zz_pXY b;
    trunc_x(b, a, i);
    return b;
}

/**  Sets `b` to be `a` mod Y^`i`.   */ 
void trunc_y(zz_pXY& b, const zz_pXY &a, long i);

/**  Returns `a` mod Y^`i`.   */ 
inline zz_pXY trunc_y(const zz_pXY &a, long i)
{
    zz_pXY b;
    trunc_y(b, a, i);
    return b;
}

/** @}  */       // doxygen truncation


/*------------------------------------------------------------*/
/* setting coefficients                                       */
/*------------------------------------------------------------*/
/** @name Accessing coefficients 
 *  @{
 */ 

/** Sets `a` to be the `i`th coefficient of `c`   */
void SetCoeff(zz_pXY& a, long i, const zz_pX& c);

/** @}  */       // doxygen coefficients


/*------------------------------------------------------------*/
/* resultant with respect to Y                                */
/*------------------------------------------------------------*/
/** @name Resultant
 *  Computes the resultant of `f` and `g` in Y. 
 *  Returns 0 if the computation failed.
 *  @{
 */

/** Sets `res` = resultant of `f` and `g` in Y, computed by evaluation / interpolation.
 *  Returns 0 if the computation failed.
 */
long resultant(zz_pX& res, const zz_pXY& f, const zz_pXY& g);

/** Sets `res` = resultant of `f` and `g` in Y, computed by Villard's 2018 algorithm.
 *  Returns 0 if the computation failed.
 */
long resultant_villard(zz_pX& res, const zz_pXY& f, const zz_pXY& g);

/** @}  */       // doxygen resultant


#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
