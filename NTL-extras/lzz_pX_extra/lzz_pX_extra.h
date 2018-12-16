#ifndef __LZZ_PX_EXTRA__H
#define __LZZ_PX_EXTRA__H

/** Additional functionalities for polynomials (lzz_pX).
 *
 * \file lzz_pX_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-15
 *
 */

#include <memory> // for unique_ptr
#include <NTL/lzz_pX.h>

NTL_CLIENT

/** Returns `true` if the polynomial `a` is monic, and `false` otherwise */
inline bool is_monic(const zz_pX & a)
{ return (a.rep.length() > 0 && IsOne(a.rep[deg(a)])); }


/** @name Power series division
 * \anchor PowerSeriesDivision
 *
 *  Given two polynomials `a` and `b` and an integer `d`, with `a` invertible
 *  modulo `x^d` (that is, invertible constant coefficient), these functions
 *  compute the truncated power series division `f = b * a^{-1} mod x^d`.
 */
//@{

/** Computes `f = b * a^{-1} mod x^d`, requiring that the constant coefficient
 * of `a` is invertible (see @ref PowerSeriesDivision). Uses threshold to
 * choose between the naive algorithm and Newton iteration. The OUT parameter
 * `f` may alias the IN parameters `a` or `b`. */
void InvTruncMul(zz_pX & f, const zz_pX & b, const zz_pX & a, long d);

/** Computes and returns `b * a^{-1} mod x^d`, requiring that the constant
 * coefficient of `a` is invertible (see @ref PowerSeriesDivision). Uses
 * threshold to choose between the naive algorithm and Newton iteration. */
inline zz_pX InvTruncMul(const zz_pX & b, const zz_pX & a, long d)
{ zz_pX f; InvTruncMul(f, b, a, d); return f; }

//@} // doxygen group: Power series division

/** An abstract class for Taylor shift */
class zz_pX_shift
{
public:
    /** Empty constructor sets the bound to `-1` */
    zz_pX_shift() : d(-1) { };
    /** Constructor for a given bound */
    zz_pX_shift(long d) : d(d) { };

    /** Destructor */
    virtual ~zz_pX_shift(){}

    /** Computes the Taylor shift `g = f(x+c)` (`c` will be an attribute in
     * derived classes); `g` may alias `f` */
    virtual void shift(zz_pX & g, const zz_pX & f) const = 0;

protected:
    long d; /**< non-strict upper bound on the degree of the polynomials `f` to be shifted */
};


/** A class for Taylor shift with a divide and conquer algorithm, making no
 * assumption on the base ring `zz_p` */
class zz_pX_shift_DAC : public zz_pX_shift
{
public:  
    /** Empty constructor, sets the degree bound to `-1` */
    zz_pX_shift_DAC() { };

    /** Main constructor. Initializes the attributes `d` and `c` with the given
     * parameters, and initializes the other attributes `cc`, `c3`, `precomp`
     * as indicated in their description. */
    zz_pX_shift_DAC(long d, const zz_p& c);

    /** Computes the Taylor shift `g = f(x+c)`. The OUT parameter `g` may alias
     * the IN parameter `f` */
    void shift(zz_pX& g, const zz_pX& f) const;

private:
    Vec<zz_pX> precomp; /**< contains `(x+c)^{2**k}` for all integers `k` such
                          that `2**k` is at least `4` and strictly less than
                          `(d+1)/2` */
    zz_p c; /**< the shifting parameter */
    zz_p cc; /**< `cc = c*c`, the square of `c` */
    zz_p c3; /**< `c3 = 3*c` */
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class that does Taylor shift                             */
/* requires 1,...,d units mod p                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_shift_large_characteristic : public zz_pX_shift
{
public:  
    /*------------------------------------------------------------*/
    /* constructor inits a few arrays                             */
    /* d is a (nonstrict) upper bound on the input degree         */
    /*------------------------------------------------------------*/
    zz_pX_shift_large_characteristic(long d, const zz_p& c);
    zz_pX_shift_large_characteristic()
    {
        d = -1;
    }

    /*------------------------------------------------------------*/
    /* g = f(x+c)                                                 */
    /* output can alias input                                     */
    /*------------------------------------------------------------*/
    void shift(zz_pX& g, const zz_pX& f) const;

private:
    Vec<zz_p> fact, ifact;
    zz_pX v;
};

/*------------------------------------------------------------*/
/* returns a zz_pX_shift of the right type                    */
/*------------------------------------------------------------*/
std::unique_ptr<zz_pX_shift> get_shift(long d, const zz_p& c);

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/* creates and discards a shift object                        */
/*------------------------------------------------------------*/
void shift(zz_pX& g, const zz_pX& f, const zz_p& c);

inline zz_pX shift(const zz_pX& f, const zz_p& c)
{
    zz_pX g;
    shift(g, f, c);
    return g;
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
