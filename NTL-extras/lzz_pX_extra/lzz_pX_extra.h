#ifndef __LZZ_PX_EXTRA__H
#define __LZZ_PX_EXTRA__H

#include <memory>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* an abstract class that does Taylor shift                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_shift
{
public:
    virtual ~zz_pX_shift(){}

    /*------------------------------------------------------------*/
    /* g = f(x+c)                                                 */
    /* output can alias input                                     */
    /*------------------------------------------------------------*/
    virtual void shift(zz_pX& g, const zz_pX& f) const = 0;

protected:
    long d;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class that does Taylor shift                             */
/* no assumption on the base ring, DAC algorithm              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_shift_DAC : public zz_pX_shift
{
public:  
    /*------------------------------------------------------------*/
    /* constructor inits a few arrays                             */
    /* d is a (nonstrict) upper bound on the input degree         */
    /*------------------------------------------------------------*/
    zz_pX_shift_DAC(long d, const zz_p& c);
    zz_pX_shift_DAC()
    {
        d = -1;
    }

    /*------------------------------------------------------------*/
    /* g = f(x+c)                                                 */
    /* output can alias input                                     */
    /*------------------------------------------------------------*/
    void shift(zz_pX& g, const zz_pX& f) const;

private:
    Vec<zz_pX> precomp;
    zz_p c, cc, c3;
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
