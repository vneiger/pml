#ifndef __LZZ_PX_EXTRA__H
#define __LZZ_PX_EXTRA__H

#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* a class that does Taylor shift                             */
/* TODO: check characteristic                                 */
/*------------------------------------------------------------*/
class zz_pX_shift 
{
private:
    Vec<zz_p> fact, ifact;
    zz_pX v;
    long d;

public:  
    /*------------------------------------------------------------*/
    /* constructor inits a few arrays                             */
    /* d is an upper bound on the degrees of the inputs           */
    /*------------------------------------------------------------*/
    zz_pX_shift(const zz_p& c, long d);

    /*------------------------------------------------------------*/
    /* g = f(x+c)                                                 */
    /* output can alias input                                     */
    /*------------------------------------------------------------*/
    void shift(zz_pX& g, const zz_pX& f);
};

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/* creates and discards a shift object                        */
/*------------------------------------------------------------*/
void shift(zz_pX& g, const zz_pX& f, const zz_p& c);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
