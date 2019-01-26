#include <NTL/lzz_p.h>

#include "lzz_p_extra.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a is a non-unit                                      */
/*------------------------------------------------------------*/
long order(const zz_p& a)
{
    if (a == 0)
        return -1;

    long o = 1;
    zz_p ap = a;
    while (ap != 1)
    {
        if (o == zz_p::modulus())
            return -1;

        mul(ap, ap, a);
        ++o;
    }
    return o;
}

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* a = 0 if no element was found                              */
/* does (by default) 100 trials                               */
/* by default, asks that all (a^i-1) are units, i=1..ord-1    */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord, long nb_trials, bool strong)
{
    if ((zz_p::modulus() - 1) < ord)
        LogicError("order too large with respect to field cardinality");

    long nb = 0;
    while (nb < nb_trials)
    {
        random(a);
        bool ok = true;
        zz_p prod = a;
        zz_p ap = a;

        for (long i = 1; i < ord; i++)
        {
            prod *= (1 - ap);
            if (ap == 1)
            {
                ok = false;
                break;
            }
            ap *= a;
        }

        if (ok)
        {
            if (! strong)
                return;
            if (strong && GCD(prod.LoopHole(), zz_p::modulus()) == 1)
                return;
        }
        nb++;
    }
    clear(a);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
