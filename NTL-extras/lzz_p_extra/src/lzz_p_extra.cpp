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
    {
        return -1;
    }
    long o = 1;
    zz_p ap = a;
    while (ap != 1)
    {
        if (o == zz_p::modulus())
        {
            return -1;
        }
        ap *= a;
        o++;
    }
    return o;
}

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* a = 0 if no element was found                              */
/* does (by default) 100 trials                               */
/* by default, asks that all (a^i-1) are units, i=1..ord-1    */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord, long nb_trials, long strong)
{

    long p = zz_p::modulus();
    if ((p - 1) < ord)
    {
        LogicError("order too large with respect to p");
    }

    long nb = 0;
    while (nb < nb_trials)
    {
        a = random_zz_p();
        long ok = 1;
        zz_p prod = a;
        zz_p ap = a;

        for (long i = 1; i < ord; i++)
        {
            prod *= (1 - ap);
            if (ap == 1)
            {
                ok = 0;
                break;
            }
            ap *= a;
        }

        if (ok == 1)
        {
            if (! strong)
                return;
            if (strong && GCD(prod.LoopHole(), p) == 1)
                return;
        }
        nb++;
    }
    a = 0;
}


/*------------------------------------------------------------*/
/* 1 if the current prime can be used as an FFT prime         */
/*------------------------------------------------------------*/
long is_FFT_prime()
{
    return (zz_pInfo->p_info != NULL);
}

/*------------------------------------------------------------*/
/* 1 if the current prime can be used as an FFT prime         */
/* with transforms of size up to 2^k                          */
/*------------------------------------------------------------*/
long is_FFT_ready(long k)
{
    return (is_FFT_prime() && k <= zz_pInfo->MaxRoot);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
