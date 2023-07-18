#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a shift by evaluation                               */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 500; i++)
    {
        zz_pX a, b, d;
        zz_p c, e, u, v;
        std::unique_ptr<zz_pX_shift> sPtr;

        a = random_zz_pX(i);
        c = random_zz_p();
        e = random_zz_p();
        v = eval(a, e+c);

        zz_pX_shift_DAC sDac(i-1, c);
        sDac.shift(b, a);
        u = eval(b, e);
        if (u != v)
            LogicError("Error in DAC shift");

        if (p == 0 || p > i)
        {
            zz_pX_shift_large_characteristic sLarge(i-1, c);
            sLarge.shift(b, a);
            u = eval(b, e);
            if (u != v)
                LogicError("Error in large shift");
        }

        sPtr = get_shift(i-1, c);
        sPtr->shift(b, a);
        u = eval(b, e);
        if (u != v)
            LogicError("Error in ptr shift");

        b = shift(a, c);
        u = eval(b, e);
        if (u != v)
            LogicError("Error in procedural shift");

    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
    check(23068673);
    check(2);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
