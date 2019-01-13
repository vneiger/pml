#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_extra.h"

NTL_CLIENT

void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 5; i < 1000; i++)
    {
        zz_pX a, b, c;

        random(a, i);
        random(b, i);
        add_LeftShift(c, a, b, 5);
        if (c != (a + (b<<5)))
            LogicError("Error in add_LeftShift1");

        random(a, i);
        random(b, 2*i);
        add_LeftShift(c, a, b, 12);
        if (c != (a + (b<<12)))
            LogicError("Error in add_LeftShift2");

        random(a, 2*i);
        random(b, i);
        add_LeftShift(c, a, b, 24);
        if (c != (a + (b<<24)))
            LogicError("Error in add_LeftShift3");

        random(c, 2*i);
        random(b, i);
        a = c; // save c
        add_LeftShift(c, c, b, 24);
        if (c != (a + (b<<24)))
            LogicError("Error in add_LeftShift4");

        random(c, i);
        random(b, 2*i);
        a = c; // save c
        add_LeftShift(c, c, b, 71);
        if (c != (a + (b<<71)))
            LogicError("Error in add_LeftShift5");
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
