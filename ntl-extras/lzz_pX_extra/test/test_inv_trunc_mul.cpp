#include <NTL/lzz_pX.h>
#include <assert.h>

#include "util.h"
#include "lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a power series division                             */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 200; i < 1000; i++)
    {
        zz_pX a, b, c;
        a = random_zz_pX(i);
        b = random_zz_pX(i);

        InvTruncMul(c, b, a, 2*i);
        assert (IsZero(trunc(c*a - b, 2*i)));
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(288230376151711813);
    check(0);
    check(786433);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
