#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "lzz_pXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates random bivariate polynomials                       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 50; i += 1)
    {
        for (long j = 1; j < 50; j += 1)
        {
            zz_pXY f = random_zz_pXY(i, j);
        }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
    check(786433);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
