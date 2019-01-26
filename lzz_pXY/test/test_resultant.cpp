#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "util.h"
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

    for (long i = 300; i < 1000; i += 5000)
    {
        zz_pXY f = random_zz_pXY(i, i);
        zz_pXY g = random_zz_pXY(i, i);
        zz_pX res1, res2;
        double t1, t2;

        cout << "i=" << i << endl;
        t1 = get_time();
        resultant(res1, f, g);
        t1 = get_time()-t1;

        t2 = get_time();
        resultant_villard(res2, f, g);
        t2 = get_time()-t2;

        cout << t1 << " " << t2 << endl;
        cout << deg(res1) << " " << deg(GCD(res1, res2)) << " ";
        cout << endl;
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
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
