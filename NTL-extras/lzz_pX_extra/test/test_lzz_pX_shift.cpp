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
        zz_p c = random_zz_p();
        zz_pX_shift s(c, i);
        a = random_zz_pX(i);
        s.shift(b, a);
        shift(d, a, c);

        if (b != d)
        {
            cerr << "Shift inconsistent for i=" << i << endl;
            exit(-1);
        }

        zz_p e = random_zz_p();
        zz_p u, v;
        u = eval(b, e);
        v = eval(a, e+c);
        if (u != v)
        {
            cerr << "Shift error for i=" << i << endl;
            exit(-1);
        }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    long p0 = 0;
    long p1 = 23068673;
    long p2 = 288230376151711813;

    check(p0);
    check(p1);
    check(p2);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
