#include <iomanip>
#include <assert.h>
#include <NTL/lzz_pX.h>

#include "util.h"
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

    cout << p << endl;

    for (long i = 1; i < 100; i++)
    {
        zz_pX a, b;
        zz_p c;
        zz_pX_shift_DAC sDac;
        zz_pX_shift_large_characteristic sLarge;

        a = random_zz_pX(i);
        c = random_zz_p();

        cout << i << " ";
        double t;
        long nb;
        const double thresh = 0.02;

        nb = 0;
        t = get_time();
        do
        {
            sDac = zz_pX_shift_DAC(i-1, c);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";


        nb = 0;
        t = get_time();
        do
        {
            sDac.shift(b, a);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";


        nb = 0;
        t = get_time();
        do
        {
            sLarge = zz_pX_shift_large_characteristic(i-1, c);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";


        nb = 0;
        t = get_time();
        do
        {
            sLarge.shift(b, a);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";
        
        cout << endl;
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup(); 
    check(0);
    check(23068673);
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
