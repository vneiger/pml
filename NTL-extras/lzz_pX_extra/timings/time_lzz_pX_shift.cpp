#include <iomanip>

#include "util.h"
#include "lzz_pX_extra.h"

NTL_CLIENT

/*--------------*/
/* time a shift */
/*--------------*/
void time_one(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    cout << "prime " << p << endl;
    std::cout << "i\tDACprecomp\tDACshift\tLargePrecomp\tLargeShift" << std::endl;

    for (long i = 1; i < 100; i++)
    {
        zz_pX a, b;
        zz_p c;
        zz_pX_shift_DAC sDac;
        zz_pX_shift_large_characteristic sLarge;

        a = random_zz_pX(i);
        c = random_zz_p();

        cout << i << "\t";
        double t;
        long nb;
        const double thresh = 0.1;

        nb = 0;
        t = get_time();
        do
        {
            sDac = zz_pX_shift_DAC(i-1, c);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << "\t";


        nb = 0;
        t = get_time();
        do
        {
            sDac.shift(b, a);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << "\t";


        nb = 0;
        t = get_time();
        do
        {
            sLarge = zz_pX_shift_large_characteristic(i-1, c);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << "\t";


        nb = 0;
        t = get_time();
        do
        {
            sLarge.shift(b, a);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << "\t";
        
        cout << endl;
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup(); 
    time_one(0);
    time_one(23068673);
    time_one(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
