#include <iomanip>
#include <vector>

#include "util.h"
#include "lzz_pX_extra.h"

NTL_CLIENT

void time_one(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    cout << "prime " << p << endl;
    std::cout << "d\tmain\t\tnewton\t\tplain" << std::endl;

    std::vector<long> degs = {2, 3, 5, 10, 15, 20, 30, 40, 50, 75, 100, 200, 300, 500, 1000, 5000, 10000, 20000, 50000, 100000 };

    for (long d : degs)
    {
        cout << d << "\t";

        double t = 0.0;
        double tt;
        long nb = 0;
        const double thresh = 0.1;

        // main interface
        while (t<thresh)
        {
            // random input
            zz_pX a, b;
            do
            {
                random(a, d);
            } while (coeff(a,0) == 0);  // we do our tests over a field, nonzero means invertible

            random(b, d);

            tt = GetWallTime();
            zz_pX f;
            InvTruncMul(f, b, a, d);
            t += GetWallTime() - tt;
            ++nb;
        }
        cout << t/nb << "\t";

        //// Newton
        //t = 0.0; nb = 0;
        //if (d>=NTL_zz_pX_NEWTON_CROSSOVER)
        //{
        //    while (t<thresh)
        //    {
        //        // random input
        //        zz_pX a, b;
        //        do
        //        {
        //            random(a, d);
        //        } while (coeff(a,0) == 0);  // we do our tests over a field, nonzero means invertible

        //        random(b, d);

        //        tt = GetWallTime();
        //        zz_pX f;
        //        NewtonInvTruncMul(f, b, a, d);
        //        t += GetWallTime() - tt;
        //        ++nb;
        //    }
        //    cout << t/nb << "\t";
        //}
        //else
        //    cout << "nan\t\t";

        //// plain
        //if (d<10000)
        //{
        //    t = 0.0; nb = 0;
        //    while (t<thresh)
        //    {
        //        // random input
        //        zz_pX a, b;
        //        do
        //        {
        //            random(a, d);
        //        } while (coeff(a,0) == 0);  // we do our tests over a field, nonzero means invertible

        //        random(b, d);

        //        tt = GetWallTime();
        //        zz_pX f;
        //        PlainInvTruncMul(f, b, a, d);
        //        t += GetWallTime() - tt;
        //        ++nb;
        //    }
        //    cout << t/nb << "\t";
        //}

        cout << endl;
    }
}

/*------------------------------------------------------------*/
/* main just calls time_one() for different primes            */
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
