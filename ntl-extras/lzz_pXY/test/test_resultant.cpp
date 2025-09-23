#include <NTL/lzz_pX.h>
#include <iomanip>
#include <NTL/vec_lzz_p.h>

#include "util.h"
#include "lzz_pXY.h"

PML_CLIENT

zz_pXY random_zz_pXY_tdeg(long d) // tdeg < d
{
    zz_pXY f{};
    f.rep.SetLength(d);

    for (long i = 0; i < d; i++)
        f.rep[i] = random_zz_pX(d - i);
    return f;
}

/*------------------------------------------------------------*/
/* creates random bivariate polynomials                       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 150; i < 250; i += 5)
    {
        zz_pX res1, res2, res3, res4;
        double t1, t2, t3, t4;

        {
            zz_pXY f = random_zz_pXY(i, i);
            zz_pXY g = random_zz_pXY(i, i);

            t1 = GetWallTime();
            resultant(res1, f, g);
            t1 = GetWallTime()-t1;

            t2 = GetWallTime();
            resultant_villard(res2, f, g);
            t2 = GetWallTime()-t2;
        }

        {
            zz_pXY f = random_zz_pXY_tdeg(i);
            zz_pXY g = random_zz_pXY_tdeg(i);

            t3 = GetWallTime();
            resultant(res3, f, g);
            t3 = GetWallTime()-t3;

            t4 = GetWallTime();
            resultant_villard_tdeg(res4, f, g);
            t4 = GetWallTime()-t4;
        }

        // normalize results to ease comparison
        MakeMonic(res1);
        MakeMonic(res2);
        MakeMonic(res3);
        MakeMonic(res4);

        cout << "deg = " << i << endl;
        cout << "maxdeg: time naive | Villard: " << t1 << " | " << t2 << endl;
        cout << "tdeg: time naive | Villard: " << t3 << " | " << t4 << endl;
        cout << "correct maxdeg | tdeg --> ";
        cout << ((res1 == res2) ? "yes" : "no") << " | ";
        cout << ((res3 == res4) ? "yes" : "no") << endl;
        std::cout << deg(res3) << "\t" << deg(res4) << std::endl;
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    std::cout<<std::fixed;
    std::cout<<std::setprecision(5);
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
