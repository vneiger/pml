#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>

#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void check(long p)
{
    SetNumThreads(1);

    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    std::vector<long> szs =
        {
            1, 10, 20, 30, 40, 50, 100
        };

    std::vector<long> degs =
        {
            100, 150, 200, 250, 300, 350, 400
        };

    const double thresh = 0.01;
    for (size_t si = 0; si < szs.size(); si++)
    {
        long sz = szs[si];
        for (size_t di = 0; di < degs.size(); di++)
        {
            double t;
            long nb, deg, prec;
            Mat<zz_pX> A, b, u, res;

            deg = degs[di];
            prec = 4 * deg;
            random(A, sz, sz, deg);
            random(b, sz, 1, deg);

            cout << sz << " " << deg << " ";

            t = get_time();
            nb = 0;
            do
            {
                solve_series_low_precision(u, A, b, prec);
                nb++;
            }
            while ((get_time()-t) <= thresh);
            t = (get_time()-t) / nb;
            cout << t << " ";

            t = get_time();
            nb = 0;
            do
            {
                solve_series_high_precision(u, A, b, prec);
                nb++;
            }
            while ((get_time()-t) <= thresh);
            t = (get_time()-t) / nb;
            cout << t << " ";

            cout << endl;
        }
    }
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
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
