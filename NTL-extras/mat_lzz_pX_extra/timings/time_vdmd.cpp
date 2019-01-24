#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

void one_bench_vdmd(long m, long n, long p, long deg)
{
    double t, tt;

    long nb;

    std::cout << m << "\t" << n << "\t" << p << "\t" << deg << "\t";

    // warmup with multiply
    t = 0.0;
    nb = 0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b, c;

        random(a, m, n, deg);
        random(b, n, p, deg);

        tt = GetWallTime();
        multiply(c, a, b);
        t += GetWallTime()-tt;
        ++nb;
    }
    std::cout << t/nb << "\t";

    // dense1
    t = 0.0;
    nb = 0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b, c;

        random(a, m, n, deg);
        random(b, n, p, deg);

        tt = GetWallTime();
        multiply_evaluate_dense(c, a, b);
        t += GetWallTime()-tt;
        ++nb;
    }
    std::cout << t/nb << "\t";

    // dense2
    t = 0.0;
    nb = 0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b, c;

        random(a, m, n, deg);
        random(b, n, p, deg);

        tt = GetWallTime();
        multiply_evaluate_dense2(c, a, b);
        t += GetWallTime()-tt;
        ++nb;
    }
    std::cout << t/nb << "\t";

    // via Mat<zz_p> mult
    t = 0.0;
    nb = 0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b, c;

        random(a, m, n, deg);
        random(b, n, p, deg);

        tt = GetWallTime();
        multiply_evaluate_FFT_matmul3(c, a, b);
        t += GetWallTime()-tt;
        ++nb;
    }
    std::cout << t/nb;


    cout << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    if (argc != 5)
    {
        std::cout << "Usage: ./time_vdmd m n p deg" << std::endl;
    }

    SetNumThreads(1);
    zz_p::FFTInit(0);
    std::cout << "m\tn\tp\tdeg\tmult\tdense1\tdense2\tmatmul3" << std::endl;
    one_bench_vdmd(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
