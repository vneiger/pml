#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"

#define TIME(function)                   \
    {                                    \
        t = 0.0;                         \
        nb_iter = 0;                     \
        Mat<zz_pX> a, b;                 \
        random(a, sz, sz, deg);          \
        random(b, sz, sz, deg);          \
        while (t<0.5)                    \
        {                                \
            tt = GetWallTime();          \
            Mat<zz_pX> c;                \
            function(c, a, b);           \
            t += GetWallTime()-tt;       \
            ++nb_iter;                   \
        }                                \
        t /= nb_iter;                    \
        std::cout << t/tmul << "\t";     \
    }

#define COMPARE(fn,fnt)                  \
    {                                    \
        tcomp = 0.0;                     \
        nb_iter = 0;                     \
        Mat<zz_pX> a, b;                 \
        random(a, sz, sz, deg);          \
        random(b, sz, sz, deg);          \
        while (tcomp<0.1)                \
        {                                \
            tt = GetWallTime();          \
            Mat<zz_pX> c;                \
            fn(c, a, b);                 \
            tcomp += GetWallTime()-tt;   \
            ++nb_iter;                   \
        }                                \
        tcomp /= nb_iter;                \
    }                                    \
    {                                    \
        t = 0.0;                         \
        nb_iter = 0;                     \
        Mat<zz_pX> a, b;                 \
        random(a, sz, sz, deg);          \
        random(b, sz, sz, deg);          \
        while (t<0.1)                    \
        {                                \
            tt = GetWallTime();          \
            Mat<zz_pX> c;                \
            fnt(c, a, b);                \
            t += GetWallTime()-tt;       \
            ++nb_iter;                   \
        }                                \
        t /= nb_iter;                    \
    }                                    \
    std::cout << t/tcomp << "\t";


NTL_CLIENT

void one_bench_fft(long sz, long deg)
{
    //double t, tt, tcomp;
    double t, tt;
    long nb_iter;

    cout << sz << "\t" << deg << "\t";

    //{ // warmup
    //    t=0.0;
    //    Mat<zz_pX> a, b;
    //    random(a, sz, sz, deg);
    //    random(b, sz, sz, deg);
    //    while (t<0.5)
    //    {
    //        tt = GetWallTime();
    //        Mat<zz_pX> c;
    //        multiply(c, a, b);
    //        t += GetWallTime()-tt;
    //    }
    //}

    double tmul = 0.0;
    {
        nb_iter = 0;
        Mat<zz_pX> a, b;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);
        while (tmul<0.5)
        {
            tt = GetWallTime();
            Mat<zz_pX> c;
            multiply(c, a, b);
            tmul += GetWallTime()-tt;
            ++nb_iter;
        }
        tmul /= nb_iter;
    }

    TIME(multiply)

    TIME(multiply_evaluate_FFT_matmul1)

    TIME(multiply_evaluate_FFT_matmul2)

    TIME(multiply_evaluate_FFT_matmul3)

    if (sz < 80)
        TIME(multiply_evaluate_FFT_direct_ll_type)
    else
        std::cout << "inf" << "\t";

    if (sz < 80)
        TIME(multiply_evaluate_FFT_direct)
    else
        std::cout << "inf" << "\t";

    if (deg<60)
        TIME(multiply_evaluate_dense)
    else
        std::cout << "inf" << "\t";

    if (deg<100)
        TIME(multiply_evaluate_dense2)
    else
        std::cout << "inf" << "\t";

    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    //std::vector<long> szs = { 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, };
    std::vector<long> szs = { 2, 4, 10, 20, 32, 75, 180 };
    //std::vector<long> degs =
    //{
    //    8, 10, 12, 14,
    //    16, 20, 24, 28,
    //    32, 40, 48, 56,
    //    64, 80, 96, 112,
    //    128, 160, 192, 224,
    //    256, 320, 384, 448,
    //    512, 640, 768, 896,
    //    1024, 1280, 1536, 1792,
    //    2048, 2560, 3072, 3584,
    //    4096, 5120, 6144, 7168,
    //    8192, 10240, 12288, 14336,
    //    16384, 20480, 24576, 28672,
    //    32768, 40960, 49152, 57344
    //};
    std::vector<long> degs =
    {
        8,
        32,
        128,
        512,
        2048,
        8192,
    };

    std::cout << "Bench polynomial matrix multiplication (FFT prime)" << std::endl;
    std::cout << "(ratios versus multiply)" << std::endl;
    if (nbits < 25)
    {
        zz_p::UserFFTInit(786433); // 20 bits
        cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 20 << ")" << endl;
    }
    else if (nbits < 35)
    {
        zz_p::UserFFTInit(2013265921); // 31 bits
        cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 31 << ")" << endl;
    }
    else if (nbits < 45)
    {
        zz_p::UserFFTInit(2748779069441); // 42 bits
        cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 42 << ")" << endl;
    }
    else if (nbits < 65)
    {
        zz_p::UserFFTInit(1139410705724735489); // 60 bits
        cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 60 << ")" << endl;
    }
    std::cout << "size\tdegree\tmult\tmm1\tmm2\tmm3\tdir_ll\tdir\tvdmd\tvdmd2" << std::endl;
    for (long sz : szs)
    {
        long maxdeg = 5000;
        if (sz > 128)
            maxdeg = 1600;
        if (sz > 256)
            maxdeg = 800;
        for (long d : degs)
            if (d < maxdeg)
            //if (sz*sz*sz*d < 2000000000)
                one_bench_fft(sz,d);
    }
    cout << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);
    //std::cout << std::unitbuf; // enable automatic flushing

    std::cout << "Usage: ./time_FFTx OR ./time_FFTx nbits OR ./time_FFTx sz deg" << std::endl;

    if (argc==1)
    {
        {
            std::cout << std::endl << "====NBITS: 20====" << std::endl;
            SetNumThreads(1);
            run_bench(20);
        }
        {
            std::cout << std::endl << "====NBITS: 31====" << std::endl;
            SetNumThreads(1);
            run_bench(31);
        }
        {
            std::cout << std::endl << "====NBITS: 42====" << std::endl;
            SetNumThreads(1);
            run_bench(42);
        }
        {
            std::cout << std::endl << "====NBITS: 60====" << std::endl;
            SetNumThreads(1);
            run_bench(60);
        }
    }
    if (argc==2)
    {
        SetNumThreads(1);
        warmup();
        run_bench(atoi(argv[1]));
    }
    if (argc==3)
    {
        SetNumThreads(1);
        //zz_p::FFTInit(0); // 60 bits, FFT
        //std::cout << "Bench polynomial matrix multiplication (FFT prime, 60 bits)" << std::endl;
        zz_p::UserFFTInit(786433); // FFT, 20 bits
        std::cout << "Bench polynomial matrix multiplication (FFT prime, 20 bits)" << std::endl;
        std::cout << "(ratios versus multiply)" << std::endl;
        std::cout << "size\tdegree\tmult\tmm1\tmm2\tmm3\tdir_ll\tdir\tvdmd\tvdmd2" << std::endl;
        warmup();
        one_bench_fft(atoi(argv[1]),atoi(argv[2]));
    }
    else
        throw std::invalid_argument("Usage: ./time_FFTx OR ./time_FFTx nbits OR ./time_FFTx sz deg");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
