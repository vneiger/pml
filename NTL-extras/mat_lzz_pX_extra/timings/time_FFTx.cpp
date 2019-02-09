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
        while (t<0.1)                    \
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
    double t, tt, tcomp;
    long nb_iter;

    cout << sz << "\t" << deg << "\t";

    { // warmup
        t=0.0;
        Mat<zz_pX> a, b;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);
        while (t<0.2)
        {
            tt = GetWallTime();
            Mat<zz_pX> c;
            multiply(c, a, b);
            t += GetWallTime()-tt;
        }
    }

    double tmul = 0.0;
    {
        nb_iter = 0;
        Mat<zz_pX> a, b;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);
        while (tmul<0.1)
        {
            tt = GetWallTime();
            Mat<zz_pX> c;
            multiply(c, a, b);
            tmul += GetWallTime()-tt;
            ++nb_iter;
        }
        tmul /= nb_iter;
    }

    COMPARE(multiply_evaluate_FFT_matmul1,
            multiply_evaluate_FFT_matmul1_trunc)
    COMPARE(multiply_evaluate_FFT_matmul2,
            multiply_evaluate_FFT_matmul2_trunc)
    COMPARE(multiply_evaluate_FFT_matmul3,
            multiply_evaluate_FFT_matmul3_trunc)

    //TIME(multiply_evaluate_FFT_matmul1)

    //TIME(multiply_evaluate_FFT_matmul2)

    //TIME(multiply_evaluate_FFT_matmul3)

    //if (sz < 80)
    //    TIME(multiply_evaluate_FFT_direct_ll_type)
    //else
    //    std::cout << "inf" << "\t";

    //if (sz < 80)
    //    TIME(multiply_evaluate_FFT_direct)
    //else
    //    std::cout << "inf" << "\t";

    //if (deg<400)
    //    TIME(multiply_evaluate_dense)
    //else
    //    std::cout << "inf" << "\t";

    //if (deg<500)
    //    TIME(multiply_evaluate_dense2)
    //else
    //    std::cout << "inf" << "\t";

    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    //std::vector<long> szs = { 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, };
    std::vector<long> szs = { 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, };
    std::vector<long> degs =
    {
        2, 
        4, 6,
        8, 10, 12, 14,
        16, 20, 24, 28,
        32, 40, 48, 56,
        64, 80, 96, 112,
        128, 160, 192, 224,
        256, 320, 384, 448,
        512, 640, 768, 896,
        1024, 1280, 1536, 1792,
        2048, 2560, 3072, 3584,
        4096, 5120, 6144, 7168,
        8192, 10240, 12288, 14336,
        16384, 20480, 24576, 28672,
        32768, 40960, 49152, 57344
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
    //std::cout << "size\tdegree\tmm1\tmm2\tmm3\tdir_ll\tdirect\tvdmd\tvdmd2" << std::endl;
    std::cout << "size\tdegree\tmm1\tmm1-tr\tmm2\tmm2-tr\tmm3\tmm3-tr\tdir_ll\tdirect\tvdmd\tvdmd2" << std::endl;
    for (size_t sz : szs)
        for (size_t d : degs)
            if (sz*sz*d < 100000000)
                one_bench_fft(sz,d);
    cout << endl;
}

void run_bench()
{
    std::vector<long> szs = { 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, };
    std::vector<long> degs =
    {
        2, 
        4, 6,
        8, 10, 12, 14,
        16, 20, 24, 28,
        32, 40, 48, 56,
        64, 80, 96, 112,
        128, 160, 192, 224,
        256, 320, 384, 448,
        512, 640, 768, 896,
        1024, 1280, 1536, 1792,
        2048, 2560, 3072, 3584,
        4096, 5120, 6144, 7168,
        8192, 10240, 12288, 14336,
        16384, 20480, 24576, 28672,
        32768, 40960, 49152, 57344
    };

    std::vector<long> primes =
    {
        786433,2013265921,2748779069441,1139410705724735489,   // FFT primes with 20,31,42,60 bits
    };

    std::cout << "Bench polynomial matrix multiplication (FFT prime)" << std::endl;
    std::cout << "(ratios versus multiply)" << std::endl;
    for (long p = 0; p < 4; ++p)
    {
        cout << "p = " << primes[p] << "  (FFT prime, bit length = ";
        switch (p)
        {
        case 0: cout << 20 << ")" << endl; break;
        case 1: cout << 31 << ")" << endl; break;
        case 2: cout << 42 << ")" << endl; break;
        case 3: cout << 60 << ")" << endl; break;
        }
        //std::cout << "size\tdegree\tmm1\tmm2\tmm3\tdir_ll\tdir\tvdmd\tvdmd2" << std::endl;
        std::cout << "size\tdegree\tmm1\tmm1-tr\tmm2\tmm2-tr\tmm3\tmm3-tr\tdir_ll\tdirect\tvdmd\tvdmd2" << std::endl;

        zz_p::UserFFTInit(primes[p]);
        for (size_t sz : szs)
            for (size_t d : degs)
                if (sz*sz*d < 100000000)
                    one_bench_fft(sz,d);
        cout << endl << endl;
    }
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
        zz_p::FFTInit(0); // 60 bits, FFT
        std::cout << "Bench polynomial matrix multiplication (FFT prime, 60 bits)" << std::endl;
        //zz_p::UserFFTInit(786433); // FFT, 20 bits
        //std::cout << "Bench polynomial matrix multiplication (FFT prime, 20 bits)" << std::endl;
        std::cout << "(ratios versus multiply)" << std::endl;
        std::cout << "size\tdegree\tmm1\tmm1-tr\tmm2\tmm2-tr\tmm3\tmm3-tr\tdir_ll\tdirect\tvdmd\tvdmd2" << std::endl;
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
