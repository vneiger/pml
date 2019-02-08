#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

void one_bench_fft(long sz, long deg)
{
    double t, tt;
    long nb_iter;

    cout << sz << "\t" << deg << "\t";

    { // warmup
        t=0.0;
        Mat<zz_pX> a, b, c;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);
        while (t<0.1)
        {
            tt = GetWallTime();
            multiply(c, a, b);
            t += GetWallTime()-tt;
        }
    }

    { // multiply
        t = 0.0;
        nb_iter = 0;
        Mat<zz_pX> a, b, c;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);
        while (t<0.1)
        {
            tt = GetWallTime();
            multiply(c, a, b);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        std::cout << t/nb_iter << "\t";
    }

    t = 0.0;
    nb_iter = 0;
    while (t<0.1)
    {
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        tt = GetWallTime();
        multiply_evaluate_FFT_matmul1(c, a, b);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    t = 0.0;
    nb_iter = 0;
    while (t<0.1)
    {
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        tt = GetWallTime();
        multiply_evaluate_FFT_matmul2(c, a, b);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    t = 0.0;
    nb_iter = 0;
    while (t<0.1)
    {
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        tt = GetWallTime();
        multiply_evaluate_FFT_matmul3(c, a, b);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    t = 0.0;
    nb_iter = 0;
    while (t<0.1)
    {
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        tt = GetWallTime();
        multiply_evaluate_FFT_direct_ll_type(c, a, b);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    { // direct (without ll type)
        t = 0.0;
        nb_iter = 0;
        Mat<zz_pX> a, b, c;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);
        while (t<0.1)
        {
            tt = GetWallTime();
            multiply_evaluate_FFT_direct(c, a, b);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        std::cout << t/nb_iter << "\t";
    }

    //std::cout << "==================old===================" << std::endl;
    if (deg<600)
    {
        t = 0.0;
        nb_iter = 0;
        while (t<0.1)
        {
            Mat<zz_pX> a, b, c;

            random(a, sz, sz, deg);
            random(b, sz, sz, deg);

            tt = GetWallTime();
            multiply_evaluate_dense(c, a, b);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        std::cout << t/nb_iter << "\t";
    }
    else
    {
        std::cout << "inf" << "\t";
    }

    if (deg<350)
    {
        t = 0.0;
        nb_iter = 0;
        while (t<0.1)
        {
            Mat<zz_pX> a, b, c;

            random(a, sz, sz, deg);
            random(b, sz, sz, deg);

            tt = GetWallTime();
            multiply_evaluate_dense2(c, a, b);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        std::cout << t/nb_iter << "\t";
    }
    else
    {
        std::cout << "inf" << "\t";
    }

    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    std::vector<long> szs =
    {
        2,2,2,2,2,2,2,2,2,2,2,2,2,2,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        8,8,8,8,8,8,8,8,8,8,8,8,8,
        16,16,16,16,16,16,16,16,16,16,16,16,16,
        32,32,32,32,32,32,32,32,32,32,32,32,
        64,64,64,64,64,64,64,64,64,64,
        128,128,128,128,128,128,128,128,
        256,256,256,256,256,256,
        512,512,512,512,
        1024,1024,1024,
    };
    std::vector<long> degs =
    {
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,131070,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,131070,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,
        15,30,60,120,240,480,960,1920,3840,7680,
        15,30,60,120,240,480,960,1920,
        15,30,60,120,240,480,
        15,30,60,120,
        15,30,60,
    };

    std::cout << "Bench polynomial matrix multiplication (FFT prime)" << std::endl;
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
    //std::cout << "size\tdegree\tmult.\tmm1\tmm1new\tmm2\tmm2new\tmm3\tmm3new\tdirect\tdirect2\tvdmd\tvdmd2\twinner" << std::endl;
    for (size_t i=0;i<szs.size();i++)
        one_bench_fft(szs[i],degs[i]);
    cout << endl;
}

void run_bench()
{
    std::vector<long> szs =
    {
        2,2,2,2,2,2,2,2,2,2,2,2,2,2,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        8,8,8,8,8,8,8,8,8,8,8,8,8,
        16,16,16,16,16,16,16,16,16,16,16,16,16,
        32,32,32,32,32,32,32,32,32,32,32,32,
        64,64,64,64,64,64,64,64,64,64,
        128,128,128,128,128,128,128,128,
        256,256,256,256,256,256,
        512,512,512,512,
        1024,1024,1024,
    };
    std::vector<long> degs =
    {
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,131070,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,131070,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,61440,
        15,30,60,120,240,480,960,1920,3840,7680,15360,30720,
        15,30,60,120,240,480,960,1920,3840,7680,
        15,30,60,120,240,480,960,1920,
        15,30,60,120,240,480,
        15,30,60,120,
        15,30,60,
    };
    std::vector<long> primes =
    {
        786433,2013265921,2748779069441,1139410705724735489,   // FFT primes with 20,31,42,60 bits
    };

    std::cout << "Bench polynomial matrix multiplication (FFT prime)" << std::endl;
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
        //std::cout << "size\tdegree\tmult.\tmm1\tmm1new\tmm2\tmm2new\tmm3\tmm3new\tdirect\tdirect2\tvdmd\tvdmd2\twinner" << std::endl;

        // NTLx initialize field
        zz_p::UserFFTInit(primes[p]);
        for (size_t i=0; i<szs.size(); ++i)
            one_bench_fft(szs[i],degs[i]);
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
        //std::cout << "size\tdegree\tmult.\tmm1\tmm1new\tmm2\tmm2new\tmm3\tmm3new\tdirect\tdirect2\tvdmd\tvdmd2\twinner" << std::endl;
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
