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
    double t1=0.0, t2=0.0, t3=0.0, t4=0.0, t5=0.0;

    long nb1=0, nb2=0, nb3=0, nb4=0, nb5=0;

    while (t1<0.2)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul1(c, a, b);
        t = GetWallTime()-t;
        t1 += t;
        ++nb1;
    }
    t1 = t1/nb1;

    while (t2<0.2)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul2(c, a, b);
        t = GetWallTime()-t;
        t2 += t;
        ++nb2;
    }
    t2 = t2/nb2;

    while (t3<0.2)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul3(c, a, b);
        t = GetWallTime()-t;
        t3 += t;
        ++nb3;
    }
    t3 = t3/nb3;

    while (t4<0.2)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_direct(c, a, b);
        t = GetWallTime()-t;
        t4 += t;
        ++nb4;
    }
    t4 = t4/nb4;

    if (deg<70)
    {
        while (t5<0.2)
        {
            double t;
            Mat<zz_pX> a, b, c;

            random(a, sz, sz, deg);
            random(b, sz, sz, deg);

            t = GetWallTime();
            multiply_evaluate_dense(c, a, b);
            t = GetWallTime()-t;
            t5 += t;
            ++nb5;
        }
        t5 = t5/nb5;
    }
    else
    {
        nb5=1; // to avoid div by zero
        t5=INFINITY; // to make sure this is not the best below
    }


    std::vector<double> times = {t1, t2, t3, t4, t5};

    // timings
    cout << sz << "\t" << deg << "\t";
    for (size_t i = 0; i < times.size(); ++i)
        cout << times[i] << "\t";

    // ratios
    auto min = std::min_element(times.begin(), times.end());
    for (size_t i = 0; i < times.size(); ++i)
        cout << times[i] / *min << "\t";

    // winner
    cout << (min - times.begin() + 1) << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    std::vector<long> szs =
    {
        2,2,2,2,2,2,2,2,2,2,2,2,
        4,4,4,4,4,4,4,4,4,4,4,4,
        8,8,8,8,8,8,8,8,8,8,8,8,
        16,16,16,16,16,16,16,16,16,16,16,16,
        32,32,32,32,32,32,32,32,32,32,32,
        64,64,64,64,64,64,64,64,64,
        128,128,128,128,128,128,128,
        256,256,256,256,256,
        512,512,512,
        1024,1024,
    };
    std::vector<long> degs =
    {
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,
        32,64,128,256,512,1024,2048,4096,8192,
        32,64,128,256,512,1024,2048,
        32,64,128,256,512,
        32,64,128,
        32,64,
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
    std::cout << "size\tdegree\tmatmul1\tmatmul2\tmatmul3\tdirect\tvdmd\tspeedups\t\t\t\twin" << std::endl;
    for (size_t i=0;i<szs.size();i++)
        one_bench_fft(szs[i],degs[i]);
    cout << endl;
}

void run_bench()
{
    std::vector<long> szs =
    {
        2,2,2,2,2,2,2,2,2,2,2,2,
        4,4,4,4,4,4,4,4,4,4,4,4,
        8,8,8,8,8,8,8,8,8,8,8,8,
        16,16,16,16,16,16,16,16,16,16,16,16,
        32,32,32,32,32,32,32,32,32,32,32,
        64,64,64,64,64,64,64,64,64,
        128,128,128,128,128,128,128,
        256,256,256,256,256,
        512,512,512,
        1024,
    };
    std::vector<long> degs =
    {
        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
        30,60,120,250,510,1020,2040,4090,8190,16380,32760,
        30,60,120,250,510,1020,2040,4090,8190,
        30,60,120,250,510,1020,2040,
        30,60,120,250,510,
        30,60,120,
        30,
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
        std::cout << "size\tdegree\tmatmul1\tmatmul2\tmatmul3\tdirect\tvdmd\tspeedups\t\t\t\twin" << std::endl;

        // NTLx initialize field
        zz_p::UserFFTInit(primes[p]);
        if (p==0) // FFT prime with 20 bits
        {
            for (size_t i=0; i<szs.size(); ++i)
                one_bench_fft(szs[i],degs[i]);
        }
        else // FFT prime with > 28 bits
        {
            for (size_t i=0; i<szs.size(); ++i)
                one_bench_fft(szs[i],degs[i]);
        }
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

    std::cout << "Usage: ./time_FFTx OR ./time_FFTx nbits OR ./time_FFTx sz deg" << std::endl;

    if (argc==1)
    {
        {
            std::cout << "====NTHREADS: 1====" << std::endl;
            SetNumThreads(1);
            warmup();
            run_bench();
        }
        {
            std::cout << "====NTHREADS: 2====" << std::endl;
            SetNumThreads(2);
            warmup();
            run_bench();
        }
        {
            std::cout << "====NTHREADS: 3====" << std::endl;
            SetNumThreads(3);
            warmup();
            run_bench();
        }
        {
            std::cout << "====NTHREADS: 4====" << std::endl;
            SetNumThreads(4);
            warmup();
            run_bench();
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
        zz_p::UserFFTInit(1139410705724735489); // 60 bits
        std::cout << "Bench polynomial matrix multiplication (FFT prime, 60 bits)" << std::endl;
        //zz_p::UserFFTInit(786433); // 20 bits
        //std::cout << "Bench polynomial matrix multiplication (FFT prime, 20 bits)" << std::endl;
        std::cout << "size\tdegree\tmatmul1\tmatmul2\tmatmul3\tdirect\tvdmd\tspeedups\t\t\t\twin" << std::endl;
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
