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
    double t0=0.0, t1=0.0, t2=0.0, t2bis=0.0, t2ter=0.0, t3bis=0.0, t3ter=0.0, t3=0.0, t4=0.0, t5=0.0, t6=0.0, t7=0.0;

    long nb0=0, nb1=0, nb2bis=0, nb2=0, nb2ter=0, nb3bis=0, nb3ter=0, nb3=0, nb4=0, nb5=0, nb6=0, nb7=0;

    cout << sz << "\t" << deg << "\t";

    // warmup
    while (t0<0.1)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply(c, a, b);
        t = GetWallTime()-t;
        t0 += t;
    }

    t0 = 0.0;
    while (t0<0.1)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply(c, a, b);
        t = GetWallTime()-t;
        t0 += t;
        ++nb0;
    }
    t0 = t0/nb0;

    std::cout << t0 << "\t";

    while (t1<0.1)
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
    std::cout << t1 << "\t";

    while (t2<0.1)
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
    std::cout << t2 << "\t";

    while (t2bis<0.1)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul2bis(c, a, b);
        t = GetWallTime()-t;
        t2bis += t;
        ++nb2bis;
    }
    t2bis = t2bis/nb2bis;
    std::cout << t2bis << "\t";

    while (t2ter<0.1)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul2ter(c, a, b);
        t = GetWallTime()-t;
        t2ter += t;
        ++nb2ter;
    }
    t2ter = t2ter/nb2ter;
    std::cout << t2ter << "\t";

    while (t3<0.1)
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
    std::cout << t3 << "\t";

    while (t3bis<0.1)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul3bis(c, a, b);
        t = GetWallTime()-t;
        t3bis += t;
        ++nb3bis;
    }
    t3bis = t3bis/nb3bis;
    std::cout << t3bis << "\t";

    while (t3ter<0.1)
    {
        double t;
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        multiply_evaluate_FFT_matmul3ter(c, a, b);
        t = GetWallTime()-t;
        t3ter += t;
        ++nb3ter;
    }
    t3ter = t3ter/nb3ter;
    std::cout << t3ter << "\t";

    //while (t4<0.1)
    //{
    //    double t;
    //    Mat<zz_pX> a, b, c;

    //    random(a, sz, sz, deg);
    //    random(b, sz, sz, deg);

    //    t = GetWallTime();
    //    multiply_evaluate_FFT_direct(c, a, b);
    //    t = GetWallTime()-t;
    //    t4 += t;
    //    ++nb4;
    //}
    //t4 = t4/nb4;
    std::cout << t4 << "\t";

    //while (t6<0.1)
    //{
    //    double t;
    //    Mat<zz_pX> a, b, c;

    //    random(a, sz, sz, deg);
    //    random(b, sz, sz, deg);

    //    t = GetWallTime();
    //    multiply_evaluate_FFT_direct_no_ll(c, a, b);
    //    t = GetWallTime()-t;
    //    t6 += t;
    //    ++nb6;
    //}
    //t6 = t6/nb6;
    std::cout << t6 << "\t";

    //if (deg<500)
    //{
    //    while (t5<0.1)
    //    {
    //        double t;
    //        Mat<zz_pX> a, b, c;

    //        random(a, sz, sz, deg);
    //        random(b, sz, sz, deg);

    //        t = GetWallTime();
    //        multiply_evaluate_dense(c, a, b);
    //        t = GetWallTime()-t;
    //        t5 += t;
    //        ++nb5;
    //    }
    //    t5 = t5/nb5;
    //}
    //else
    //{
    //    nb5=1; // to avoid div by zero
    //    t5=INFINITY; // to make sure this is not the best below
    //}
    std::cout << t5 << "\t";

    //if (deg<400)
    //{
    //    while (t7<0.1)
    //    {
    //        double t;
    //        Mat<zz_pX> a, b, c;

    //        random(a, sz, sz, deg);
    //        random(b, sz, sz, deg);

    //        t = GetWallTime();
    //        multiply_evaluate_dense2(c, a, b);
    //        t = GetWallTime()-t;
    //        t7 += t;
    //        ++nb7;
    //    }
    //    t7 = t7/nb7;
    //}
    //else
    //{
    //    nb7=1; // to avoid div by zero
    //    t7=INFINITY; // to make sure this is not the best below
    //}
    std::cout << t7 << "\t";

    std::vector<double> times = {t0, t1, t2, t2bis, t2ter, t3, t3bis, t3ter, t4, t6, t5, t7};

    // timings
    //cout << sz << "\t" << deg << "\t";
    //for (size_t i = 0; i < times.size(); ++i)
    //    cout << times[i] << "\t";

    // ratios versus multiply
    auto min = std::min_element(times.begin(), times.end());
    //for (size_t i = 1; i < times.size(); ++i)
    //    cout << times[i] / times[0] << "\t";
    //cout << times[i] / *min << "\t";

    // winner
    switch (min - times.begin())
    {
        case 0: cout << "multiply"; break;
        case 1: cout << "mm1"; break;
        case 2: cout << "mm2"; break;
        case 3: cout << "mm2bis"; break;
        case 4: cout << "mm2ter"; break;
        case 5: cout << "mm3"; break;
        case 6: cout << "mm3bis"; break;
        case 7: cout << "mm3ter"; break;
        case 8: cout << "direct"; break;
        case 9: cout << "direct no LL"; break;
        case 10: cout << "vandermonde"; break;
        case 11: cout << "vandermonde2"; break;
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
    std::cout << "size\tdegree\tmult.\tmm1\tmm2\tmm2bis\tmm2ter\tmm3\tmm3bis\tmm3ter\tdirect\tdirect2\tvdmd\tvdmd2\twinner" << std::endl;
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
        std::cout << "size\tdegree\tmult.\tmm1\tmm2\tmm2bis\tmm2ter\tmm3\tmm3bis\tmm3ter\tdirect\tdirect2\tvdmd\tvdmd2\twinner" << std::endl;

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
    std::cout << std::unitbuf; // enable automatic flushing

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
        zz_p::FFTInit(0); // 60 bits, FFT
        std::cout << "Bench polynomial matrix multiplication (FFT prime, 60 bits)" << std::endl;
        //zz_p::UserFFTInit(786433); // FFT, 20 bits
        //std::cout << "Bench polynomial matrix multiplication (FFT prime, 20 bits)" << std::endl;
        std::cout << "size\tdegree\tmult.\tmm1\tmm2\tmm2bis\tmm2ter\tmm3\tmm3bis\tmm3ter\tdirect\tdirect2\tvdmd\tvdmd2\twinner" << std::endl;
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
