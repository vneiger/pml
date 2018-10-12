// NTLX
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

// LinBox
#include <linbox/linbox-config.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <linbox/ring/modular.h>
#include <givaro/zring.h>
#include <recint/rint.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>

// for LinBox, create random vector (used below to create random polynomial matrices)
template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v)
{
    size_t s = v.size();				   
    for (size_t i = 0; i < s; ++i)
        r.random(v[i]); 
}

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_bench_fft(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    double t;

    random(a, sz, sz, deg);
    random(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    if (sz<=4)
    {
        t = GetWallTime();
        multiply_waksman(c1, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (naive)" << endl;

        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_FFT(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (evaluate)" << endl;
    }
    else
    {
        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_FFT(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << endl;
    }
}

/*------------------------------------------------------------*/
void one_bench_3primes(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    double t;

    random(a, sz, sz, deg);
    random(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    if (sz<=4)
    {
        t = GetWallTime();
        multiply_waksman(c1, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (naive)" << endl;

        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_3_primes(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (3 primes)" << endl;
    }
    else
    {
        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_3_primes(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << endl;
    }
}

/*------------------------------------------------------------*/
void one_bench_geometric(long sz, long deg)
{
    Mat<zz_pX> a, b, c0, c1, c2, c4, c5;
    double t;

    random(a, sz, sz, deg);
    random(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    if (sz<=4)
    {
        t = GetWallTime();
        multiply_waksman(c1, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (naive)" << endl;

        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_geometric(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (geometric)" << endl;
    }
    else
    {
        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_geometric(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << endl;
    }
}

/*------------------------------------------------------------*/
template <typename Field>
void one_bench_multiply(long sz, long deg, Field field)
{
    double t_ntlx,t_linbox;

    {
        Mat<zz_pX> a, b, c;

        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t_ntlx = GetWallTime();
        multiply(c, a, b);
        t_ntlx = GetWallTime()-t_ntlx;
    }

    {
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field> a(field,sz,sz,deg),b(field,sz,sz,deg),c(field,sz,sz,deg+deg-1);
        typename Field::RandIter Gen(field,time(NULL));
        LinBox::PolynomialMatrixDomain<Field> PMD(field);

        // Generate random matrix of polynomials
        for (long i=0;i<sz*sz;++i)
            randomVect(Gen,a(i));
        for (long i=0;i<sz*sz;++i)
            randomVect(Gen,b(i));

        t_linbox = GetWallTime();
        PMD.mul(c,a,b);
        t_linbox = GetWallTime()-t_linbox;
    }


    cout << sz << "," << deg << "," << t_ntlx << "," << t_linbox << "," << t_ntlx/t_linbox << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long test, long nbits)
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

    if (test==0 || test==1)
    {
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
        for (size_t i=0;i<szs.size();i++)
            one_bench_fft(szs[i],degs[i]);
        cout << endl;
    }

    if (test==0 || test==2)
    {
        long prime = NTL::GenPrime_long(nbits);
        zz_p::init(prime);
        std::cout << "Bench polynomial matrix multiplication (3 primes)" << std::endl;
        cout << "p = " << zz_p::modulus() << "  (normal prime, bit length = " << nbits << ")" << endl;
        for (size_t i=0;i<szs.size();i++)
            one_bench_3primes(szs[i],degs[i]);
        cout << endl;
    }

    if (test==0 || test==3)
    {
        long prime = NTL::GenPrime_long(nbits);
        zz_p::init(prime);
        std::cout << "Bench polynomial matrix multiplication (geometric)" << std::endl;
        cout << "p = " << zz_p::modulus() << "  (normal prime, bit length = " << nbits << ")" << endl;
        for (size_t i=0;i<szs.size();i++)
            one_bench_geometric(szs[i],degs[i]);
        cout << endl;
    }

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
        20,30,40,50,60 // normal primes with 20, 30, 40, 50, 60 bits
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
        std::cout << "size,degree,ntlx,linbox,ratio(ntlx/linbox)" << std::endl;

        // NTLx initialize field
        zz_p::UserFFTInit(primes[p]);
        if (p==0) // FFT prime with 20 bits
        {
            Givaro::Modular<double> field((int32_t)primes[p]);
            for (size_t i=0; i<szs.size(); ++i)
                one_bench_multiply(szs[i],degs[i],field);
        }
        else // FFT prime with > 28 bits
        {
            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(primes[p]);
            for (size_t i=0; i<szs.size(); ++i)
                one_bench_multiply(szs[i],degs[i],field);
        }
        cout << endl << endl;
    }

    for (size_t p = 4; p < primes.size(); ++p)
    {
        long prime = NTL::GenPrime_long(primes[p]);
        std::cout << "Bench polynomial matrix multiplication (normal prime)" << std::endl;
        cout << "p = " << prime << "  (normal prime, bit length = " << primes[p] << ")" << endl;
        std::cout << "size,degree,ntlx,linbox,ratio(ntlx/linbox)" << std::endl;
        // NTLx initialize field
        zz_p::init(prime);
        if (p==0) // normal prime with <29 bits
        {
            // LinBox initialize field
            Givaro::Modular<double> field((int32_t)prime);
            for (size_t i=0;i<szs.size();i++)
                one_bench_multiply(szs[i],degs[i],field);
        }
        else // FFT prime with >= 29 bits
        {
            // LinBox initialize field
            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(prime);
            for (size_t i=0;i<szs.size();i++)
                one_bench_multiply(szs[i],degs[i],field);
        }
        cout << endl << endl;
    }
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(4);
    //SetNumThreads(4);

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc==1)
    {
        warmup();
        run_bench();
    }
    if (argc==2 || argc==3)
    {
        long test = 0; // default: run all benchs
        if (argc==3)
            test = atoi(argv[2]);
        long nbits = atoi(argv[1]);
        warmup();
        run_bench(test,nbits);
    }
    if (argc>3)
        throw std::invalid_argument("Usage: ./time_multiply_comparelinbox OR ./time_multiply_comparelinbox nbits OR ./time_multiply_comparelinbox nbits test");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
