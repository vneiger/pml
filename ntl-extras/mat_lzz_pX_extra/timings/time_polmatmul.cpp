/** \file time_polmatmul.cpp
 *
 * This file benchmarks the PML polynomial matrix multiplication, and several
 * variants of it to check the thresholds are relatively finely chosen.
 *
 */
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

#define TIME(function)                   \
    {                                    \
        t = 0.0;                         \
        nb_iter = 0;                     \
        Mat<zz_pX> a, b;                 \
        random(a, rdim, idim, deg);      \
        random(b, idim, cdim, deg);      \
        while (t<0.5)                    \
        {                                \
            tt = GetWallTime();          \
            Mat<zz_pX> c;                \
            function(c, a, b);           \
            t += GetWallTime()-tt;       \
            ++nb_iter;                   \
        }                                \
        t /= nb_iter;                    \
        std::cout << t << "\t";          \
    }

NTL_CLIENT

/** Bench of polmatmul for given dimension and degree, modulo FFT prime.
 *
 *  This function takes the parameters (matrix dimensions, degree bound)
 *  creates two random matrices a and b, and measures the time it takes to
 *  compute c = a*b. The prime field must already have been initialized
 *  with an FFT prime.
 *
 * \param rdim row dimension
 * \param idim inner dimension
 * \param cdim column dimension
 * \param deg degree bound
 * \return void
 */
void benchmark_polmatmul_fftprime(long rdim, long idim, long cdim, long deg)
{
    double t, tt;
    long nb_iter;

    cout << rdim << "\t" << idim << "\t" << cdim << "\t" << deg << "\t";

    { // warmup
        t=0.0;
        Mat<zz_pX> a, b;
        random(a, rdim, idim, deg);
        random(b, idim, cdim, deg);
        while (t<0.5)
        {
            tt = GetWallTime();
            Mat<zz_pX> c;
            multiply(c, a, b);
            t += GetWallTime()-tt;
        }
    }

    TIME(multiply)

    TIME(multiply_evaluate_FFT_matmul1)

    TIME(multiply_evaluate_FFT_matmul2)

    TIME(multiply_evaluate_FFT_matmul3)

    if (rdim*idim*cdim < 80*80*80)
        TIME(multiply_evaluate_FFT_direct_ll_type)
    else
        std::cout << "inf" << "\t";

    if (rdim*idim*cdim < 80*80*80)
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

/** Bench of polmatmul for given dimension and degree.
 *
 *  This function takes the parameters (matrix dimensions, degree bound)
 *  creates two random matrices a and b, and measures the time it takes to
 *  compute c = a*b. The prime field must already have been initialized.
 *
 * \param rdim row dimension
 * \param idim inner dimension
 * \param cdim column dimension
 * \param deg degree bound
 * \return void
 */
void benchmark_polmatmul(long rdim, long idim, long cdim, long deg)
{
    double t, tt;
    long nb_iter;

    cout << rdim << "\t" << idim << "\t" << cdim << "\t" << deg << "\t";

    { // warmup
        t=0.0;
        Mat<zz_pX> a, b;
        random(a, rdim, idim, deg);
        random(b, idim, cdim, deg);
        while (t<0.5)
        {
            tt = GetWallTime();
            Mat<zz_pX> c;
            multiply(c, a, b);
            t += GetWallTime()-tt;
        }
    }

    TIME(multiply)

    TIME(multiply_3_primes);

    TIME(multiply_evaluate_geometric);

    if (deg<100)
        TIME(multiply_evaluate_dense)
    else
        std::cout << "inf" << "\t";

    if (deg<150)
        TIME(multiply_evaluate_dense2)
    else
        std::cout << "inf" << "\t";

    cout << endl;
}


/** Launches a series of benchmarks for a given prime size.
 *
 *  Launches benchmark for many matrix dimensions and degrees for a given
 *  bitlength for the prime defining the field of coefficients. If fftprime is
 *  true, this picks an FFT prime not too far from the given bitlength.
 *
 * \param nbits bitlength of the prime modulus
 * \param fftprime boolean, whether to take an FFT prime
 * \return void
 */
void benchmark_nbits(long nbits, bool fftprime)
{
    std::vector<long> dims = { 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, };
    //std::vector<long> dims = { 2, 4, 10, 20, 32, 75, 180 };
    std::vector<long> degs =
    {
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
    //std::vector<long> degs =
    //{
    //    8,
    //    32,
    //    128,
    //    512,
    //    2048,
    //    8192,
    //};

    if (fftprime)
    {
        std::cout << "Bench square polynomial matrix multiplication (FFT prime):" << std::endl;
        if (nbits < 25)
            zz_p::UserFFTInit(786433); // 20 bits
        else if (nbits < 35)
            zz_p::UserFFTInit(2013265921); // 31 bits
        else if (nbits < 45)
            zz_p::UserFFTInit(2748779069441); // 42 bits
        else if (nbits < 65)
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
        std::cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << NumBits(zz_p::modulus()) << ")" << std::endl;
        std::cout << "dim\tdim\tdim\tdegree\tmult\tmm1\tmm2\tmm3\tdir_ll\tdir\tvdmd\tvdmd2" << std::endl;
        for (long dim : dims)
            for (long d : degs)
                if (dim*dim*dim*d < 30000000000)
                    benchmark_polmatmul_fftprime(dim,dim,dim,d);
        std::cout << std::endl;
    }
    else
    {
        std::cout << "Bench square polynomial matrix multiplication:" << std::endl;
        zz_p::init(GenPrime_long(nbits));
        std::cout << "p = " << zz_p::modulus() << "  (prime, bit length = " << NumBits(zz_p::modulus()) << ")" << std::endl;
        std::cout << "dim\tdim\tdim\tdegree\tmult\t3pri\tev-geo\tvdmd\tvdmd2" << std::endl;
        for (long dim : dims)
            for (long d : degs)
                if (dim*dim*dim*d < 30000000000)
                    benchmark_polmatmul(dim,dim,dim,d);
        std::cout << std::endl;
    }
}

/** Launches single benchmark for specified parameters.
 *
 *  Launches one benchmark for the given square matrix dimensions, matrix
 *  degree, and bitlength for the prime defining the field of coefficients. If
 *  fftprime is true, this picks an FFT prime not too far from the given
 *  bitlength.
 *
 * \param nbits bitlength of the prime modulus
 * \param dim square matrix dimension
 * \param deg matrix degree
 * \param fftprime boolean, whether to take an FFT prime
 * \return void
 */
void benchmark_nbits_dim_deg(long nbits, long dim, long deg, bool fftprime)
{
    if (fftprime)
    {
        std::cout << "Bench square polynomial matrix multiplication (FFT prime):" << std::endl;
        if (nbits < 25)
            zz_p::UserFFTInit(786433); // 20 bits
        else if (nbits < 35)
            zz_p::UserFFTInit(2013265921); // 31 bits
        else if (nbits < 45)
            zz_p::UserFFTInit(2748779069441); // 42 bits
        else if (nbits < 65)
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
        std::cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << NumBits(zz_p::modulus()) << ")" << std::endl;
        std::cout << "dim\tdim\tdim\tdegree\tmult\tmm1\tmm2\tmm3\tdir_ll\tdir\tvdmd\tvdmd2" << std::endl;
        benchmark_polmatmul_fftprime(dim,dim,dim,deg);
    }
    else
    {
        std::cout << "Bench square polynomial matrix multiplication:" << std::endl;
        zz_p::init(GenPrime_long(nbits));
        std::cout << "p = " << zz_p::modulus() << "  (prime, bit length = " << NumBits(zz_p::modulus()) << ")" << std::endl;
        std::cout << "dim\tdim\tdim\tdegree\tmult\t3pri\tev-geo\tvdmd\tvdmd2" << std::endl;
        benchmark_polmatmul(dim,dim,dim,deg);
    }
}

int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);
    //std::cout << std::unitbuf; // enable automatic flushing

    if (argc!=3 && argc!=5)
    {
        std::cout << "Usage: " << argv[0] << " nbits fftprime OR " << argv[0] << " nbits dim deg fftprime" << std::endl;
        return 1;
    }

    if (argc==3)
    {
        long nbits = atoi(argv[1]);
        bool fftprime = (atoi(argv[2]) == 1);
        benchmark_nbits(nbits,fftprime);
        return 0;
    }

    if (argc==5)
    {
        long nbits = atoi(argv[1]);
        long dim = atoi(argv[2]);
        long deg = atoi(argv[3]);
        bool fftprime = (atoi(argv[4]) == 1);
        benchmark_nbits_dim_deg(nbits,dim,deg,fftprime);
    }

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
