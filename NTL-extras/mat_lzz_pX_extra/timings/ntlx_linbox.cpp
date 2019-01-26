// NTLX
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>
#include <iomanip>

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
#include <linbox/algorithms/polynomial-matrix/order-basis.h>


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
/* compare polynomial matrix multiplication                   */
/*------------------------------------------------------------*/
template <typename Field>
void one_bench_multiply(long sz, long deg, Field field)
{
    double t_ntlx,t_linbox;
    long nb_iter;
    double t;

    nb_iter=0; t_ntlx=0.0;
    while (t_ntlx<0.2)
    {
        Mat<zz_pX> a, b;
        random(a, sz, sz, deg);
        random(b, sz, sz, deg);

        t = GetWallTime();
        Mat<zz_pX> c;
        multiply(c, a, b);
        t_ntlx += GetWallTime()-t;

        ++nb_iter;
    }
    t_ntlx /= nb_iter;

    nb_iter=0; t_linbox=0.0;
    while (t_linbox<0.2)
    {
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field> a(field,sz,sz,deg),b(field,sz,sz,deg);
        typename Field::RandIter Gen(field,time(NULL));
        LinBox::PolynomialMatrixDomain<Field> PMD(field);
        // Generate random matrix of polynomials
        for (long i=0;i<sz*sz;++i)
            randomVect(Gen,a(i));
        for (long i=0;i<sz*sz;++i)
            randomVect(Gen,b(i));

        t = GetWallTime();
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field> c(field,sz,sz,2*deg-1);
        PMD.mul(c,a,b);
        t_linbox += GetWallTime() - t;

        ++nb_iter;
    }
    t_linbox /= nb_iter;

    cout << sz << "\t" << deg << "\t" << t_ntlx << "\t" << t_linbox << "\t" << t_ntlx/t_linbox << endl;
}

/*------------------------------------------------------------*/
/* compare divide and conquer approximant basis               */
/* (2n x n matrix, uniform order, uniform shift)              */
/*------------------------------------------------------------*/
template <typename Field>
void one_bench_pmbasis(long n, long order, Field field)
{
    //double t_ntlx,t_linbox,t_ntlx_generic;
    double t_ntlx,t_linbox;
    long nb_iter;
    double t;

    nb_iter=0; t_ntlx=0.0;
    VecLong shift_ntl(2*n,0);
    while (t_ntlx<0.2)
    {
        Mat<zz_pX> mat;
        random(mat, 2*n, n, order);

        t = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift_ntl);
        pmbasis(appbas, mat, order, rdeg);
        t_ntlx += GetWallTime()-t;

        ++nb_iter;
    }
    t_ntlx /= nb_iter;

    nb_iter=0; t_linbox=0.0;
    std::vector<unsigned long> shift_linbox(2*n,0);
    while (t_linbox<0.2)
    {
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field>
        mat(field,2*n,n,order);
        typename Field::RandIter Gen(field,time(NULL));
        LinBox::OrderBasis<Field> OB(field);
        // Generate random matrix of polynomials
        for (long i=0;i<2*n*n;++i)
            randomVect(Gen,mat(i));

        t = GetWallTime();
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field>
        appbas(field,2*n,2*n,order);
		OB.PM_Basis(appbas, mat, order, shift_linbox);
        t_linbox += GetWallTime() - t;

        ++nb_iter;
    }
    t_linbox /= nb_iter;

    //nb_iter=0; t_ntlx_generic=0.0;
    //while (t_ntlx_generic<0.2)
    //{
    //    Mat<zz_pX> mat;
    //    random(mat, 2*n, n, order);

    //    t = GetWallTime();
    //    Mat<zz_pX> appbas;
    //    pmbasis_generic_2n_n(appbas, mat, order);
    //    t_ntlx_generic += GetWallTime()-t;

    //    ++nb_iter;
    //}
    //t_ntlx_generic /= nb_iter;

    //cout << 2*n << "\t" << n << "\t" << order << "\t" << t_ntlx << "\t" << t_linbox << "\t" << t_ntlx_generic << "\t" << t_ntlx/t_linbox << "\t" << t_ntlx_generic / t_ntlx << endl;
    cout << 2*n << "\t" << n << "\t" << order << "\t" << t_ntlx << "\t" << t_linbox << "\t" << t_ntlx/t_linbox << endl;
}


/*------------------------------------------------------------*/
/* compare iterative approximant basis                        */
/* (2n x n matrix, uniform order, uniform shift)              */
/*------------------------------------------------------------*/
template <typename Field>
void one_bench_mbasis(long n, long order, Field field)
{
    double t_ntlx,t_linbox;
    long nb_iter;
    double t;

    nb_iter=0; t_ntlx=0.0;
    VecLong shift_ntl(2*n,0);
    while (t_ntlx<0.2)
    {
        Mat<zz_pX> mat;
        random(mat, 2*n, n, order);

        t = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift_ntl);
        mbasis(appbas, mat, order, rdeg);
        t_ntlx += GetWallTime()-t;

        ++nb_iter;
    }
    t_ntlx /= nb_iter;

    nb_iter=0; t_linbox=0.0;
    std::vector<unsigned long> shift_linbox(2*n,0);
    while (t_linbox<0.2)
    {
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field>
        mat(field,2*n,n,order);
        typename Field::RandIter Gen(field,time(NULL));
        LinBox::OrderBasis<Field> OB(field);
        // Generate random matrix of polynomials
        for (long i=0;i<2*n*n;++i)
            randomVect(Gen,mat(i));

        t = GetWallTime();
        LinBox::PolynomialMatrix<LinBox::PMType::polfirst,LinBox::PMStorage::plain,Field>
        appbas(field,2*n,2*n,order);
		OB.M_Basis(appbas, mat, order, shift_linbox);
        t_linbox += GetWallTime() - t;

        ++nb_iter;
    }
    t_linbox /= nb_iter;

    cout << 2*n << "\t" << n << "\t" << order << "\t" << t_ntlx << "\t" << t_linbox << "\t" << t_ntlx/t_linbox << endl;
}

/*------------------------------------------------------------*/
/* run bench on collection of size/degree                     */
/*------------------------------------------------------------*/
void run_bench(long test, long nbits, bool fftprime)
{
    VecLong szs =
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
    VecLong degs =
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

    if (fftprime)
    {
        if (nbits < 25)
        {
            zz_p::UserFFTInit(786433); // 20 bits
            cout << "Prime p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 20 << ")" << endl;
        }
        else if (nbits < 35)
        {
            zz_p::UserFFTInit(2013265921); // 31 bits
            cout << "Prime p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 31 << ")" << endl;
        }
        else if (nbits < 45)
        {
            zz_p::UserFFTInit(2748779069441); // 42 bits
            cout << "Prime p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 42 << ")" << endl;
        }
        else if (nbits < 65)
        {
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
            cout << "Prime p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 60 << ")" << endl;
        }
    }
    else // normal prime
    {
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << "Prime p = " << zz_p::modulus() << "  (bit length = " << nbits << ")" << endl;
    }

    if ((test==0) | (test==1) )
    {
        std::cout << "Bench polynomial matrix multiplication" << std::endl;
        std::cout << "size\tdegree\tntlx\tlinbox\tratio(ntlx/linbox)" << std::endl;

        if (NumBits(zz_p::modulus()) < 26)
        {
            // LinBox initialize field
            Givaro::Modular<double> field((int32_t)zz_p::modulus());
            for (size_t i=0;i<szs.size();i++)
                one_bench_multiply(szs[i],degs[i],field);
        }
        else // nbits >= 26
        {
            // LinBox initialize field
            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(zz_p::modulus());
            for (size_t i=0;i<szs.size();i++)
                one_bench_multiply(szs[i],degs[i],field);
        }
        cout << endl << endl;
    }

    if ((test==0) | (test==2))
    {
        std::cout << "Bench pmbasis (divide and conquer minimal approximant basis)" << std::endl;
        //std::cout << "rowdim\tcoldim\tdegree\tntlx\tlinbox\tntlxgen\tntl/lb\tntlxgen/ntlx" << std::endl;
        std::cout << "rowdim\tcoldim\tdegree\tntlx\tlinbox\tntl/lb" << std::endl;

        if (NumBits(zz_p::modulus()) < 26)
        {
            // LinBox initialize field
            Givaro::Modular<double> field((int32_t)zz_p::modulus());
            for (size_t i=0;i<szs.size();i++)
                one_bench_pmbasis(szs[i],degs[i],field);
        }
        else // nbits >= 26
        {
            // LinBox initialize field
            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(zz_p::modulus());
            for (size_t i=0;i<szs.size();i++)
                one_bench_pmbasis(szs[i],degs[i],field);
        }
        cout << endl << endl;
    }

    if ((test==0) | (test==3))
    {
        std::cout << "Bench mbasis (iterative minimal approximant basis)" << std::endl;
        //std::cout << "rowdim\tcoldim\tdegree\tntlx\tlinbox\tntlxgen\tntl/lb\tntlxgen/ntlx" << std::endl;
        std::cout << "rowdim\tcoldim\tdegree\tntlx\tlinbox\tntl/lb" << std::endl;

        if (NumBits(zz_p::modulus()) < 26)
        {
            // LinBox initialize field
            Givaro::Modular<double> field((int32_t)zz_p::modulus());
            for (size_t i=0;i<szs.size();i++)
                if (degs[i] <= 128)
                    one_bench_mbasis(szs[i],degs[i],field);
        }
        else // nbits >= 26
        {
            // LinBox initialize field
            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(zz_p::modulus());
            for (size_t i=0;i<szs.size();i++)
                if (degs[i] <= 128)
                    one_bench_mbasis(szs[i],degs[i],field);
        }
        cout << endl << endl;
    }
}

//void run_bench()
//{
//    VecLong szs =
//    {
//        2,2,2,2,2,2,2,2,2,2,2,2,
//        4,4,4,4,4,4,4,4,4,4,4,4,
//        8,8,8,8,8,8,8,8,8,8,8,8,
//        16,16,16,16,16,16,16,16,16,16,16,16,
//        32,32,32,32,32,32,32,32,32,32,32,
//        64,64,64,64,64,64,64,64,64,
//        128,128,128,128,128,128,128,
//        256,256,256,256,256,
//        512,512,512,
//        1024,
//    };
//    VecLong degs =
//    {
//        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
//        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
//        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
//        30,60,120,250,510,1020,2040,4090,8190,16380,32760,131070,
//        30,60,120,250,510,1020,2040,4090,8190,16380,32760,
//        30,60,120,250,510,1020,2040,4090,8190,
//        30,60,120,250,510,1020,2040,
//        30,60,120,250,510,
//        30,60,120,
//        30,
//    };
//    VecLong primes =
//    {
//        786433,2013265921,2748779069441,1139410705724735489,   // FFT primes with 20,31,42,60 bits
//        20,30,40,50,60 // normal primes with 20, 30, 40, 50, 60 bits
//    };
//
//    std::cout << "Bench polynomial matrix multiplication (FFT prime)" << std::endl;
//    for (long p = 0; p < 4; ++p)
//    {
//        cout << "p = " << primes[p] << "  (FFT prime, bit length = ";
//        switch (p)
//        {
//            case 0: cout << 20 << ")" << endl; break;
//            case 1: cout << 31 << ")" << endl; break;
//            case 2: cout << 42 << ")" << endl; break;
//            case 3: cout << 60 << ")" << endl; break;
//        }
//        std::cout << "size\tdegree\tntlx\tlinbox\tratio(ntlx/linbox)" << std::endl;
//
//        // NTLx initialize field
//        zz_p::UserFFTInit(primes[p]);
//        if (p==0) // FFT prime with 20 bits
//        {
//            Givaro::Modular<double> field((int32_t)primes[p]);
//            for (size_t i=0; i<szs.size(); ++i)
//            {
//                one_bench_multiply(szs[i],degs[i],field);
//                one_bench_multiply(szs[i],3*degs[i]/2,field);
//            }
//        }
//        else // FFT prime with > 28 bits
//        {
//            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(primes[p]);
//            for (size_t i=0; i<szs.size(); ++i)
//            {
//                one_bench_multiply(szs[i],degs[i],field);
//                one_bench_multiply(szs[i],3*degs[i]/2,field);
//            }
//        }
//        cout << endl << endl;
//    }
//
//    for (size_t p = 4; p < primes.size(); ++p)
//    {
//        long prime = NTL::GenPrime_long(primes[p]);
//        std::cout << "Bench polynomial matrix multiplication (normal prime)" << std::endl;
//        cout << "p = " << prime << "  (normal prime, bit length = " << primes[p] << ")" << endl;
//        std::cout << "size\tdegree\tntlx\tlinbox\tratio(ntlx/linbox)" << std::endl;
//        // NTLx initialize field
//        zz_p::init(prime);
//        if (p==0) // normal prime with <29 bits
//        {
//            // LinBox initialize field
//            Givaro::Modular<double> field((int32_t)prime);
//            for (size_t i=0;i<szs.size();i++)
//            {
//                one_bench_multiply(szs[i],degs[i],field);
//                one_bench_multiply(szs[i],3*degs[i]/2,field);
//            }
//        }
//        else // FFT prime with >= 29 bits
//        {
//            // LinBox initialize field
//            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(prime);
//            for (size_t i=0;i<szs.size();i++)
//            {
//                one_bench_multiply(szs[i],degs[i],field);
//                one_bench_multiply(szs[i],3*degs[i]/2,field);
//            }
//        }
//        cout << endl << endl;
//    }
//
//    std::cout << "Bench pmbasis (divide and conquer minimal approximant basis)" << std::endl;
//    for (long p = 0; p < 4; ++p)
//    {
//        cout << "p = " << primes[p] << "  (FFT prime, bit length = ";
//        switch (p)
//        {
//            case 0: cout << 20 << ")" << endl; break;
//            case 1: cout << 31 << ")" << endl; break;
//            case 2: cout << 42 << ")" << endl; break;
//            case 3: cout << 60 << ")" << endl; break;
//        }
//        std::cout << "rowdim\tcoldim\tdegree\tntlx\tlinbox\tntlxgen\tntl/lb\tntlxgen/ntlx" << std::endl;
//
//        // NTLx initialize field
//        zz_p::UserFFTInit(primes[p]);
//        if (p==0) // FFT prime with 20 bits
//        {
//            Givaro::Modular<double> field((int32_t)primes[p]);
//            for (size_t i=0; i<szs.size(); ++i)
//            {
//                one_bench_pmbasis(szs[i],degs[i],field);
//                one_bench_pmbasis(szs[i],3*degs[i]/2,field);
//            }
//        }
//        else // FFT prime with > 28 bits
//        {
//            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(primes[p]);
//            for (size_t i=0; i<szs.size(); ++i)
//            {
//                one_bench_pmbasis(szs[i],degs[i],field);
//                one_bench_pmbasis(szs[i],3*degs[i]/2,field);
//            }
//        }
//        cout << endl << endl;
//    }
//
//    for (size_t p = 4; p < primes.size(); ++p)
//    {
//        long prime = NTL::GenPrime_long(primes[p]);
//        std::cout << "Bench pmbasis (divide and conquer minimal approximant basis)" << std::endl;
//        cout << "p = " << prime << "  (normal prime, bit length = " << primes[p] << ")" << endl;
//        std::cout << "rowdim\tcoldim\tdegree\tntlx\tlinbox\tntlxgen\tntl/lb\tntlxgen/ntlx" << std::endl;
//        // NTLx initialize field
//        zz_p::init(prime);
//        if (p==0) // normal prime with <29 bits
//        {
//            // LinBox initialize field
//            Givaro::Modular<double> field((int32_t)prime);
//            for (size_t i=0;i<szs.size();i++)
//            {
//                one_bench_pmbasis(szs[i],degs[i],field);
//                one_bench_pmbasis(szs[i],3*degs[i]/2,field);
//            }
//        }
//        else // FFT prime with >= 29 bits
//        {
//            // LinBox initialize field
//            Givaro::Modular<RecInt::ruint128,RecInt::ruint256> field(prime);
//            for (size_t i=0;i<szs.size();i++)
//            {
//                one_bench_pmbasis(szs[i],degs[i],field);
//                one_bench_pmbasis(szs[i],3*degs[i]/2,field);
//            }
//        }
//        cout << endl << endl;
//    }
//}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(4);

    // NTL threads
    SetNumThreads(1);
    //SetNumThreads(4);

    // OpenMP threads
    omp_set_num_threads(1);
    //omp_set_num_threads(4);

    if (argc!=4)
    {
        std::cout << "Usage: ./ntlx_linbox test nbits fftprime" << std::endl;
        std::cout << "  test = {0: all; 1: polynomial matrix multiplication; 2: pmbasis, 3: mbasis}" << std::endl;
        std::cout << "  nbits : integer from 2 to 63, number of bits in prime defining the base field" << std::endl;
        std::cout << "  fftprime = {0,1} : whether to use an FFT prime, if possible" << std::endl;
        return 0;
    }

    long test = atoi(argv[1]);
    long nbits = atoi(argv[2]);
    bool fftprime = (atoi(argv[3])==1);

    if (test < 0 || test > 3)
    {
        std::cout << "test = {0: all; 1: polynomial matrix multiplication; 2: pmbasis, 3: mbasis}" << std::endl;
        return 0;
    }
    if (nbits < 3 || nbits > 60)
    {
        std::cout << "nbits : integer from 3 to 63, number of bits in prime defining the base field" << std::endl;
        return 0;
    }

    run_bench(test,nbits,fftprime);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
