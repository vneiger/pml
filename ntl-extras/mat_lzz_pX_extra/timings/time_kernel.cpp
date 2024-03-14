// NTLX
#include <iomanip>

#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_kernel.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* runs one bench                                             */
/*------------------------------------------------------------*/
void one_bench_kernel(long rdim, long cdim, long deg)
{
    double t,tt;
    long nb_iter;

    std::cout << rdim << "\t" << cdim << "\t" << deg << "\t";

    bool reductionshift = false;
    bool hermiteshift = false;

    // kernel direct via approx
    t = 0.0; nb_iter = 0;
    while (t < 0.2)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, deg);
        VecLong shift(rdim); // uniform shift
        if (reductionshift)
            for (long i = rdim/2; i < rdim; ++i) // basis reduction-like shift
                shift[i] = cdim*deg;
        if (hermiteshift)
            for (long i = 1; i < rdim; ++i) // Hermite shift
                shift[i] = shift[i-1] + cdim*deg;
        tt = GetWallTime();
        Mat<zz_pX> kerbas;
        VecLong pivind;
        kernel_basis_via_approximation(kerbas,pivind,pmat,shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    // kernel direct via interpolation
    t = 0.0; nb_iter = 0;
    while (t < 0.2)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, deg);
        VecLong shift(rdim); // uniform shift
        if (reductionshift)
            for (long i = rdim/2; i < rdim; ++i) // basis reduction-like shift
                shift[i] = cdim*deg;
        if (hermiteshift)
            for (long i = 1; i < rdim; ++i) // Hermite shift
                shift[i] = shift[i-1] + cdim*deg;
        tt = GetWallTime();
        Mat<zz_pX> kerbas;
        VecLong pivind;
        kernel_basis_via_interpolation(kerbas,pivind,pmat,shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    // kernel ZLS via approx
    t = 0.0; nb_iter = 0;
    while (t < 0.2)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, deg);
        VecLong shift(rdim); // uniform shift
        if (reductionshift)
            for (long i = rdim/2; i < rdim; ++i) // basis reduction-like shift
                shift[i] = cdim*deg;
        if (hermiteshift)
            for (long i = 1; i < rdim; ++i) // Hermite shift
                shift[i] = shift[i-1] + cdim*deg;
        tt = GetWallTime();
        Mat<zz_pX> kerbas;
        kernel_basis_zls_via_approximation(kerbas,pmat,shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    // kernel ZLS via interp
    t = 0.0; nb_iter = 0;
    while (t<0.2)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, deg);
        VecLong shift(rdim); // uniform shift
        if (reductionshift)
            for (long i = rdim/2; i < rdim; ++i) // basis reduction-like shift
                shift[i] = cdim*deg;
        if (hermiteshift)
            for (long i = 1; i < rdim; ++i) // Hermite shift
                shift[i] = shift[i-1] + cdim*deg;
        tt = GetWallTime();
        Mat<zz_pX> kerbas;
        kernel_basis_zls_via_interpolation(kerbas,pmat,shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter;

    std::cout << std::endl;
}

/*------------------------------------------------------------*/
/* runs bench on variety of parameters                        */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime, long rdim=-1, long cdim=-1, long degree=-1)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench kernel, FFT prime p = ";
        if (nbits <= 20)
        {
            zz_p::UserFFTInit(786433); // 20 bits
            cout << zz_p::modulus() << ", bit length = " << 20 << endl;
        }
        else if (nbits <= 31)
        {
            zz_p::UserFFTInit(2013265921); // 31 bits
            cout << zz_p::modulus() << ", bit length = " << 31 << endl;
        }
        else if (nbits <= 42)
        {
            zz_p::UserFFTInit(2748779069441); // 42 bits
            cout << zz_p::modulus() << ", bit length = " << 42 << endl;
        }
        else if (nbits <= 60)
        {
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
            cout << zz_p::modulus() << ", bit length = " << 60 << endl;
        }
        else
        {
            std::cout << "FFT prime with more than 60 bits is not supported" << std::endl;
            return;
        }
    }
    else
    {
        cout << "Bench kernel, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    cout << "rdim\tcdim\tdeg\tdirect-app\tdirect-int\tzls-app\t\tzls-int" << endl;

    if (rdim==-1) // means cdim==-1 and degree==-1 as well
    {
        VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256};

        for (size_t i=0; i<szs.size(); ++i)
        {
            long interval = ceil( (double)szs[i] / 4);
            //for (long j=1; 2*j<3*szs[i]; j+=interval)
            for (long j=1; j<szs[i]; j+=interval) // only rdim<cdim for the moment
            {
                long max_order=4096;
                for (long k=2; k<=max_order; k=2*k)
                    one_bench_kernel(szs[i],j,k);
            }
        }
        cout << endl;
    }
    else
        one_bench_kernel(rdim,cdim,degree);
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 && argc!=6)
        throw std::invalid_argument("Usage: ./time_kernel nbits fftprime (rdim cdim degree)");
    // assume rdim>0 , cdim>0, degree>0

    const long nbits = atoi(argv[1]);
    const bool fftprime = (atoi(argv[2])==1);

    if (argc==6)
    {
        const long rdim = atoi(argv[3]);
        const long cdim = atoi(argv[4]);
        const long degree = atoi(argv[5]);
        run_bench(1,nbits,fftprime,rdim,cdim,degree);
    }
    else
        run_bench(1,nbits,fftprime);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
