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

    // kernel direct via approx
    t = 0.0; nb_iter = 0;
    while (t < 0.2)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, deg);
        tt = GetWallTime();
        VecLong shift(rdim); // uniform shift
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
        tt = GetWallTime();
        VecLong shift(rdim); // uniform shift
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
        tt = GetWallTime();
        VecLong shift(rdim); // uniform shift
        Mat<zz_pX> kerbas;
        VecLong pivind;
        kernel_basis_zls_via_approximation(kerbas,pivind,pmat,shift);
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
        tt = GetWallTime();
        Mat<zz_pX> kerbas;
        VecLong shift(rdim);
        VecLong pivind;
        kernel_basis_zls_via_interpolation(kerbas,pivind,pmat,shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter;

    std::cout << std::endl;
}

/*------------------------------------------------------------*/
/* runs bench on variety of parameters                        */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime)
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

    VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256};

    cout << "rdim\tcdim\tdeg\tdirect-app\tdirect-int\tzls-app\t\tzls-int" << endl;
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

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3)
        throw std::invalid_argument("Usage: ./time_kernel nbits fftprime");

    const long nbits = {atoi(argv[1])};
    const bool fftprime = {(atoi(argv[2])==1) ? true : false};

    warmup();

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
