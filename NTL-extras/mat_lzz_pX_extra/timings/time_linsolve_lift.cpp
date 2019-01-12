#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>
#include <iomanip>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_linsolve.h"

NTL_CLIENT

void one_bench(long sz, long deg)
{
    Mat<zz_pX> A;
    Vec<zz_pX> b, u;
    zz_pX den;

    // to record timings
    double t, tt;
    long nb_iter;

    // to record the best option
    double tmin = 1000000.0;
    long idx=-1;

    for (long nb = 1; nb < 6; ++nb)
    {
        nb_iter = 0; t=0.0;
        while (t<0.2)
        {
            random(A, sz, sz, deg);
            random(b, sz, deg);
            tt = get_time();
            linsolve_via_series(u, den, A, b, nb);
            t += get_time()-tt;
            ++nb_iter;
        }
        t /= nb_iter;
        cout << t << "\t";
        if (t < tmin)
        {
            tmin = t;
            idx = nb;
        }
    }

    cout << idx << endl;
}

void run_bench(long nthreads, long nbits, bool fftprime)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench linsolve, FFT prime p = ";
        if (nbits < 25)
        {
            zz_p::UserFFTInit(786433); // 20 bits
            cout << zz_p::modulus() << ", bit length = " << 20 << endl;
        }
        else if (nbits < 35)
        {
            zz_p::UserFFTInit(2013265921); // 31 bits
            cout << zz_p::modulus() << ", bit length = " << 31 << endl;
        }
        else if (nbits < 45)
        {
            zz_p::UserFFTInit(2748779069441); // 42 bits
            cout << zz_p::modulus() << ", bit length = " << 42 << endl;
        }
        else if (nbits < 61)
        {
            zz_p::FFTInit(0); // 60 bits
            cout << zz_p::modulus() << ", bit length = " << NumBits(zz_p::modulus()) << endl;
        }
        else
        {
            std::cout << "Asking for FFT prime with too large bitsize (> 60). Exiting." << std::endl;
            return;
        }
    }
    else
    {
        cout << "Bench mbasis, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    VecLong szs =
    {
        5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100
    };

    VecLong degs =
    {
        15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400,
        500, 750, 1000
    };

    std::cout << "sz\tdeg\t" << std::endl;

    for (long sz : szs)
        for (long deg : degs)
            if (sz * deg < 15000)
            {
                std::cout << sz << "\t" << deg << "\t";
                one_bench(sz, deg);
            }

}

int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3)
        throw std::invalid_argument("Usage: ./time_linsolve nbits fftprime");

    const long nbits = atoi(argv[1]);
    const bool fftprime = (atoi(argv[2])==1);

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
