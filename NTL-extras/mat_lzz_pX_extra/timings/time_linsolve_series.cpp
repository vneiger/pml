#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>

#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_linsolve.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void one_bench(long sz, long degmat, long degvec, long prec)
{
    std::cout << sz << "\t" << degmat << "\t" << degvec << "\t" << prec << "\t";
    Mat<zz_pX> A, b, u, res;

    // for timings
    double t, tt;
    long nb_iter;

    t = 0.0; nb_iter = 0;
    while (t<0.2)
    {
        random(A, sz, sz, degmat);
        random(b, sz, 1, degvec);

        tt=GetWallTime();
        solve_series_low_precision(u, A, b, prec);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    t /= nb_iter;
    cout << t << "\t";

    t=0.0; nb_iter = 0;
    while (t<0.2)
    {
        random(A, sz, sz, degmat);
        random(b, sz, 1, degvec);

        tt=GetWallTime();
        solve_series_high_precision(u, A, b, prec);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    t /= nb_iter;
    cout << t << "\t";

    if (degvec<degmat)
    {
        t=0.0; nb_iter = 0;
        while (t<0.2)
        {
            random(A, sz, sz, degmat);
            random(b, sz, 1, degvec);

            tt=GetWallTime();
            solve_series_high_order_lifting(u, A, b, prec);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        t /= nb_iter;
        cout << t << "\t";
    }
    else
        std::cout << "(degvec>=degmat)";
    cout << endl;
}

void run_bench(long nthreads, long nbits, bool fftprime, long dim=-1, long degmat=-1, long degvec=-1, long prec=-1)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench linsolve_series, FFT prime p = ";
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
        cout << "Bench linsolve_series, random prime p = ";
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

    std::cout << "sz\tdegmat\tdegvec\tprec\tlow-prec\thigh-prec\thigh-ord-lift" << std::endl;

    if (dim==-1) // then degmat==-1 && degvec==-1, default case
    {
        for (long sz : szs)
            for (long deg : degs)
                if (sz * deg < 15000)
                {
                    one_bench(sz, deg, deg-1, deg+1);
                    if (sz>5)
                        one_bench(sz, deg, deg-1, 4*deg+1);
                    one_bench(sz, deg, deg-1, sz*deg+1);
                    one_bench(sz, deg, sz*deg, sz*deg+1);
                }
    }
    else
    {
        one_bench(dim, degmat, degvec, prec);
    }
}


/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 && argc!=7)
        throw std::invalid_argument("Usage: ./time_linsolve_series nbits fftprime (dim degmat degvec prec)");
    // assume dim degmat degvec prec all positive

    const long nbits = atoi(argv[1]);
    const bool fftprime = (atoi(argv[2])==1);

    if (argc==7)
    {
        const long dim = atoi(argv[3]);
        const long degmat = atoi(argv[4]);
        const long degvec = atoi(argv[5]);
        const long prec = atoi(argv[6]);
        run_bench(1,nbits,fftprime,dim,degmat,degvec,prec);
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
