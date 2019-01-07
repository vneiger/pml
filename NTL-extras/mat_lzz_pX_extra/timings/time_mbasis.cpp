#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

NTL_CLIENT

/*------------------------------------------------------------*/
/* run one bench for specified rdim,cdim,order                */
/*------------------------------------------------------------*/
void one_bench_mbasis(long rdim, long cdim, long order)
{
    VecLong shift(rdim,0);

    double t1,t2;

    long nb_iter=0;

    double t_mbasis_rescomp=0.0;
    while (t_mbasis_rescomp<0.1)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, order);

        t1 = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift);
        mbasis_rescomp(appbas,pmat,order,rdeg);
        t2 = GetWallTime();

        t_mbasis_rescomp += t2-t1;
        ++nb_iter;
    }
    t_mbasis_rescomp /= nb_iter;

    double t_mbasis_resupdate=0.0;
    nb_iter=0;
    while (t_mbasis_resupdate<0.1)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, order);
        VecLong pivdeg;

        t1 = GetWallTime();
        Mat<zz_pX> appbas;
        VecLong rdeg(shift);
        mbasis_resupdate(appbas,pmat,order,rdeg);
        t2 = GetWallTime();

        t_mbasis_resupdate += t2-t1;
        ++nb_iter;
    }
    t_mbasis_resupdate /= nb_iter;

    double t_mbasis_int_rescomp=0.0;
    nb_iter=0;
    while (t_mbasis_int_rescomp<0.1)
    {
        Vec<zz_p> pts(INIT_SIZE, order);
        random(pts, order);
        Vec<Mat<zz_p>> evals(INIT_SIZE, order);
        for (long i = 0; i < order; ++i)
            random(evals[i], rdim, cdim);

        t1 = GetWallTime();
        Mat<zz_pX> intbas;
        mbasis_rescomp(intbas,evals,pts,shift,0,order);
        t2 = GetWallTime();

        t_mbasis_int_rescomp += t2-t1;
        ++nb_iter;
    }
    t_mbasis_int_rescomp /= nb_iter;

    double t_mbasis_int_resupdate=0.0;
    nb_iter=0;
    while (t_mbasis_int_resupdate<0.1)
    {
        Vec<zz_p> pts(INIT_SIZE, order);
        random(pts, order);
        Vec<Mat<zz_p>> evals(INIT_SIZE, order);
        for (long i = 0; i < order; ++i)
            random(evals[i], rdim, cdim);

        t1 = GetWallTime();
        Mat<zz_pX> intbas;
        mbasis_resupdate(intbas,evals,pts,shift,0,order);
        t2 = GetWallTime();

        t_mbasis_int_resupdate += t2-t1;
        ++nb_iter;
    }
    t_mbasis_int_resupdate /= nb_iter;


    cout << rdim << "\t" << cdim << "\t" << order;
    cout << "\t" << t_mbasis_rescomp << "\t" << t_mbasis_resupdate;
    cout << "\t" << t_mbasis_int_rescomp << "\t" << t_mbasis_int_resupdate;
    cout << endl;
}

/*------------------------------------------------------------*/
/* run bench on variety of parameters                         */
/*------------------------------------------------------------*/
void run_bench(long nthreads, long nbits, bool fftprime)
{
    SetNumThreads(nthreads);

    if (fftprime)
    {
        cout << "Bench mbasis, FFT prime p = ";
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

    VecLong szs = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

    cout << "rdim\tcdim\torder\tapp-rcmp\tapp-rupd\tint-rcmp\tint-rupd" << endl;
    for (size_t i=0; i<szs.size(); ++i)
    {
        long interval = ceil( (double)szs[i] / 20);
        for (long j=1; j<szs[i]; j+=interval)
        {
            long max_order=128;
            if (szs[i]==512)
                max_order=64;
            else if (szs[i]==1024)
                max_order=32;
            for (long k=2; k<=max_order; k=2*k)
                one_bench_mbasis(szs[i],j,k);
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
        throw std::invalid_argument("Usage: ./time_mbasis nbits fftprime");

    long nbits = atoi(argv[1]);
    bool fftprime = (atoi(argv[2])==1);

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
