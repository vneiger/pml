// NTLX
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s)
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
    std::vector<long> shift(rdim,0);

    double t1,t2;
    double t_mbasis=0.0;

    long nb_iter=0;

    while (t_mbasis<0.1)
    {
        Mat<zz_pX> pmat;
        random_mat_zz_pX(pmat, rdim, cdim, order);
        std::vector<long> pivdeg;

        t1 = GetWallTime();
        Mat<zz_pX> appbas;
        pivdeg = mbasis(appbas,pmat,order,shift);
        t2 = GetWallTime();

        t_mbasis += t2-t1;
        ++nb_iter;
    }

    t_mbasis /= nb_iter;

    cout << rdim << "," << cdim << "," << order << "," << AvailableThreads() << "," << t_mbasis << endl;
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
        else if (nbits < 65)
        {
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
            cout << zz_p::modulus() << ", bit length = " << 60 << endl;
        }
    }
    else
    {
        cout << "Bench mbasis, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    std::vector<long> szs = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

    cout << "rdim,cdim,order,nthreads,time" << endl;
    for (size_t i=0;i<szs.size();i++)
    {
        long interval = ceil( (double)szs[i] / 20);
        for (long j=1; j<szs[i]; j+=interval)
        {
            long max_order=128;
            if (szs[i]==512)
                max_order=64;
            else if (szs[i]==1024)
                max_order=32;
            for (long k=1; k<=max_order; k=2*k)
            {
                one_bench_mbasis(szs[i],j,k);
            }
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

    std::vector<long> nthreads = {1,2,3,4};
    std::vector<long> nbits = {20,30,40,60};
    std::vector<bool> fftprime = {true, false};

    if (argc>=2)
        nthreads = {atoi(argv[1])};
    if (argc>=3)
        nbits = {atoi(argv[2])};
    if (argc==4)
        fftprime = {(atoi(argv[3])==1) ? true : false};
    if (argc>4)
        throw std::invalid_argument("Usage: ./time_mbasis OR ./time_mbasis nthreads OR ./time_mbasis nthreads nbits fftprime");

    for (size_t i = 0; i < nthreads.size(); ++i)
        for (size_t j = 0; j < nbits.size(); ++j)
            for (size_t k = 0; k < fftprime.size(); ++k)
                run_bench(nthreads[i],nbits[j],fftprime[k]);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
