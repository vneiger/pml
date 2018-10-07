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
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_bench_mbasis(long rdim, long cdim, long order)
{
    std::vector<long> shift(rdim,0);

    double t1w,t2w;

    // build random matrix
    Mat<zz_pX> pmat;
    t1w = GetWallTime();
    random_mat_zz_pX(pmat, rdim, cdim, order);
    t2w = GetWallTime();

    std::vector<long> pivdeg;

    { // 1 thread, normal mbasis
        SetNumThreads(1);
        t1w = GetWallTime();
        Mat<zz_pX> appbas;
        pivdeg = mbasis(appbas,pmat,order,shift);
        t2w = GetWallTime();
        cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << ", (" << AvailableThreads() << " threads)" << endl;
    }

    { // 2 threads, normal mbasis
        SetNumThreads(2);
        t1w = GetWallTime();
        Mat<zz_pX> appbas;
        pivdeg = mbasis(appbas,pmat,order,shift);
        t2w = GetWallTime();
        cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << ", (" << AvailableThreads() << " threads)" << endl;
    }

    { // 3 threads, normal mbasis
        SetNumThreads(3);
        t1w = GetWallTime();
        Mat<zz_pX> appbas;
        pivdeg = mbasis(appbas,pmat,order,shift);
        t2w = GetWallTime();
        cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << ", (" << AvailableThreads() << " threads)" << endl;
    }

    { // 4 threads, normal mbasis
        SetNumThreads(4);
        t1w = GetWallTime();
        Mat<zz_pX> appbas;
        pivdeg = mbasis(appbas,pmat,order,shift);
        t2w = GetWallTime();
        cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << ", (" << AvailableThreads() << " threads)" << endl;
    }

}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    std::vector<long> szs = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 };

    // TODO try with non-FFT primes
    std::cout << "Bench pm-basis (FFT prime)" << std::endl;
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
    {
        one_bench_mbasis(szs[i],1,szs[i]);
        one_bench_mbasis(szs[i],szs[i]/2,1);
        one_bench_mbasis(szs[i],szs[i]/2,8);
        one_bench_mbasis(szs[i],szs[i]/2,16);
        one_bench_mbasis(szs[i],szs[i]/2,32);
        if (szs[i] < 1025)
            one_bench_mbasis(szs[i],szs[i]/2,64);
        if (szs[i] < 513)
            one_bench_mbasis(szs[i],szs[i]/2,128);
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
    
    cout << "rdim,cdim,order,time" << endl;

    if (argc==1)
    {
        warmup();
        run_bench(60);
    }
    if (argc==2)
    {
        long nbits = atoi(argv[1]);
        warmup();
        run_bench(nbits);
    }
    if (argc>3)
        throw std::invalid_argument("Usage: ./time_mbasis OR ./time_mbasis nbits");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
