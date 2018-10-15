// NTLX
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "magma_output.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"

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
void one_bench_pmbasis(long sz, long deg, long nbits)
{
    long rdim = sz*2;
    long cdim = sz;
    long order = 2*deg;

    std::vector<long> shift(rdim,0);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    double t1w,t2w;

    // build random matrix
    Mat<zz_pX> pmat;
    t1w = GetWallTime();
    random(pmat, rdim, cdim, order);
    t2w = GetWallTime();

    std::vector<long> pivdeg;

    // GCD computation, for reference
    if (rdim==2 && cdim==1)
    {
        long deg_gcd = (order>>1);
        zz_pX a,b,g,u,v; 
        random(a, deg_gcd);
        random(b, deg_gcd);
        t1w = GetWallTime();
        NTL::XGCD(g, u, v, a, b);
        t2w = GetWallTime();
        cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << " (ntl xgcd)" << endl;
    }

    // pmbasis
    {
        //std::cout << "~~~Testing pmbasis~~~" << std::endl;
        t1w = GetWallTime();
        Mat<zz_pX> appbas;
        pivdeg = pmbasis(appbas,pmat,order,shift);
        t2w = GetWallTime();

        cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << endl;
        //cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << " (pmbasis)" << endl;
    }
    // popov_pmbasis
    //{
    //    //std::cout << "~~~Testing popov_pmbasis~~~" << std::endl;
    //    t1w = GetWallTime(); t1 = GetTime();
    //    Mat<zz_pX> appbas;
    //    pivdeg = popov_pmbasis(appbas,pmat,order,shift);
    //    t2w = GetWallTime(); t2 = GetTime();

    //    cout << rdim << "," << cdim << "," << order << "," << (t2w-t1w) << "(povov_pmbasis)" << endl;   
    //}
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    std::vector<long> szs =
    {
        1,1,1,1,1,1,1,1,1,1,1,1,
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
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,16384,32768,
        32,64,128,256,512,1024,2048,4096,8192,
        32,64,128,256,512,1024,2048,
        32,64,128,256,512,
        32,64,128,
        32,64,
    };

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
        one_bench_pmbasis(szs[i],degs[i],nbits);
    cout << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(1);

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
        throw std::invalid_argument("Usage: ./time_pmbasis OR ./time_pmbasis nbits");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
