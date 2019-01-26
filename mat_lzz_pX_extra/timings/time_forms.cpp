#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>

#include "util.h"
#include "mat_lzz_pX_extra.h"

#define TIME(a)         \
    tt = get_time();    \
    a                   \
    tt = get_time()-tt; \
    cout << tt << "\t";

#define SMALL_SUITE


NTL_CLIENT

/******************************************************/
/* Times the degree functions for polynomial matrices */
/******************************************************/

using namespace std;

void run_bench(long m, long n, long d)
{
    double tt;

    TIME(
         Mat<zz_pX> pmat;
         random(pmat, m, n, d);
        )

}


int main()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);

#ifdef SMALL_SUITE
    VecLong rdims = {5,200,1000,};
    VecLong cdims = {5,200,1000,};
    VecLong degs = {5,500,5000,};
    VecLong nbits = {20,60};
    VecLong fftprimes = {786433,1139410705724735489,}; // 20, 60 bits
#else
    VecLong rdims = {5,10,20,40,70,100,150,200,400,1000,};
    VecLong cdims = {5,10,20,40,70,100,150,200,400,1000,};
    VecLong degs = {5,10,20,40,70,100,150,200,400,1000,2000,4000,8000,16000,};
    VecLong nbits = {20,31,42,60};
    VecLong fftprimes = {786433,2013265921,2748779069441,1139410705724735489,}; // 20, 31, 42, 60 bits
#endif // SMALL_SUITE

    warmup();

    for (size_t b = 0; b < nbits.size(); ++b)
    {
        std::cout << "Bench utils, nbits = " << nbits[b];
        long p = NTL::GenPrime_long(nbits[b]);
        std::cout << " (normal prime p = " << p << ", fft prime p = " << fftprimes[b] << ")" << std::endl;

        std::cout << "info\trdim\tcdim\tdeg\trand\t" << std::endl;

        for (size_t i = 0; i < rdims.size(); ++i)
            for (size_t j = 0; j < cdims.size(); ++j)
                for (size_t k = 0; k < degs.size(); ++k)
                {
                    std::cout << "\t" << rdims[i] << "\t" << cdims[j] << "\t" << degs[k] << "\t";
                    zz_p::init(p);
                    run_bench(rdims[i], cdims[j], degs[k]);
                    std::cout << std::endl << "fft\t" << rdims[i] << "\t" << cdims[j] << "\t" << degs[k] << "\t";
                    zz_p::UserFFTInit(fftprimes[b]);
                    run_bench(rdims[i], cdims[j], degs[k]);
                    std::cout << std::endl;
                }
    }

    std::cout << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
