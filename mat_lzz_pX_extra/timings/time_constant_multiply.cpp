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

NTL_CLIENT

/************************************************************************/
/* Times the multiplication of a polynomial matrix by a constant matrix */
/************************************************************************/

using namespace std;

void run_bench(long m, long n, long k, long d)
{
    double tt;

    TIME(
         Mat<zz_pX> pmat1;
         Mat<zz_p> pmat2;
         random(pmat1, m, n, d);
         random(pmat2, n, k);
        )

    TIME(
         Mat<zz_pX> pmat3;
         mul(pmat3, pmat1, pmat2);
        )

    TIME(
         Vec<Mat<zz_p>> pmat1_conv;
         conv(pmat1_conv, pmat1);
        )

    TIME(
         Vec<Mat<zz_p>> pmat3_conv;
         pmat3_conv.SetLength(d);
         for (long dd = 0; dd < d; ++dd)
            mul(pmat3_conv[dd], pmat1_conv[dd], pmat2);
        )

    TIME(
         Mat<zz_p> pmat1_const;
         pmat1_const.SetDims(m*d,n);
         for (long dd = 0; dd < d; ++dd)
            for (long i = 0; i < m; ++i)
                for (long j = 0; j < n; ++j)
                    pmat1_const[dd*m+i][j] = pmat1[i][j][dd];
        )

    TIME(
         Mat<zz_p> pmat3_const;
         mul(pmat3_const, pmat1_const, pmat2);
        )
}


int main(int argc, char* argv[])
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    VecLong rdims = {5,50,100,200,1000};
    VecLong cdims = {5,50,100,200,1000};
    VecLong rhsdims = {5,50,100,200,1000};
    VecLong degs = {5,10,20,50,100,200,500,1000};
    VecLong nbits = {20,30,40,50,60};

    if ((argc<6 && argc>1) || argc > 6)
    {
        std::cout << "ERROR. Usage: './time_constant_multiply' or './time_constant_multiply rdim cdim rhsdim deg nbits'" << std::endl;
        return 0;
    }
    if (argc==6)
    {
        rdims = {atoi(argv[1])};
        cdims = {atoi(argv[2])};
        rhsdims = {atoi(argv[3])};
        degs = {atoi(argv[4])};
        nbits = {atoi(argv[5])};
    }

    warmup();

    for (size_t b = 0; b < nbits.size(); ++b)
    {
        std::cout << "Bench poly-mat x constant-mat multiplication, nbits = " << nbits[b];
        long p = NTL::GenPrime_long(nbits[b]);
        std::cout << " (normal prime p = " << p << ")" << std::endl;

        std::cout << "info\trdim\tcdim\trhsdim\tdeg\tnbits\trand\tmul\tconv\tcvmul\tconst\tcstmul" << std::endl;

        for (size_t i = 0; i < rdims.size(); ++i)
            for (size_t j = 0; j < cdims.size(); ++j)
                for (size_t k = 0; k < rhsdims.size(); ++k)
                    for (size_t d = 0; d < degs.size(); ++d)
                    {
                        std::cout << "\t" << rdims[i] << "\t" << cdims[j] << "\t" << rhsdims[k] << "\t" << degs[d] << "\t" << nbits[b] << "\t";
                        zz_p::init(p);
                        run_bench(rdims[i], cdims[j], rhsdims[k], degs[d]);
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
