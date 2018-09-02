#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/********************************************
 *  tests the determinant algorithms        *
 ********************************************/

/*------------------------------------------------------------*/
/* tests the determinant algorithms                           */
/*------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    SetNumThreads(4);

    bool verify=false;

    if (argc!=4 && argc!=5)
        throw std::invalid_argument("Usage: ./test_determinant rdim degree nbits (verify)");

    long rdim = atoi(argv[1]);
    long degree = atoi(argv[2]);
    long nbits = atoi(argv[3]);
    if (argc==5)
        verify = (atoi(argv[4])==1);

    if (nbits==0)
        zz_p::FFTInit(0);
        //zz_p::UserFFTInit(23068673);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Testing determinant with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus();
    if (nbits==0) std::cout << "  (FFT prime)";
    std::cout << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--degree =\t" << degree << std::endl;

    double t1,t2;

    // build random matrix
    Mat<zz_pX> pmat;
    t1 = GetTime();
    random_mat_zz_pX(pmat, rdim, rdim, degree+1);
    t2 = GetTime();

    std::cout << "Time(random mat creation): " << (t2-t1) << "\n";

    std::cout << "warming up..." << std::endl;
    //warmup();

    // generic case
    {
        std::cout << "~~~computing determinant, algorithm for generic case with known degree~~~" << std::endl;
        t1 = GetTime();
        zz_pX det;
        bool b = determinant_generic_knowing_degree(det, pmat, rdim*degree);
        t2 = GetTime();

        std::cout << "Time(determinant computation): " << (t2-t1) << "s\n";
        std::cout << "Issue detected during computation? " << (b?"yes":"no") << std::endl;

        if (verify)
        {
            std::cout << "Verifying computed determinant:" << std::endl;
            t1 = GetTime();
            bool correct = verify_determinant(det, pmat, true, true);
            t2 = GetTime();
            std::cout << "Time(verification): " << (t2-t1) << "\n";
            std::cout << "Correctness: " << (correct?"correct":"wrong") << std::endl;

            if (rdim*rdim*degree<200)
            {
                std::cout << "Input matrix:" << std::endl;
                std::cout << pmat << std::endl;
                std::cout << "Output determinant:" << std::endl;
                std::cout << det << std::endl;
            }
        }
    }

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
