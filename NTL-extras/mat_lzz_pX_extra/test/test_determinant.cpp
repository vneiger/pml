#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>


#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests the determinant algorithms                           */
/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{

    bool verify=false;

    if (argc!=6)
        throw std::invalid_argument("Usage: ./test_determinant rdim degree nbits nthreads verify");

    long rdim = atoi(argv[1]);
    long degree = atoi(argv[2]);
    long nbits = atoi(argv[3]);
    SetNumThreads(atoi(argv[4]));
    verify = (atoi(argv[5])==1);

    if (nbits==0)
        zz_p::FFTInit(0);
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
    t1 = GetWallTime();
    random(pmat, rdim, rdim, degree+1);
    t2 = GetWallTime();

    std::cout << "Time(random mat creation): " << (t2-t1) << "\n";

    std::cout << "warming up..." << std::endl;
    warmup();

    // generic case
    {
        std::cout << "~~~computing determinant, algorithm for generic case with known degree~~~" << std::endl;
        t1 = GetWallTime();
        zz_pX det;
        bool b = determinant_generic_knowing_degree(det, pmat, rdim*degree);
        t2 = GetWallTime();

        std::cout << "Time(determinant-triangular): " << (t2-t1) << "s\n";
        std::cout << "Issue detected during computation? " << (b?"no":"yes") << std::endl;

        if (verify)
        {
            std::cout << "Verifying computed determinant:" << std::endl;
            t1 = GetWallTime();
            bool correct = verify_determinant(det, pmat, true, true);
            t2 = GetWallTime();
            std::cout << "Time(verification): " << (t2-t1) << "\n";
            std::cout << "Correctness: " << (correct?"correct":"wrong") << std::endl;

            if (rdim*rdim*degree<200)
            {
                std::cout << "Input matrix:" << std::endl;
                //std::cout << pmat << std::endl;
                std::cout << "Output determinant:" << std::endl;
                MakeMonic(det);
                std::cout << det << std::endl;
            }
        }
    }

    // via random system, high prec
    {
        std::cout << "~~~computing determinant, via random system solving~~~" << std::endl;
        double t_det=0.0;
        t1 = GetWallTime();
        zz_pX det;
        determinant_via_linsolve(det, pmat);
        t2 = GetWallTime();
        t_det += t2-t1;

        std::cout << "Time(determinant-linsolve: total): " << t_det << "s\n";

        if (verify)
        {
            std::cout << "Verifying computed determinant:" << std::endl;
            t1 = GetWallTime();
            bool correct = verify_determinant(det, pmat, true, true);
            t2 = GetWallTime();
            std::cout << "Time(verification): " << (t2-t1) << "\n";
            std::cout << "Correctness: " << (correct?"correct":"wrong") << std::endl;

            if (rdim*rdim*degree<200)
            {
                std::cout << "Input matrix:" << std::endl;
                //std::cout << pmat << std::endl;
                std::cout << "Output determinant:" << deg(det) << std::endl;
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
