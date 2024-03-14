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
#include "mat_lzz_pX_approximant.h"
#include "test_examples.h"

//#define VERBOSE

NTL_CLIENT

/********************************************
 *  tests the approximant basis algorithms  *
 ********************************************/

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

int main(int argc, char *argv[])
{
    if (argc!=3)
        throw std::invalid_argument("Usage: ./test_utils nbits nthreads");

    long nbits = atoi(argv[1]);
    SetNumThreads(atoi(argv[2]));

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    // build couple (test_matrices, test_shifts)
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples();

    std::cout << "Testing util functions. (For the moment, just run them on various matrices, to verify it doesn't crash.)" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    long i=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
        long rdim = pmat->NumRows();
        long cdim = pmat->NumCols();
        long d = deg(*pmat);
#ifdef VERBOSE
        std::cout << "=======================================================" << std::endl;
        std::cout << i << std::endl;
        std::cout << "--rdim =\t" << rdim << std::endl;
        std::cout << "--cdim =\t" << cdim << std::endl;
        std::cout << "--deg =\t" << d << std::endl;
#endif

        Mat<zz_pX> pmat2;
        random(pmat2, rdim, cdim, d+1);
#ifdef VERBOSE
        std::cout << "Random: ok." << std::endl;
#endif // VERBOSE

        ident(pmat2, rdim);

        pmat2 = ident_mat_zz_pX(rdim);

        IsZero(*pmat);

#ifdef VERBOSE
        std::cout << "Zero: ok." << std::endl;
#endif // VERBOSE

        IsIdent(*pmat);

        IsIdent(*pmat, rdim);

#ifdef VERBOSE
        std::cout << "Identity: ok." << std::endl;
#endif // VERBOSE

        long dd = deg(*pmat);
        dd += 1; // to avoid warning "unused variable"

        Mat<zz_p> c1;
        GetCoeff(c1, *pmat, dd/2);

        Mat<zz_p> c2;
        c2 = coeff(*pmat, dd/2);

        SetCoeff(*pmat, dd/4, c1);

#ifdef VERBOSE
        std::cout << "Degree, Get/SetCoeff: ok." << std::endl;
#endif // VERBOSE

        transpose(pmat2, *pmat);

        trunc(pmat2, pmat2, 10);

        pmat2 = *pmat;

        trunc(pmat2,pmat2,dd/2);

        trunc(pmat2,pmat2,5*dd);

        pmat2 = *pmat;

        RightShift(pmat2,pmat2,2);

        pmat2 = *pmat;

        RightShift(pmat2,pmat2,dd-3);

        pmat2 = *pmat;

        LeftShift(pmat2,pmat2,2);

        pmat2 = *pmat;

        LeftShift(pmat2,pmat2,dd-4);

        reverse(pmat2, *pmat, dd/2);

        reverse(pmat2, *pmat);

#ifdef VERBOSE
        std::cout << "Transpose, trunc, shift, reverse, eval: ok." << std::endl;
#endif // VERBOSE

        eval(c1, *pmat, random_zz_p());

        Vec<Mat<zz_p>> coeffs;
        conv(coeffs, *pmat);

        conv(pmat2, coeffs);
#ifdef VERBOSE
        std::cout << "Conversion: ok." << std::endl;
#endif // VERBOSE
    }
    std::cout << "End of test." << std::endl;
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
