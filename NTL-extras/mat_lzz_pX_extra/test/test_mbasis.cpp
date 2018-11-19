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
        throw std::invalid_argument("Usage: ./test_mbasis nbits nthreads");

    long nbits = atoi(argv[1]);
    SetNumThreads(atoi(argv[2]));

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    if (false)
    {
        zz_p::init(3);
        Mat<zz_pX> pmat;
        pmat.SetDims(3,2);
        //[[[] [2 2]]
        //[[1 1] [0 1]]
        //[[1 2] [1 2]]
        SetCoeff(pmat[1][0], 0, 1);
        SetCoeff(pmat[1][0], 1, 1);
        SetCoeff(pmat[2][0], 0, 1);
        SetCoeff(pmat[2][0], 1, 2);
        SetCoeff(pmat[0][1], 0, 2);
        SetCoeff(pmat[0][1], 1, 2);
        SetCoeff(pmat[1][1], 1, 1);
        SetCoeff(pmat[2][1], 0, 1);
        SetCoeff(pmat[2][1], 1, 2);

        VecLong shift = {0,1,2};
        long order = 5;
        Mat<zz_pX> appbas;
        VecLong pivdeg = mbasis_rescomp(appbas, pmat, order, shift);

        std::cout << "pivdeg: " << pivdeg << std::endl;
        std::cout << "appbas: " << std::endl << appbas << std::endl;
//]
//[[[1 1 2 1 1] [2 1] [1 1]]
//[[0 0 1 0 2] [0 0 2 2 1] [0 0 1]]
//[[0 1 1 2 1] [0 2 1] [0 1 1]]
//]
        is_approximant_basis(appbas, pmat, order, shift, ORD_WEAK_POPOV, true);
        return 0;
    }


    // build couple (test_matrices, test_shifts)
    std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
    test_examples = build_test_examples();

    std::cout << "Testing approximant basis computation (mbasis)." << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    VecLong pivdeg; 
    Mat<zz_pX> appbas;

    size_t i=0;
    for (auto pmat = test_examples.first.begin(); pmat!= test_examples.first.end(); ++pmat, ++i)
    {
#ifdef VERBOSE
        std::cout << i << std::endl;
#endif // VERBOSE
        
        long rdim = pmat->NumRows();
        long cdim = pmat->NumCols();
        long d = deg(*pmat);
        VecLong orders;
        if (d <= 1)
            orders = {1, 2, 5}; // recall the order must be (strictly) positive
        else
            orders = {d/2,d+1,2*d};

        for (long order : orders)
        {
            for (VecLong shift : test_examples.second[i])
            {
#ifdef VERBOSE
                std::cout << "--rdim =\t" << rdim << std::endl;
                std::cout << "--cdim =\t" << cdim << std::endl;
                std::cout << "--deg =\t" << d << std::endl;
                std::cout << "--order =\t" << order << std::endl;
                std::cout << "--shift =\t" << shift << std::endl;
#endif // VERBOSE

//#ifdef VERBOSE
//                std::cout << "Computation popov_mbasis1... ";
//#endif // VERBOSE
//                // popov_mbasis1
//                Mat<zz_p> kerbas;
//                pivdeg = popov_mbasis1(kerbas,coeff(*pmat,0),shift);
//
//                // build approx basis from kerbas
//                clear(appbas);
//                appbas.SetDims(rdim,rdim);
//                long row=0;
//                for (long r = 0; r < rdim; ++r)
//                {
//                    if (pivdeg[r]==0)
//                    {
//                        for (long j = 0; j < rdim; ++j)
//                            appbas[r][j] = kerbas[row][j];
//                        ++row;
//                    }
//                    else
//                        SetX(appbas[r][r]);
//                }
//#ifdef VERBOSE
//                std::cout << "OK. Testing... ";
//#endif // VERBOSE
//
//                if (not is_approximant_basis(appbas,*pmat,1,shift,ORD_WEAK_POPOV,true))
//                {
//                    std::cout << "Error in popov_mbasis1." << std::endl;
//                    std::cout << "--rdim =\t" << rdim << std::endl;
//                    std::cout << "--cdim =\t" << cdim << std::endl;
//                    std::cout << "--deg =\t" << d << std::endl;
//                    std::cout << "--order =\t" << order << std::endl;
//                    std::cout << "--shift =\t" << shift << std::endl;
//                    std::cout << zz_p::modulus() << std::endl;
//                    std::cout << coeff(*pmat,0) << std::endl;
//                    std::cout << appbas << std::endl;
//                    std::cout << kerbas << std::endl;
//                    std::cout << pivdeg << std::endl;
//                    return 0;
//                }
//#ifdef VERBOSE
//                std::cout << "OK." << std::endl;
//#endif // VERBOSE
//
//                // plain mbasis
//#ifdef VERBOSE
//                std::cout << "Computation mbasis_plain... ";
//#endif // VERBOSE
//                pivdeg = mbasis_plain(appbas,*pmat,order,shift);
//#ifdef VERBOSE
//                std::cout << "OK. Testing... ";
//#endif // VERBOSE
//                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
//                {
//                    std::cout << "Error in mbasis_plain." << std::endl;
//                    std::cout << "--rdim =\t" << rdim << std::endl;
//                    std::cout << "--cdim =\t" << cdim << std::endl;
//                    std::cout << "--deg =\t" << d << std::endl;
//                    std::cout << "--order =\t" << order << std::endl;
//                    std::cout << "--shift =\t" << shift << std::endl;
//                    std::cout << zz_p::modulus() << std::endl;
//                    std::cout << *pmat << std::endl;
//                    std::cout << appbas << std::endl;
//                    return 0;
//                }
//#ifdef VERBOSE
//                std::cout << "OK." << std::endl;
//#endif // VERBOSE

                // mbasis Vec<Mat<zz_p>>, rescomp
#ifdef VERBOSE
                std::cout << "Computation mbasis_rescomp... ";
#endif // VERBOSE
                pivdeg = mbasis_rescomp(appbas,*pmat,order,shift);
#ifdef VERBOSE
                std::cout << "OK. Testing... ";
#endif // VERBOSE
                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
                {
                    std::cout << "Error in mbasis_rescomp." << std::endl;
                    std::cout << "--rdim =\t" << rdim << std::endl;
                    std::cout << "--cdim =\t" << cdim << std::endl;
                    std::cout << "--deg =\t" << d << std::endl;
                    std::cout << "--order =\t" << order << std::endl;
                    std::cout << "--shift =\t" << shift << std::endl;
                    std::cout << zz_p::modulus() << std::endl;
                    std::cout << *pmat << std::endl;
                    std::cout << appbas << std::endl;
                    return 0;
                }
#ifdef VERBOSE
                std::cout << "OK." << std::endl;
#endif // VERBOSE

//                // mbasis_resupdate version 1
//#ifdef VERBOSE
//                std::cout << "Computation mbasis_resupdate... ";
//#endif // VERBOSE
//                pivdeg = mbasis_resupdate(appbas,*pmat,order,shift);
//#ifdef VERBOSE
//                std::cout << "OK. Testing...";
//#endif // VERBOSE
//                if (not is_approximant_basis(appbas,*pmat,order,shift,ORD_WEAK_POPOV,true))
//                {
//                    std::cout << "Error in mbasis_resupdate." << std::endl;
//                    std::cout << "--rdim =\t" << rdim << std::endl;
//                    std::cout << "--cdim =\t" << cdim << std::endl;
//                    std::cout << "--deg =\t" << d << std::endl;
//                    std::cout << "--order =\t" << order << std::endl;
//                    std::cout << "--shift =\t" << shift << std::endl;
//                    std::cout << zz_p::modulus() << std::endl;
//                    std::cout << *pmat << std::endl;
//                    std::cout << appbas << std::endl;
//                    return 0;
//                }
//#ifdef VERBOSE
//                std::cout << "OK." << std::endl;
//#endif // VERBOSE
//
//                // popov_mbasis
//#ifdef VERBOSE
//                std::cout << "Computation popov_mbasis... ";
//#endif // VERBOSE
//                pivdeg = popov_mbasis(appbas,*pmat,order,shift);
//#ifdef VERBOSE
//                std::cout << "OK. Testing...";
//#endif // VERBOSE
//                if (not is_approximant_basis(appbas,*pmat,order,shift,POPOV,true))
//                {
//                    std::cout << "Error in popov_mbasis." << std::endl;
//                    std::cout << "--rdim =\t" << rdim << std::endl;
//                    std::cout << "--cdim =\t" << cdim << std::endl;
//                    std::cout << "--deg =\t" << d << std::endl;
//                    std::cout << "--order =\t" << order << std::endl;
//                    std::cout << "--shift =\t" << shift << std::endl;
//                    std::cout << zz_p::modulus() << std::endl;
//                    std::cout << *pmat << std::endl;
//                    std::cout << appbas << std::endl;
//                    return 0;
//                }
//#ifdef VERBOSE
//                std::cout << "OK." << std::endl;
//#endif // VERBOSE
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
