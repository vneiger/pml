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
#include "mat_lzz_pX_partial_linearization.h"

NTL_CLIENT


/********************************************
 *  tests the approximant basis algorithms  *
 ********************************************/

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

int main(int argc, char *argv[])
{
    if (argc!=7)
        throw std::invalid_argument("Usage: ./test_appbas_mbasis rdim cdim order nbits verify nthreads");

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long order = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    bool verify = (atoi(argv[5])==1);
    SetNumThreads(atoi(argv[6]));

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    // declare shifts
    Shift shift1(rdim,0); // uniform [0,...,0]
    Shift shift2(rdim); // increasing [0,1,2,..,rdim-1]
    std::iota(shift2.begin(), shift2.end(),0);
    Shift shift3(rdim); // decreasing [rdim,..,3,2,1]
    for (long i = 0; i < rdim; ++i)
        shift3[i] = rdim - i;
    Shift shift4(rdim); // random shuffle of [0,1,...,rdim-1]
    std::iota(shift4.begin(), shift4.end(),0);
    std::shuffle(shift4.begin(), shift4.end(), std::mt19937{std::random_device{}()});
    Shift shift5(rdim); // Hermite shift
    for (long i = 0; i < rdim; ++i)
        shift5[i] = rdim*cdim*order*i;
    Shift shift6(rdim); // reverse Hermite shift
    for (long i = 0; i < rdim; ++i)
        shift6[i] = rdim*cdim*order*(rdim-1-i);
    Shift shift7(rdim);
    for (long i = 0; i < rdim; ++i)
        if (i>=rdim/2)
            shift7[i] = rdim*cdim*order;

    std::vector<Shift> shifts = {shift1, shift2, shift3, shift4, shift5, shift6, shift7};

    std::cout << "Testing approximant basis computation (pmbasis) with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    // build random matrix
    double t1,t2,t1w,t2w;
    Mat<zz_pX> pmat;
    t1w = GetWallTime(); t1 = GetTime();
    random(pmat, rdim, cdim, order);
    pmat.SetDims(rdim,cdim);
    t2w =  GetWallTime(); t2 = GetTime();
    //std::cout << "Time(random mat creation): " << (t2w-t1w) <<  "s,  " << (t2-t1) << "s\n";

    DegVec pivdeg;

    warmup();

    // GCD computation, for reference
    if (rdim==2 && cdim==1)
    {
        long deg_gcd = (order>>1);
        std::cout << "For reference, timings for GCD computation (degree " << deg_gcd << "):" << std::endl;
        {
            zz_pX a,b,g;
            random(a, deg_gcd);
            random(b, deg_gcd);
            t1w = GetWallTime();
            NTL::GCD(g, a, b);
            t2w = GetWallTime();
            std::cout << "\t GCD --> " << (t2w-t1w) << std::endl;
        }
        {
            zz_pX a,b,g,u,v; 
            random(a, deg_gcd);
            random(b, deg_gcd);
            t1w = GetWallTime();
            NTL::XGCD(g, u, v, a, b);
            t2w = GetWallTime();
            std::cout << "\tXGCD --> " << (t2w-t1w) << std::endl;
        }
    }

    for (Shift shift : shifts)
    {
        { // pmbasis
            std::cout << "~~~Testing pmbasis~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            pivdeg = pmbasis(appbas,pmat,order,shift);
            t2w = GetWallTime(); t2 = GetTime();
            std::cout << degree_matrix(appbas) << std::endl;

            std::cout << "Time(pmbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        //{ // pmbasis generic
        //    std::cout << "~~~Testing pmbasis generic~~~" << std::endl;
        //    t1w = GetWallTime(); t1 = GetTime();
        //    Mat<zz_pX> appbas;
        //    pivdeg = pmbasis_generic(appbas,pmat,order,shift);
        //    t2w = GetWallTime(); t2 = GetTime();

        //    std::cout << "Time(pmbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

        //    if (verify)
        //    {
        //        std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
        //        t1w = GetWallTime(); t1 = GetTime();
        //        bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true,false);
        //        t2w = GetWallTime(); t2 = GetTime();
        //        std::cout << (verif?"correct":"wrong") << std::endl;
        //        std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";


        //        if (std::max(rdim,cdim)<33) {
        //            Mat<long> degmat;
        //            degree_matrix(degmat,appbas,shift,true);
        //            std::cout << "Print degree matrix of approx basis..." << std::endl;
        //            std::cout << degmat << std::endl;
        //        }
        //    }
        //}

        //{ // popov_pmbasis
        //    std::cout << "~~~Testing popov_pmbasis~~~" << std::endl;
        //    t1w = GetWallTime(); t1 = GetTime();
        //    Mat<zz_pX> appbas;
        //    pivdeg = popov_pmbasis(appbas,pmat,order,shift);
        //    t2w = GetWallTime(); t2 = GetTime();

        //    std::cout << "Time(popov_pmbasis computation): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

        //    if (verify)
        //    {
        //        std::cout << "Verifying Popov approximant basis..." << std::endl;
        //        t1w = GetWallTime(); t1 = GetTime();
        //        bool verif = is_approximant_basis(appbas,pmat,order,shift,POPOV,true,false);
        //        t2w = GetWallTime(); t2 = GetTime();
        //        std::cout << (verif?"correct":"wrong") << std::endl;
        //        std::cout << "Time(verification): " << (t2w-t1w) << "s,  " << (t2-t1) << "s\n";

        //        if (std::max(rdim,cdim)<33) {
        //            Mat<long> degmat;
        //            degree_matrix(degmat,appbas,shift,true);
        //            std::cout << "Print degree matrix of approx basis..." << std::endl;
        //            std::cout << degmat << std::endl;
        //        }
        //    }
        //}
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
