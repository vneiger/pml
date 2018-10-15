#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"

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
        throw std::invalid_argument("Usage: ./test_appbas rdim cdim order nbits verify nthreads");

    long rdim   = atoi(argv[1]);
    long cdim   = atoi(argv[2]);
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
    Shift shift7(rdim); // "big step" shift
    for (long i = 0; i < rdim; ++i)
        if (i>=rdim/2)
            shift7[i] = rdim*cdim*order;

    std::vector<Shift> shifts = {shift1, shift2, shift3, shift4, shift5, shift6, shift7};

    std::cout << "Testing approximant basis computation (iterative) with random input matrix" << std::endl;
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
    t2w =  GetWallTime(); t2 = GetTime();
    //std::cout << "Time(random mat creation): " << (t2w-t1w) <<  "s,  " << (t2-t1) << "s\n";

    warmup();

    for (Shift shift : shifts)
    {
        if (rdim<40)
            std::cout << "--shift =\t" << shift << std::endl;
        else
            std::cout << "--shift =\t" << "<length " << rdim << ">" << std::endl;

        {
            std::cout << "~~~Iterative approximant basis (order-wise)~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            std::vector<long> orders(cdim,order);
            std::vector<long> pivdeg = appbas_iterative(appbas,pmat,orders,shift,true);
            t2w =    GetWallTime(); t2 = GetTime();

            std::cout << "Time(appbas computation): " << (t2w-t1w) << "s,    " << (t2-t1) <<    "s\n";

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas,pmat,orders,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,    " << (t2-t1) << "s\n";

                Mat<zz_pX> residual;
                multiply_naive(residual,appbas,pmat);

                if (rdim*cdim*order < 100)
                {
                    std::cout << "Print output approx basis..." << std::endl;
                    std::cout << appbas << std::endl;
                    std::cout << "Print final residual..." << std::endl;
                    std::cout << residual << std::endl;
                }

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
        }

        {
            std::cout << "~~~Iterative approximant basis (column-wise)~~~" << std::endl;
            t1w = GetWallTime(); t1 = GetTime();
            Mat<zz_pX> appbas;
            std::vector<long> orders(cdim,order);
            std::vector<long> pivdeg = appbas_iterative(appbas,pmat,orders,shift,false);
            t2w =    GetWallTime(); t2 = GetTime();

            std::cout << "Time(appbas computation): " << (t2w-t1w) << "s,    " << (t2-t1) <<    "s\n";

            if (verify)
            {
                std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
                t1w = GetWallTime(); t1 = GetTime();
                bool verif = is_approximant_basis(appbas,pmat,orders,shift,ORD_WEAK_POPOV,true,false);
                t2w = GetWallTime(); t2 = GetTime();
                std::cout << (verif?"correct":"wrong") << std::endl;
                std::cout << "Time(verification): " << (t2w-t1w) << "s,    " << (t2-t1) << "s\n";

                Mat<zz_pX> residual;
                multiply_naive(residual,appbas,pmat);

                if (rdim*cdim*order < 100)
                {
                    std::cout << "Print output approx basis..." << std::endl;
                    std::cout << appbas << std::endl;
                    std::cout << "Print final residual..." << std::endl;
                    std::cout << residual << std::endl;
                }

                if (std::max(rdim,cdim)<33) {
                    Mat<long> degmat;
                    degree_matrix(degmat,appbas,shift,true);
                    std::cout << "Print degree matrix of approx basis..." << std::endl;
                    std::cout << degmat << std::endl;
                }
            }
            std::cout << std::endl;
        }


        //{
        //    std::cout << "~~~Iterative Popov approximant basis~~~" << std::endl;
        //    t1w = GetWallTime(); t1 = GetTime();
        //    Mat<zz_pX> appbas;
        //    std::vector<long> orders(cdim,order);
        //    std::vector<long> pivdeg = popov_appbas_iterative(appbas,pmat,orders,shift,true);
        //    t2w =    GetWallTime(); t2 = GetTime();

        //    std::cout << "Time(appbas computation): " << (t2w-t1w) << "s,    " << (t2-t1) <<    "s\n";

        //    if (verify)
        //    {
        //        std::cout << "Verifying Popov approximant basis..." << std::endl;
        //        t1w = GetWallTime(); t1 = GetTime();
        //        bool verif = is_approximant_basis(appbas,pmat,orders,shift,POPOV,true,false);
        //        t2w =    GetWallTime(); t2 = GetTime();
        //        std::cout << (verif?"correct":"wrong") << std::endl;
        //        std::cout << "Time(verification): " << (t2w-t1w) << "s,    " << (t2-t1) << "s\n";

        //        Mat<zz_pX> residual;
        //        multiply_naive(residual,appbas,pmat);

        //        if (rdim*cdim*order < 100)
        //        {
        //            std::cout << "Print output approx basis..." << std::endl;
        //            std::cout << appbas << std::endl;
        //            std::cout << "Print final residual..." << std::endl;
        //            std::cout << residual << std::endl;
        //        }

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
