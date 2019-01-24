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
    if (argc!=8)
        throw std::invalid_argument("Usage: ./test_mbasisgen rdim cdim order nbits verify verbose nthreads");

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long order = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    bool verify = (atoi(argv[5])==1);
    bool verbose = (atoi(argv[6])==1);
    SetNumThreads(atoi(argv[7]));

    VecLong shift(rdim,0);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Testing mbasis_generic_2n_n_rescomp with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    // build random matrix
    double t1,t2,tt;
    tt=0.0;
    Mat<zz_pX> pmat;
    random(pmat, rdim, cdim, order);

    // mbasis_generic_2n_n_rescomp
    Mat<zz_pX> appbas;
    std::cout << "~~~Testing mbasis_generic_2n_n_rescomp~~~" << std::endl;
    long nb_iter=0;
    while (tt<0.5)
    {
        clear(appbas);
        t1 = GetWallTime();
        mbasis_generic_2n_n_rescomp(appbas,pmat,order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    std::cout << "Time(mbasis_generic_2n_n_rescomp): " << tt/nb_iter << std::endl;

    if (verify)
    {
        std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
        bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct":"wrong") << std::endl;
    }

    if (verbose)
    {
        Mat<long> degmat;
        degree_matrix(degmat,appbas);
        std::cout << "Print degree matrix of approx basis..." << std::endl;
        std::cout << degmat << std::endl;
    }

    // mbasis_generic_2n_n_resupdate
    std::cout << "~~~Testing mbasis_generic_2n_n_resupdate~~~" << std::endl;
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        clear(appbas);
        t1 = GetWallTime();
        mbasis_generic_2n_n_resupdate(appbas,pmat,order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    std::cout << "Time(mbasis_generic_2n_n_resupdate): " << tt/nb_iter << std::endl;

    if (verify)
    {
        std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
        bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct":"wrong") << std::endl;
    }

    if (verbose)
    {
        Mat<long> degmat;
        degree_matrix(degmat,appbas);
        std::cout << "Print degree matrix of approx basis..." << std::endl;
        std::cout << degmat << std::endl;
    }

    // timing mbasis_rescomp, for comparison
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        clear(appbas);
        t1 = GetWallTime();
        mbasis_rescomp(appbas,pmat,order,shift);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }

    if (verify)
    {
        std::cout << "mbasis_rescomp: ";
        bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct, ":"wrong, ");
    }

    std::cout << "time(mbasis_rescomp): " << tt/nb_iter << std::endl;

    // timing mbasis_resupdate, for comparison
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        clear(appbas);
        t1 = GetWallTime();
        mbasis_resupdate(appbas,pmat,order,shift);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }

    if (verify)
    {
        std::cout << "mbasis_resupdate: ";
        bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct, ":"wrong, ");
    }

    std::cout << "time(mbasis_resupdate): " << tt/nb_iter << std::endl;

    // timing mbasis_plain, for comparison
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        clear(appbas);
        t1 = GetWallTime();
        mbasis_plain(appbas,pmat,order,shift);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }

    if (verify)
    {
        std::cout << "mbasis_plain: ";
        bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct, ":"wrong, ");
    }

    std::cout << "time(mbasis_plain): " << tt/nb_iter << std::endl;
 
    // timing mbasis_plain, for comparison
    tt=0.0; nb_iter=0;
    {
        clear(appbas);
        t1 = GetWallTime();
        pmbasis_generic_2n_n(appbas,pmat,order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
    if (verify)
        std::cout << (verif?"correct, ":"wrong, ");
    std::cout << "time(pmbasis_gen): " << tt/nb_iter << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
