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
    if (argc!=7)
        throw std::invalid_argument("Usage: ./test_matpadegen rdim order nbits verify verbose nthreads");

    long rdim = atoi(argv[1]);
    long order = atoi(argv[2]);
    long nbits = atoi(argv[3]);
    bool verify = (atoi(argv[4])==1);
    bool verbose = (atoi(argv[5])==1);
    SetNumThreads(atoi(argv[6]));

    VecLong shift(2*rdim,0);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Testing matrix_pade_generic_iterative with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    // build random matrix
    double t1,t2,tt;
    tt=0.0;
    long nb_iter=0;
    Mat<zz_pX> pmat;
    random(pmat, rdim, rdim, order);
    Mat<zz_pX> pmat_approx;
    pmat_approx.SetDims(2*rdim, rdim);
    for (long i = 0; i < rdim; ++i)
        pmat_approx[i] = pmat[i];
    for (long i = 0; i < rdim; ++i)
        pmat_approx[rdim+i][i] = zz_p(-1);


    // via mbasis_generic_2n_n_resupdate
    std::cout << "~~~Testing matrix Pade via mbasis_generic_2n_n_resupdate~~~" << std::endl;
    Mat<zz_pX> appbas;
    while (tt<0.5)
    {
        clear(appbas);
        t1 = GetWallTime();
        mbasis_generic_2n_n_resupdate(appbas,pmat_approx,order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    std::cout << "Time(mbasis_generic_2n_n_resupdate): " << tt/nb_iter << std::endl;

    if (verify)
    {
        std::cout << "Verifying ordered weak Popov approximant basis..." << std::endl;
        bool verif = is_approximant_basis(appbas,pmat_approx,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct":"wrong") << std::endl;
    }

    if (verbose)
    {
        Mat<long> degmat;
        degree_matrix(degmat,appbas);
        std::cout << "Print degree matrix of approx basis..." << std::endl;
        std::cout << degmat << std::endl;
    }


    // timing matrix_pade_generic
    tt=0.0; nb_iter=0;
    Mat<zz_pX> den1, den2;
    while (tt<0.5)
    {
        clear(den1); clear(den2);
        t1 = GetWallTime();
        matrix_pade_generic_iterative(den1, den2, pmat, order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }

    //if (verify)
    //{
    //    std::cout << "mbasis_plain: ";
    //    bool verif = is_approximant_basis(appbas,pmat,order,shift,ORD_WEAK_POPOV,true);
    //    std::cout << (verif?"correct, ":"wrong, ");
    //}

    std::cout << "time(matrix_pade_generic_iterative): " << tt/nb_iter << std::endl;
    //std::cout << "appbas" << std::endl;
    //std::cout << appbas << std::endl;
    //std::cout << "den1" << std::endl;
    //std::cout << den1 << std::endl;
    //std::cout << "den2" << std::endl;
    //std::cout << den2 << std::endl;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
