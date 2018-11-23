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
#include "mat_lzz_pX_linearization.h"
#include "mat_lzz_pX_approximant.h"

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
    if (argc!=5)
        throw std::invalid_argument("Usage: ./test_appbas_mbasis rdim cdim order nbits");

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long order = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    SetNumThreads(1);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    VecLong shift(rdim, 0);

    std::cout << "Testing approximant basis computation (pmbasis) with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;
    std::cout << "--nthreads =\t" << AvailableThreads() << std::endl;

    // build random matrix
    double t=0.0;
    double tt;
    long nb_iter=0;

    // GCD computation, for reference
    if (rdim==2 && cdim==1)
    {
        long deg_gcd = order/2;
        {
            while (t<0.5)
            {
                zz_pX a,b,g;
                random(a, deg_gcd);
                random(b, deg_gcd);
                tt = GetWallTime();
                NTL::GCD(g, a, b);
                t += GetWallTime() - tt;
                ++nb_iter;
            }

            std::cout << "ntl-gcd \t" << t/nb_iter << std::endl;
        }
        {
            t = 0.0;
            nb_iter=0;
            while (t<0.5)
            {
                zz_pX a,b,g,u,v; 
                random(a, deg_gcd);
                random(b, deg_gcd);
                tt = GetWallTime();
                NTL::XGCD(g, u, v, a, b);
                t += GetWallTime() - tt;
                ++nb_iter;
            }
            std::cout << "ntl-xgcd\t" << t/nb_iter << std::endl;
        }
    }

    // pmbasis
    //t = 0.0;
    //nb_iter=0;
    //while (t<0.5)
    //{
    //    Mat<zz_pX> pmat;
    //    random(pmat, rdim, cdim, order);
    //    tt = GetWallTime();
    //    Mat<zz_pX> appbas;
    //    pmbasis(appbas,pmat,order,shift);
    //    t += GetWallTime() - tt;
    //    ++nb_iter;
    //}

    std::cout << "pmbasis   \t" << t/nb_iter << std::endl;

    t = 0.0;
    nb_iter=0;
    while (t<0.5)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, order);
        tt = GetWallTime();
        Mat<zz_pX> appbas;
        pmbasis_generic_2n_n(appbas,pmat,order);
        t += GetWallTime() - tt;
        ++nb_iter;
    }

    std::cout << "pmbasis_gen\t" << t/nb_iter << std::endl;

    t = 0.0;
    nb_iter=0;
    while (t<0.5)
    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, order);
        tt = GetWallTime();
        Mat<zz_pX> appbas;
        pmbasis_generic_2n_n_top_rows(appbas,pmat,order);
        t += GetWallTime() - tt;
        ++nb_iter;
    }

    std::cout << "pmbasis gen toprows\t" << t/nb_iter << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
