#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

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
    if (argc!=6)
        throw std::invalid_argument("Usage: ./test_matpadegen rdim order nbits verify iter");

    long rdim = atoi(argv[1]);
    long order = atoi(argv[2]);
    long nbits = atoi(argv[3]);
    bool verify = (atoi(argv[4])==1);
    bool iter = (atoi(argv[5])==1);

    VecLong shift(2*rdim,0);

    if (nbits==0)
        zz_p::FFTInit(0);
        //zz_p::UserFFTInit(786433);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Testing matrix_pade_generic with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--order <\t" << order << std::endl;

    long deg_gcd = (order>>1);

    // GCD computation, for reference
    if (rdim==1) // 2 x 1 pmbasis below
    {
        double t, tt;
        long nb_iter=0;
        std::cout << "For reference, timings for GCD computation (degree " << deg_gcd << "):" << std::endl;
        //{
        //    t = 0.;
        //    while (t < 0.2)
        //    {
        //        zz_pX a,b,g;
        //        random(a, deg_gcd);
        //        random(b, deg_gcd);
        //        tt = GetWallTime();
        //        NTL::GCD(g, a, b);
        //        t += GetWallTime()-tt;
        //        ++nb_iter;
        //    }
        //    std::cout << "\t GCD --> " << t/nb_iter << std::endl;
        //}
        //{
        //    nb_iter=0;
        //    t = 0.;
        //    while (t < 0.2)
        //    {
        //        zz_pX a,b,g;
        //        random(a, deg_gcd);
        //        random(b, deg_gcd);
        //        tt = GetWallTime();
        //        NTL::PlainGCD(g, a, b);
        //        t += GetWallTime()-tt;
        //        ++nb_iter;
        //    }
        //    std::cout << "\t GCD --> " << t/nb_iter << std::endl;
        //}
        {
            nb_iter=0;
            t = 0.;
            while (t < 0.2)
            {
                zz_pX a,b,g,u,v; 
                random(a, deg_gcd);
                random(b, deg_gcd);
                tt = GetWallTime();
                NTL::XGCD(g, u, v, a, b);
                t += GetWallTime()-tt;
                ++nb_iter;
            }
            std::cout << "\tXGCD --> " << t/nb_iter << std::endl;
        }
    }


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


    Mat<zz_pX> appbas;

    if (iter)
    {
        // via mbasis_generic_2n_n_resupdate
        while (tt<0.5)
        {
            t1 = GetWallTime();
            mbasis_generic_2n_n_resupdate(appbas,pmat_approx,order);
            t2 = GetWallTime();
            tt += t2-t1;
            ++nb_iter;
        }
        std::cout << "Time(mbasis_generic_2n_n_resupdate): " << tt/nb_iter << ", ";

        if (verify)
        {
            bool verif = is_approximant_basis(appbas,pmat_approx,order,shift,ORD_WEAK_POPOV,true);
            std::cout << (verif?"correct":"wrong");
        }
        std::cout << std::endl;
    }

    // via general pmbasis
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        t1 = GetWallTime();
        pmbasis(appbas,pmat_approx,order,shift);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    std::cout << "Time(pmbasis): " << tt/nb_iter << ", ";

    if (verify)
    {
        bool verif = is_approximant_basis(appbas,pmat_approx,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct":"wrong");
    }
    std::cout << std::endl;


    // via pmbasis_generic_2n_n
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        t1 = GetWallTime();
        pmbasis_generic_2n_n(appbas,pmat_approx,order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    std::cout << "Time(pmbasis_generic_2n_n): " << tt/nb_iter << ", ";

    if (verify)
    {
        bool verif = is_approximant_basis(appbas,pmat_approx,order,shift,ORD_WEAK_POPOV,true);
        std::cout << (verif?"correct":"wrong");
    }
    std::cout << std::endl;

    if (iter)
    {
        // timing matrix_pade_generic_iterative
        tt=0.0; nb_iter=0;
        Mat<zz_pX> den1, den2;
        while (tt<0.5)
        {
            t1 = GetWallTime();
            matrix_pade_generic_iterative(den1, den2, pmat, order);
            t2 = GetWallTime();
            tt += t2-t1;
            ++nb_iter;
        }

        std::cout << "time(matrix_pade_generic_iterative): " << tt/nb_iter << ", ";

        if (verify)
        {
            bool verif = true;
            for (long i = 0; i < rdim; ++i)
            {
                for (long j = 0; j < rdim; ++j)
                    if (appbas[i][j] != den1[i][j])
                    {
                        verif = false;
                        break;
                    }
                if (not verif)
                    break;
            }
            std::cout << (verif?"correct":"wrong");
        }
        std::cout << std::endl;
    }

    // timing matrix_pade_generic
    tt=0.0; nb_iter=0;
    Mat<zz_pX> den;
    while (tt<0.5)
    {
        t1 = GetWallTime();
        matrix_pade_generic(den, pmat, order);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }

    std::cout << "time(matrix_pade_generic): " << tt/nb_iter << ", ";

    if (verify)
    {
        bool verif = true;
        for (long i = 0; i < rdim; ++i)
        {
            for (long j = 0; j < rdim; ++j)
                if (appbas[i][j] != den[i][j])
                {
                    verif = false;
                    break;
                }
            if (not verif)
                break;
        }
        std::cout << (verif?"correct":"wrong");
    }
    std::cout << std::endl;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
