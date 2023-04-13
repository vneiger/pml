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
    if (argc!=3)
        throw std::invalid_argument("Usage: ./test_vs_gcd deg nbits");

    long degree = atoi(argv[1]);
    long nbits = atoi(argv[2]);

    VecLong shift(2,0);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    // GCD computation
    double t, tt;
    long nb_iter=0;
    std::cout << "Timings for NTL's GCD" << degree << std::endl;
    //{
    //    t = 0.;
    //    while (t < 0.2)
    //    {
    //        zz_pX a,b,g;
    //        random(a, degree);
    //        random(b, degree);
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
    //        random(a, degree);
    //        random(b, degree);
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
            random(a, degree);
            random(b, degree);
            tt = GetWallTime();
            NTL::XGCD(g, u, v, a, b);
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        std::cout << "\tXGCD --> " << t/nb_iter << std::endl;
    }

    // build random matrix
    double t1,t2;
    tt=0.0;
    Mat<zz_pX> pmat;
    random(pmat, 2, 1, degree);
    Mat<zz_pX> appbas;

    // via general pmbasis
    tt=0.0; nb_iter=0;
    while (tt<0.5)
    {
        t1 = GetWallTime();
        pmbasis(appbas,pmat,2*degree-1,shift);
        t2 = GetWallTime();
        tt += t2-t1;
        ++nb_iter;
    }
    std::cout << "Time(pmbasis): " << tt/nb_iter << std::endl;

//    // slightly optimized pmbasis
//    Mat<zz_pX> polys;
//    random(polys, 2, 1, degree);
//    tt=0.0; nb_iter=0;
//    while (tt<0.5)
//    {
//        t1 = GetWallTime();
//
//        Mat<zz_pX> appbas1,appbas2,residual;
//        zz_pX uu,vv;
//
//        // first recursive call, with 'pmat' and 'shift'
//        pmbasis(appbas1,polys,degree,shift);
//        // shift is now the shifted row degree of appbas,
//        // which is the shift for second call
//
//        // residual = (appbas * pmat * X^-order1) mod X^order2
//        multiply(residual, appbas1, polys);
//        RightShift(residual, residual, degree);
//
//        // second recursive call, with 'residual' and 'rdeg'
//        pmbasis(appbas2,residual,degree-1,shift);
//
//        uu = appbas2[0][0] * appbas1[0][0] + appbas2[0][1] * appbas1[1][0];
//        vv = appbas2[0][0] * appbas1[0][1] + appbas2[0][1] * appbas1[1][1];
//
//        t2 = GetWallTime();
//        tt += t2-t1;
//        ++nb_iter;
//    }
    std::cout << "Time(modified-pmbasis): " << tt/nb_iter << std::endl;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
