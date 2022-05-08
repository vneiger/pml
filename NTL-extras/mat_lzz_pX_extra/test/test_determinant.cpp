#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_determinant.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests the determinant algorithms                           */
/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc!=3 && argc!=4)
        throw std::invalid_argument("Usage: ./test_determinant rdim degree (fftprime?)");

    const long rdim = atoi(argv[1]);
    const long degree = atoi(argv[2]);
    const bool fftprime = (argc==3) ? false : (atoi(argv[3])==1);

    if (fftprime)
        zz_p::FFTInit(0);
    else
        //zz_p::init(NTL::GenPrime_long(60));
        zz_p::init(3);

    std::cout << "Testing determinant with random input matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus() << (fftprime?"  (FFT prime)":"") << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--degree =\t" << degree << std::endl;

    double t,tt;
    long nb_iter;

    { // generic case
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_generic_knowing_degree(det, pmat, rdim*degree);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(triangular):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    { // via random linear system
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_linsolve(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(system-solve):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (rdim<=50)
    { // via evaluation, general points
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_evaluation_general(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(ev-general):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (rdim<=50)
    { // via evaluation, geometric points
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_evaluation_geometric(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(ev-geometric):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (fftprime)
    { // via evaluation, FFT points
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_evaluation_FFT(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(ev-FFT):\t\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (rdim<=6)
    { // naive
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            determinant_expansion_by_minors_rec(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(minors-rec):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (rdim<=4)
    { // naive
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<0.2)
        {
            Mat<zz_pX> pmat;
            random(pmat, rdim, rdim, degree+1);
            tt = GetWallTime();
            zz_pX det;
            determinant_expansion_by_minors(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(minors):\t\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
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
