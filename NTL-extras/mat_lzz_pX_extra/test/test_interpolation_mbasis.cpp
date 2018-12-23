#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

#include "util.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"

NTL_CLIENT

/********************************************
 *  tests the interpolant basis algorithms  *
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
    if (argc!=5 && argc!=6)
        throw std::invalid_argument("Usage: ./test_intbas_mbasis rdim cdim npoints nbits verify");

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long npoints = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    bool verify = (atoi(argv[5])==1);

    VecLong shift(rdim,0);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "Testing interpolant basis (mbasis) with random points/matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus();
    if (nbits==0) std::cout << "  (FFT prime)";
    std::cout << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--npoints =\t" << npoints << std::endl;
    std::cout << "--shift = uniform\t" << std::endl;

    double tt,t;
    long nb_iter;

    // will contain random matrix and random points
    Vec<Mat<zz_p>> evals;
    Vec<zz_p> pts;

    warmup();

    { // plain mbasis
        std::cout << "plain mbasis, time:\t";
        nb_iter=0; t=0.0;
        Mat<zz_pX> intbas;
        while (t<0.5)
        {
            evals.SetLength(npoints);
            for (long pt = 0; pt < npoints; ++pt)
                random(evals[pt], rdim, cdim);
            random(pts, npoints);
            tt = GetWallTime();
            mbasis(intbas,evals,pts,shift);
            t += GetWallTime()-tt;
            ++nb_iter;
        }

        std::cout << t/nb_iter;

        if (verify) // checks the last iteration of the while loop
        {
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false);
            std::cout << (verif?", correct":", wrong");
        }
        std::cout << std::endl;
    }

    { // mbasis_rescomp
        std::cout << "mbasis_rescomp, time:\t";
        nb_iter=0; t=0.0;
        Mat<zz_pX> intbas;
        while (t<0.5)
        {
            evals.SetLength(npoints);
            for (long pt = 0; pt < npoints; ++pt)
                random(evals[pt], rdim, cdim);
            random(pts, npoints);
            tt = GetWallTime();
            mbasis_rescomp(intbas,evals,pts,shift);
            t += GetWallTime()-tt;
            ++nb_iter;
        }

        std::cout << t/nb_iter;

        if (verify)
        {
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false);
            std::cout << (verif?", correct":", wrong");
        }
        std::cout << std::endl;
    }


    { // mbasis_resupdate
        std::cout << "mbasis_resupdate, time:\t";
        nb_iter=0; t=0.0;
        Mat<zz_pX> intbas;
        while (t<0.5)
        {
            evals.SetLength(npoints);
            for (long pt = 0; pt < npoints; ++pt)
                random(evals[pt], rdim, cdim);
            random(pts, npoints);
            tt = GetWallTime();
            mbasis_resupdate(intbas,evals,pts,shift);
            t += GetWallTime()-tt;
            ++nb_iter;
        }

        std::cout << t/nb_iter;

        if (verify)
        {
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false);
            std::cout << (verif?", correct":", wrong");
        }
        std::cout << std::endl;
    }

    { // mbasis_rescomp, eval rep
        std::cout << "mbasis_rescomp (eval), time:\t";
        nb_iter=0; t=0.0;
        Vec<Mat<zz_p>> intbas_ev;
        while (t<0.5)
        {
            evals.SetLength(npoints);
            for (long pt = 0; pt < npoints; ++pt)
                random(evals[pt], rdim, cdim);
            random(pts, npoints);
            tt = GetWallTime();
            mbasis_rescomp(intbas_ev,evals,pts,shift);
            t += GetWallTime()-tt;
            ++nb_iter;
        }

        std::cout << t/nb_iter;

        if (verify)
        {
            Mat<zz_pX> intbas;
            intbas.SetDims(rdim, rdim);
            for (long i = 0; i < rdim; ++i)
                for (long j = 0; j < rdim; ++j)
                {
                    Vec<zz_p> evs(INIT_SIZE, npoints);
                    for (long k = 0; k < npoints; ++k)
                        evs[k] = intbas_ev[k][i][j];
                    interpolate(intbas[i][j], pts, evs);
                }
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false);
            std::cout << (verif?", correct":", wrong");
        }
        std::cout << std::endl;
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
