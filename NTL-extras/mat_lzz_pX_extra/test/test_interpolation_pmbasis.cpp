#include <iomanip>

#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "lzz_pX_CRT.h"
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
    SetNumThreads(1);

    if (argc!=6)
        throw std::invalid_argument("Usage: ./test_intbas_pmbasis rdim cdim npoints nbits verify");

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

    std::cout << "Testing interpolant basis (pmbasis) with random points/matrix" << std::endl;
    std::cout << "--prime =\t" << zz_p::modulus();
    if (nbits==0) std::cout << "  (FFT prime)";
    std::cout << std::endl;
    std::cout << "--rdim =\t" << rdim << std::endl;
    std::cout << "--cdim =\t" << cdim << std::endl;
    std::cout << "--npoints =\t" << npoints << std::endl;
    std::cout << "--shift = uniform" << std::endl;

    double t1w,t2w;

    {
        // pmbasis for interpolants
        zz_p r;
        random(r);
        zz_pX_Multipoint_Geometric eval(r,zz_p(1),npoints);
        Vec<zz_p> pts(INIT_SIZE, npoints);
        for (long i = 0; i < npoints; ++i)
            eval.get_point(pts[i], i);

        // set up evaluations of pmat
        Vec<Mat<zz_p>> evals;
        evals.SetLength(npoints);
        for (long d = 0; d < npoints; d++)
            random(evals[d], rdim, cdim);

        Vec<Mat<zz_p>> copy_evals(evals); // since this one changes its input
        t1w = GetWallTime();
        Mat<zz_pX> intbas;
        VecLong pivdeg = pmbasis_geometric(intbas,copy_evals,pts,r,shift,0,npoints);
        t2w = GetWallTime();
        std::cout << "pmbasis-geometric, time:\t" << (t2w-t1w);
        if (verify)
        {
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false);
            std::cout << (verif?", correct":", wrong");
        }
        std::cout << std::endl;
    }

    {
        Vec<zz_p> pts;
        random(pts, npoints);
        Vec<Mat<zz_p>> evals;
        evals.SetLength(npoints);
        for (long d = 0; d < npoints; d++)
            random(evals[d], rdim, cdim);

        Vec<Mat<zz_p>> copy_evals(evals); // since this one changes its input
        t1w = GetWallTime();
        Mat<zz_pX> intbas;
        VecLong pivdeg = pmbasis(intbas,copy_evals,pts,shift,0,npoints);
        t2w = GetWallTime();

        std::cout << "pmbasis-general, time:\t\t" << (t2w-t1w);
        if (verify)
        {
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,false);
            std::cout << (verif?", correct":", wrong");
        }
        std::cout << std::endl;
    }

    {
        Mat<zz_pX> pmat;
        random(pmat, rdim, cdim, npoints);
        zz_p r = random_zz_p();
        Vec<zz_p> pts;
        t1w = GetWallTime();
        Mat<zz_pX> intbas;
        VecLong pivdeg = pmbasis_geometric(intbas,pmat,r,npoints,shift,pts);
        t2w = GetWallTime();

        std::cout << "pmbasis-geom, time:\t\t" << (t2w-t1w);
        if (verify)
        {
            zz_pX_Multipoint_Geometric eval(r, npoints);
            Vec<Mat<zz_p>> evals;
            eval.evaluate_matrix(evals, pmat);

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
