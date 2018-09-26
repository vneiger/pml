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
 *  tests the interpolant basis algorithms  *
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
    SetNumThreads(4);

    bool verify=false;

    if (argc!=5 && argc!=6)
        throw std::invalid_argument("Usage: ./test_intbas_pmbasis rdim cdim npoints nbits (verify)");

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long npoints = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    if (argc==6)
        verify = (atoi(argv[5])==1);

    std::vector<long> shift(rdim,0);
    //std::vector<long> shift {0,1,0,1};
    //std::vector<long> shift {4,1,0,1};
    //std::iota(shift.begin(), shift.end(),0);
    //std::shuffle(shift.begin(), shift.end(), std::mt19937{std::random_device{}()});

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
    std::cout << "--shift =\t";
    if (shift.size()<50)
        std::cout << shift << std::endl; 
    else
        std::cout << "length " << shift.size() << std::endl;

    double t1w,t2w;

    // build random matrix
    Vec<Mat<zz_p>> evals;
    Vec<zz_p> pts;
    t1w = GetWallTime();
    evals.SetLength(npoints);
    for (long pt = 0; pt < npoints; ++pt)
        random(evals[pt], rdim, cdim);
    random(pts, npoints);
    t2w = GetWallTime();

    std::cout << "Time(random matrix/points): " << (t2w-t1w) << std::endl;

    if (pts.length()<50 && zz_p::modulus() < 100)
        std::cout << "Points: " << pts << std::endl;

    // to warm up
    std::cout << "warming up..." << std::endl;
    warmup();

    // generic uniform shift pmbasis for interpolants
    {
        std::vector<long> pivdeg;
        std::cout << "~~~Testing pmbasis (generic, uniform shift)~~~" << std::endl;
        t1w = GetWallTime();
        Mat<zz_pX> intbas;
        pivdeg = pmbasis(intbas,evals,pts,shift);
        t2w = GetWallTime();

        std::cout << "Time(pmbasis-interpolation): " << (t2w-t1w) << std::endl;

        if (verify)
        {
            std::cout << "Verifying ordered weak Popov interpolant basis..." << std::endl;
            t1w = GetWallTime();
            bool verif = is_interpolant_basis(intbas,evals,pts,shift,ORD_WEAK_POPOV,true,false);
            t2w = GetWallTime();
            std::cout << (verif?"correct":"wrong") << std::endl;
            std::cout << "Time(verification): " << (t2w-t1w) << std::endl;

            if (rdim*cdim*npoints < 100)
            {
                std::cout << "Print output interpolant basis..." << std::endl;
                std::cout << intbas << std::endl;
            }

            if (std::max(rdim,cdim)<33) {
                Mat<long> degmat;
                degree_matrix(degmat,intbas,shift,true);
                std::cout << "Print degree matrix of interpolant basis..." << std::endl;
                std::cout << degmat << std::endl;
            }
        }
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
