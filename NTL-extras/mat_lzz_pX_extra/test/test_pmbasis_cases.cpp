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

void intbas_test(const long rdim, const long cdim, const long deg, 
        const long order, const DegVec &shift){
    cout << "STARTING INTBAS TEST: " << rdim << " "
         << cdim << " " << deg << " " << order << endl;

    Mat<zz_pX> pmat;
    bool verify = true;
    auto npoints = order;
    random(pmat, rdim, cdim, deg+1);
    zz_p r;
    random(r);
    zz_pX_Multipoint_Geometric eval(r,zz_p(1),npoints);

    // set up pts
    Vec<zz_p> pts;
    zz_pX x;
    SetCoeff(x,1,1);
    eval.evaluate(pts, x); // just gets powers of r

    // set up evaluations of pmat
    Vec<Mat<zz_p>> evals;
    evals.SetLength(npoints);
    for (long d = 0; d < npoints; d++)
        evals[d].SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < pmat.NumRows(); r++)
    {
        for (long c = 0; c < pmat.NumCols(); c++)
        {
            Vec<zz_p> vals;
            eval.evaluate(vals, pmat[r][c]);
            for (long d = 0; d < npoints; d++)
                evals[d][r][c] = vals[d];
        }
    }

    double t1w,t2w;
    std::vector<long> pivdeg;
    std::cout << "~~~Testing pmbasis-geometric (generic, uniform shift)~~~" << std::endl;
    t1w = GetWallTime();
    Mat<zz_pX> intbas;
    pivdeg = pmbasis_geometric(intbas,evals,pts,r,shift);
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

int main(){
    //zz_p::init(NTL::GenPrime_long(16));
    zz_p::init(9001);

    /** INTBAS TESTS ****************************/
    {
        // test non-uniform shifts
        long rdim = 10;
        long cdim = 9;
        long deg = 11;
        long order = 48;
        long t[10] = {11,12,13,14,15,16,17,18,19,20};
        DegVec shift;
        for (long i = 0; i < 10; i++)
            shift.emplace_back(t[i]);
        intbas_test(rdim,cdim,deg,order,shift);
    }
}














