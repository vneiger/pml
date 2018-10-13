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

void check(long m, long d1, long d2)
{
    double t1, t2;

    Mat<zz_pX> pmat;
    random(pmat,m,m,d1+1);
    cout << "pmat: " << endl << degree_matrix(pmat) << endl;

    Mat<zz_pX> bmat;
    random(bmat,1,m,d2+1);
    Vec<zz_pX> b = bmat[0];

    cout << "b:\n " << degree_matrix(bmat) << endl;

    t1 = GetWallTime();
    Vec<zz_pX> a;
    zz_pX denom;
    // should compute a and denom such that a * pmat = denom * b
    linsolve_via_kernel(a,denom,pmat,b);
    t2 = GetWallTime();
    std::cout << "Time (linsolve_kernel): " << (t2-t1) << std::endl;

    Mat<zz_pX> amat;
    amat.SetDims(1,m);
    amat[0] = a;

    std::cout << "Degrees in solution vector\n" << degree_matrix(amat) << std::endl;
    std::cout << "Degree denominator\n" << deg(b) << std::endl;

    // test equality: res = a * pmat, is it denom * b?
    Mat<zz_pX> res;
    multiply(res,amat,pmat);
    for (long i = 0; i < m; i++)
        b[i] = b[i] * denom;

    cout << "Verification: " << (res[0] == b ? "correct" : "wrong") << endl;
}

int main(int argc, char ** argv)
{
    SetNumThreads(4);

    zz_p::init(NTL::GenPrime_long(60));

    if (argc==4)
        check(atoi(argv[1]), atoi(argv[2]),atoi(argv[3]));
    else
        throw std::invalid_argument("Usage: ./test_solve_kernel rdim degmat degvec");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
