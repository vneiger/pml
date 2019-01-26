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
    if (m<25)
        cout << "pmat: " << endl << degree_matrix(pmat) << endl;

    // random rhs
    Mat<zz_pX> bmat;
    random(bmat,1,m,d2+1);
    Vec<zz_pX> b = bmat[0];
     std::cout << "random right-hand side of degree " << d2 << std::endl;

    // rhs in row space, constant combination
    //Vec<zz_pX> b; b.SetLength(m);
    //for (long i = 0; i < m; ++i)
    //{
    //    zz_p a; random(a);
    //    for (long j = 0; j < m; ++j)
    //        b[j] += a*pmat[i][j];
    //}
    //std::cout << "Right-hand side: random constant linear combination of system matrix" << std::endl;

    // rhs in row space, degree d2/2 linear combination
    //Vec<zz_pX> b; b.SetLength(m);
    //for (long i = 0; i < m; ++i)
    //{
    //    zz_pX a; random(a, d2/2+1);
    //    for (long j = 0; j < m; ++j)
    //        b[j] += a*pmat[i][j];
    //}
    //std::cout << "Right-hand side: random linear combination of system matrix, cofactors of degree " << (d2/2) << std::endl;

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

    if (m<20)
    {
        std::cout << "Degrees in solution vector\n" << degree_matrix(amat) << std::endl;
        std::cout << "Degree denominator\n" << deg(denom) << std::endl;
    }

    // test equality: res = a * pmat, is it denom * b?
    Mat<zz_pX> res;
    multiply(res,amat,pmat);
    for (long i = 0; i < m; i++)
        b[i] = b[i] * denom;

    cout << "Product verification: " << (res[0] == b ? "correct" : "wrong") << endl;

    zz_pX det;
    t1 = GetWallTime();
    determinant_generic_knowing_degree(det, pmat, m*d1);
    t2 = GetWallTime();
    std::cout << "Time (determinant): " << (t2-t1) << std::endl;
    MakeMonic(det);
    std::cout << "Denom is determinant? "  << (det == denom ? "yes" : "no") << std::endl;
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
