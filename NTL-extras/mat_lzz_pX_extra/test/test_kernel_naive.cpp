#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

void check(long m, long n, long d){
    Mat<zz_pX> pmat;
    random(pmat, m, n, d+1);

    //cout << "pmat: " << endl << degree_matrix(pmat) << endl;

    // Uniform shift:
    Shift s(m, d);
    // non-uniform shift:
    //Shift s;
    //for (long i = 0; i < m; ++i)
    //    s.emplace_back(i);

    double t1w,t2w;

    Mat<zz_pX> kerbas;
    t1w = GetWallTime();
    kernel_basis_via_approximation(kerbas, pmat, s);
    t2w = GetWallTime();
    cout << "kerbas (appbas): " << endl << degree_matrix(kerbas) << endl;
    cout << "time (appbas): " << t2w-t1w << endl;

    t1w = GetWallTime();
    bool correct = is_kernel_basis(kerbas, pmat, s, ORD_WEAK_POPOV, true, false);
    t2w = GetWallTime();
    std::cout << "Verification: " << (correct ? "correct" : "wrong") << ", time " << (t2w-t1w) << std::endl;

    //kerbas = Mat<zz_pX>();
    //t1w = GetWallTime();
    //kernel_basis_via_interpolation(kerbas, pmat, s);
    //t2w = GetWallTime();
    //
    ////cout << "kerbas (intbas): " << endl << degree_matrix(kerbas) << endl;
    //cout << "time (intbas): " << t2w-t1w << endl;

    //multiply(res, kerbas, pmat);
    //cout << "kerbas zero?: " << IsZero(kerbas) << endl;
    //cout << "product should be zero: " << IsZero(res) << endl;
    //cout << "Dims: " << kerbas.NumRows() << ", " << kerbas.NumCols() << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(4);

    zz_p::init(NTL::GenPrime_long(60));

    if (argc==4)
    {
        check(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    }
    else
        throw std::invalid_argument("Usage: ./test_kernel_naive rdim cdim degree");


    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
