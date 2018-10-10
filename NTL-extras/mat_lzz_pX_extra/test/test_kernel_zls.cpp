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
    random_mat_zz_pX(pmat, m, n, d+1);

    //cout << "pmat: " << endl << degree_matrix(pmat) << endl;

    // Uniform shift:
    // Shift s(m, d);
    // uniform shift:
    Shift s;
    for (long i = 0; i < m; ++i)
        s.emplace_back(d+i);
    double t1w,t2w;

    Mat<zz_pX> kerbas;
    t1w = GetWallTime();
    kernel_basis(kerbas, pmat, s);
    t2w = GetWallTime();
    //cout << "kerbas (appbas): " << endl << degree_matrix(kerbas) << endl;
    cout << "time (appbas): " << t2w-t1w << endl;

    Mat<zz_pX> res;
    multiply(res, kerbas, pmat);
    cout << "Dims: " << kerbas.NumRows() << ", " << kerbas.NumCols() << endl;
    cout << "kerbas shouldn't be zero: " << IsZero(kerbas) << endl;
    cout << "product should be zero: " << boolalpha << IsZero(res) << endl << endl;

    kerbas = Mat<zz_pX>();
    t1w = GetWallTime();
    kernel_basis_intbas(kerbas, pmat, s);
    t2w = GetWallTime();
    
    //cout << "kerbas (intbas): " << endl << degree_matrix(kerbas) << endl;
    cout << "time (intbas): " << t2w-t1w << endl;

    multiply(res, kerbas, pmat);
    cout << "kerbas shouldn't be zero: " << IsZero(kerbas) << endl;
    cout << "product should be zero: " << IsZero(res) << endl;
    cout << "Dims: " << kerbas.NumRows() << ", " << kerbas.NumCols() << endl;
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


    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
