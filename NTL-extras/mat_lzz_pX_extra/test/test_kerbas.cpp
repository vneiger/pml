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
    random_mat_zz_pX(pmat, m, n, d);
    
    cout << "pmat: " << degree_matrix(pmat) << endl;
    
    Shift s;
    for (long i = 0; i < m; i++)
        s.emplace_back(d);
    
    Mat<zz_pX> kerbas;
    kernel_basis(kerbas, pmat, s);
    
    cout << "kerbas (appbas): " << degree_matrix(kerbas) << endl;
    
    Mat<zz_pX> res;
    multiply(res, kerbas, pmat);
    cout << degree_matrix(res) << endl;
    
    kerbas = Mat<zz_pX>();
    kernel_basis_intbas(kerbas, pmat, s);
    
    cout << "kerbas (intbas): " << degree_matrix(kerbas) << endl;
    
    multiply(res, kerbas, pmat);
    cout << degree_matrix(res) << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    //SetNumThreads(4);
    
    zz_p::init(NTL::GenPrime_long(60));

    if (argc==4)
    {
        check(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    }
    

    return 0;
}
