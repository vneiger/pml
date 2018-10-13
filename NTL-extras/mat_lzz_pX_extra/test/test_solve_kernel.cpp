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

void check(long m, long n, long d1, long d2)
{
    Mat<zz_pX> pmat;
    random(pmat,m,n,d1+1);
    cout << "pmat: " << endl << degree_matrix(pmat) << endl;
    
    Mat<zz_pX> bmat;
    random(bmat,1,n,d2+1);
    Vec<zz_pX> b = bmat[0];
    
    cout << "b:\n " << degree_matrix(bmat) << endl;
    
    Vec<zz_pX> a;
    zz_pX denom;
    linsolve_via_kernel(a,denom,pmat,b);
		// should compute a and denom such that a * pmat = denom * b
    
    Mat<zz_pX> amat;
    amat.SetDims(1,n);
    amat[0] = a;
    
		// test equality: res = a * pmat, is it denom * b?
    Mat<zz_pX> res;
    multiply(res,amat,pmat);
    for (long i = 0; i < n; i++)
        b[i] = b[i] * denom;
    
    cout << "Verification: " << (res[0] == b ? "correct" : "wrong") << endl;
}

int main(int argc, char ** argv)
{
    SetNumThreads(4);

    //zz_p::init(NTL::GenPrime_long(60));
    zz_p::init(NTL::GenPrime_long(60));

    if (argc==5)
    {
        check(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),atoi(argv[4]));
    }else
    {
        cout << "enter rows, cols, deg (pmat), deg (b)" << endl;
        return 1;
    }


    return 0;
}
