#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <NTL/BasicThreadPool.h>
#include <cmath>

//#define SAFETY_CHECKS

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"

NTL_CLIENT

int main(int argc, char *argv[])
{   
    if (argc!=7)
        throw std::invalid_argument("Usage: ./test_appbas_mbasis rdim cdim deg nbits verify nthreads");

    long rdim = atoi(argv[1]);
    long cdim = atoi(argv[2]);
    long deg = atoi(argv[3]);
    long nbits = atoi(argv[4]);
    bool verify = (atoi(argv[5])==1);
    SetNumThreads(atoi(argv[6]));

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));
    
    Mat<zz_pX> A;
    Mat<zz_pX> B;    
    random(A,rdim,cdim,deg);
    random(B,rdim,rdim,deg);
    
    Mat<zz_pX> Q;
    Mat<zz_pX> R;
    
    quo_rem(Q,R,A,B);
    
    if (verify)
    {
        Mat<zz_pX> T;
        multiply(T,B,Q);
        T = T+R;
        if (T == A) cout << "correct" << endl;
        else cout << "wrong" << endl;
    }
}












