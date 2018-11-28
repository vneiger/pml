#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <NTL/BasicThreadPool.h>
#include <cmath>

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_sequence.h"
#include "structured_lzz_p.h"

NTL_CLIENT

int main(int argc, char *argv[])
{
    if (argc != 3) throw std::invalid_argument("usage: ./time_compare m d");
    zz_p::FFTInit(0);

    long m = atoi(argv[1]);
    long d = atoi(argv[2]);

    // pmbasis
    Mat<zz_pX> pmat;
    random(pmat,2*m, m, d);
    Mat<zz_pX> appbas;
    VecLong shift(2*m,0);

    double t = GetWallTime();
    pmbasis(appbas, pmat, d+1, shift);
    cout << "time pmbasis: " << GetWallTime() - t << endl;

    // structured matrix
    long len = 2*d-1; //num of elements needed dxd blocks
    Vec<Vec<hankel_lzz_p>> h_vec; //set to m blocks by m blocks
    h_vec.SetLength(m);
    for (long i = 0; i < m ; i++)
    {
        h_vec[i].SetLength(m);
    }
    
    for (long i = 0; i < h_vec.length(); i++)
    {
        Vec<zz_p> rvec;
        rvec.SetLength(len);
        for (long j = 0; j < h_vec[0].length(); j++)
        {
            for (long t = 0; t < len; t++)
                random(rvec[t]);
            h_vec[i][j] = hankel_lzz_p(rvec,d,d);
        }
    }
    mosaic_hankel_lzz_p mh{h_vec};
    Vec<zz_p> b;
    b.SetLength(m*d);
    for (long i = 0; i < m*d; i++)
        random(b[i]);

    Vec<zz_p> x;
    t = GetWallTime();
    mh.solve(x,b);
    cout << "time mh solve: " << GetWallTime() - t << endl;
    
}

