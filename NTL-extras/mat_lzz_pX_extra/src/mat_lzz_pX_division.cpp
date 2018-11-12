#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

void quo_rem(Mat<zz_pX> &Q, 
             Mat<zz_pX> &R, 
             const Mat<zz_pX> &A,
             const Mat<zz_pX> &B)
{
    const long m = B.NumRows();
    const long n = A.NumCols();

    // step 0: find parameter d (delta in Neiger-Vu17)
    DegVec rdegA;
    DegVec rdeg;
    rdegA.resize(m);
    rdeg.resize(m);
    row_degree(rdegA, A);
    row_degree(rdeg, B);
    long d = rdegA[0] - rdeg[0] + 1;
    for (long i = 1; i < m; i++)
    {
        if (rdegA[i] - rdeg[i] + 1 > d)
            d = rdegA[i] - rdeg[i] + 1;
    }
    if (d <= 0) // A already reduce mod B, Q = 0
    {
        R = A;
        Q = Mat<zz_pX>();
        Q.SetDims(m,n);
    }

    // step 1: reverse input matrices
    Mat<zz_pX> Arev, Brev;
    DegVec Arevdeg;
    for (long i = 0; i < m; i++)
        Arevdeg.emplace_back(d+rdeg[i]-1);
    reverse(Arev, A, Arevdeg);
    reverse(Brev, B, rdeg);
    
    // step 2: compute quotient
    solve_series_high_precision(Q, Brev, Arev, d);
    reverse(Q,Q,d-1);
    
    // step 3: deduce remainder
    Mat<zz_pX> T;
    multiply(T, B, Q);
    R = A - T;
}

























