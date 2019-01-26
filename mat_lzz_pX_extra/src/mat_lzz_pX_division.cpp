#include "mat_lzz_pX_extra.h"

NTL_CLIENT

// Left division: 
// B must be row reduced
// computes Q,R such that A = BQ + R
// Follows classical algorithm: reverse and truncated inverse
// (e.g. detailed in Algorithm 1 in [Neiger-Vu ISSAC 17])
void quo_rem(
             Mat<zz_pX> & Q,
             Mat<zz_pX> & R,
             const Mat<zz_pX> & A,
             const Mat<zz_pX> & B
            )
{
    const long m = B.NumRows();
    const long n = A.NumCols();

    // step 0: find parameter d = max(rdeg(A) - rdeg(B))
    // (corresponds to delta-1 in [Neiger-Vu ISSAC 17]
    VecLong rdegA;
    row_degree(rdegA, A);
    VecLong rdegB;
    row_degree(rdegB, B);
    long d = rdegA[0] - rdegB[0];
    for (long i = 1; i < m; ++i)
        if (rdegA[i] - rdegB[i] > d)
            d = rdegA[i] - rdegB[i];

    // if rdeg(A) < rdeg(B), then A is already reduced
    // modulo B: quotient is 0, remainder is A
    if (d < 0)
    {
        R = A;
        Q.SetDims(m,n);
        clear(Q);
        return;
    }

    // step 1: reverse input matrices
    // B is reversed w.r.t. its row degree rdegB
    Mat<zz_pX> Brev;
    row_reverse(Brev, B, rdegB);
    // A is reversed w.r.t. d + rdegB
    // we directly modify rdegB += d, since it is not used afterwards
    std::for_each(rdegB.begin(), rdegB.end(), [&](long & r){r+=d;});
    Mat<zz_pX> buf;
    row_reverse(buf, A, rdegB);

    // step 2: compute quotient
    // Qrev = Brev^{-1} R mod X^{d+1}
    // (we use R to store the reversed quotient Qrev)
    solve_series(R, Brev, buf, d+1);
    reverse(Q, R, d);

    // step 3: deduce remainder
    // R = A - B*Q
    multiply(buf, B, Q);
    sub(R, A, buf);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
