#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_partial_linearization.h"
#include "lzz_pX_CRT.h"

DegVec kernel_basis(
                    Mat<zz_pX> & kerbas,
                    const Mat<zz_pX> & pmat,
                    const Shift & shift
                   )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    
    long rho = 0;
    for (long i = m-n; i < m; i++)
        rho += shift[i];
    long lambda = ceil((rho*1.0)/n);
    long order = 3*lambda;
    
    cout << "rho: " << rho << endl;
    cout << "order: " << order << endl;
    
    Mat<zz_pX> P;
    auto dvec = pmbasis(P, pmat, order, shift);
    
    // find row degrees
    DegVec rdegP;
    rdegP.resize(m);
    cout << "P: " << degree_matrix(P) << endl;
    row_degree(rdegP,P,shift);
    
    // partition
    Mat<zz_pX> P1,P2;
    DegVec rdegP1, rdegP2;
    
    long row_P1 = 0;
    for (auto &i : rdegP)
        if (i < order) row_P1++;
    P1.SetDims(row_P1, P.NumCols());
    rdegP1.resize(row_P1);
    P2.SetDims(m-row_P1, P.NumCols());
    rdegP2.resize(m-row_P1);
    
    // copy into P1,P2
    long r1 = 0;
    long r2 = 0;
    for (long r = 0; r < m; r++)
    {
        if (rdegP[r] < order) // copy into P1
        {
            rdegP1[r1] = rdegP[r];
            for (long c = 0; c < P.NumCols(); c++)
                P1[r1][c] = P[r][c];
            r1++;
        }else // copy into P2
        {
            rdegP2[r2] = rdegP[r];
            for (long c = 0; c < P.NumCols(); c++)
                P2[r2][c] = P[r][c];
            r2++;
        }
    }
    
    cout << "P1: " << degree_matrix(P1) << endl;
    cout << "P2: " << degree_matrix(P2) << endl;
    
    if (n == 1)
    {
        kerbas = P1;
        return rdegP1;
    }
    
    // set up the recursive calls
    for (unsigned long i = 0; i < rdegP2.size(); i++)
        rdegP2[i] -= order; // set rdegP2 = t from paper
    Mat<zz_pX> G;
    multiply(G,P2,pmat);
    RightShift(G,G,order);
    
    // split G
    Mat<zz_pX> G1,G2;
    long n1 = ceil(n/2);
    long n2 = n-n1;
    G1.SetDims(P2.NumRows(), n1);
    G2.SetDims(P2.NumRows(), n2);
    for (long r = 0; r < P2.NumRows(); r++)
    {
        long c_at = 0;
        for (long c = 0; c < n1; c++, c_at++)
            G1[r][c] = G[r][c_at];
        for (long c = 0; c < n2; c++, c_at++)
            G2[r][c] = G[r][c_at];
    }
    
    // recursive calls
    Mat<zz_pX> N1, N2;
    DegVec u = kernel_basis(N1, G1, rdegP2);
    multiply(G2, N1, G2);
    DegVec v = kernel_basis(N2, G2, u);
    
    // collect output
    multiply(G1,N2,N1);
    if (G1.NumRows() != 0)
        multiply(G1,G1,P2);
    kerbas.SetDims(P1.NumRows()+G1.NumRows(), P1.NumCols());
    long r_at = 0;
    for (long r = 0; r < P1.NumRows(); r++, r_at++)
    {
        for (long c = 0; c < kerbas.NumCols(); c++)
        {
            kerbas[r_at][c] = P1[r][c];
        }
    } 
    for (long r = 0; r < G1.NumRows(); r++, r_at++)
    {
        for (long c = 0; c < kerbas.NumCols(); c++)
        {
            kerbas[r_at][c] = G1[r][c];
        }
    }
    for (auto &i: v)
        rdegP1.emplace_back(i);
    return rdegP1;
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

















