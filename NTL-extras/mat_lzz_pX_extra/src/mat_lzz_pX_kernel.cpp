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

NTL_CLIENT

// TODO: issues when m <= n !! go directly to divide and conquer
// and if we know the first approximant basis will not give any kernel vector, should we also directly divide and conquer?
// TODO: doc mentioning requirement: entries of shift should bound row degrees of pmat
DegVec kernel_basis(
                    Mat<zz_pX> & kerbas,
                    const Mat<zz_pX> & pmat,
                    const Shift & shift
                   )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // find parameter: sum of the m-n largest entries of shift
    // TODO assumes m<n ?
    Shift sorted_shift(shift);
    std::sort(sorted_shift.begin(), sorted_shift.end());
    long rho = 0;
    for (long i = m-n; i < m; i++)
        rho += sorted_shift[i];

    // order for call to approximation
    // TODO threshold ( 3* ?) to determine
    long order = 3 * ceil( (double)rho / n);

    // compute approximant basis
    Mat<zz_pX> appbas;
    DegVec rdeg = pmbasis(appbas, pmat, order, shift);

    // rdeg is now the shift-pivot degree of appbas; deduce shift-row degree
    // which is the componentwise addition pivot degree + shift
    std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::plus<long>());

    // identify submatrix of some rows of appbas which are in the kernel
    // (not necessarily all of them; but in most cases yes)
    // note the criterion: since shift >= rdeg(pmat), we have
    // rdeg >= rdeg(appbas*pmat) and therefore rows with rdeg[i] < order
    // are such that appbas[i] * pmat = 0.
    std::vector<long> ker_rows;
    std::vector<long> other_rows;
    for (long i=0; i<m; ++i)
        if (rdeg[i] < order)
            ker_rows.emplace_back(i);
        else
            other_rows.emplace_back(i);
    long m1 = ker_rows.size();
    long m2 = other_rows.size();
    Mat<zz_pX> P2;
    P2.SetDims(m2, m);
    DegVec rdegP1(m1), rdegP2(m2);

    // copy into P2
    for (long i = 0; i < m2; ++i)
        rdegP2[i] = rdeg[other_rows[i]];
    for (long i = 0; i < m2; ++i)
        P2[i] = appbas[other_rows[i]]; // FIXME could use swap or something, since appbas will be destroyed?

    std::cout << "degrees P" << std::endl << degree_matrix(appbas) << std::endl;
    std::cout << "degrees P2" << std::endl << degree_matrix(P2) << std::endl;

    // copy into P1
    for (long i = 0; i < m1; ++i)
        rdegP1[i] = rdeg[ker_rows[i]];

    // TODO make sure we don't bother about copying P2 in this case
    if (n == 1)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; i++)
            kerbas[i] = appbas[ker_rows[i]];  // (FIXME cf above could use swap?)
        return rdegP1;
    }

    // set up the recursive calls
    for (long i = 0; i < m2; ++i)
        rdegP2[i] -= order; // set rdegP2 = t from paper
    Mat<zz_pX> G;
    multiply(G,P2,pmat);
    RightShift(G,G,order);
    // TODO middleproduct

    // split G
    Mat<zz_pX> G1,G2;
    long n1 = ceil(n/2);
    long n2 = n-n1;
    G1.SetDims(m2, n1);
    G2.SetDims(m2, n2);
    for (long r = 0; r < m2; ++r)
    {
        long c = 0;
        for (; c < n1; ++c)
            G1[r][c] = G[r][c];
        for (; c < n2; ++c)
            G2[r][c] = G[r][n1+c];
    }

    // recursive calls
    Mat<zz_pX> N1, N2;
    DegVec u = kernel_basis(N1, G1, rdegP2);
    multiply(G2, N1, G2);
    DegVec v = kernel_basis(N2, G2, u);
    rdegP1.reserve(m1+v.size());
    for (auto &i: v)
        rdegP1.emplace_back(i);

    // collect output
    multiply(G1,N2,N1);
    if (G1.NumRows() != 0) // FIXME can it be zero??
        multiply(G1,G1,P2);

    kerbas.SetDims(m1+G1.NumRows(), m);
    for (long i = 0; i < m1; ++i)
        kerbas[i] = appbas[ker_rows[i]];  // (FIXME cf above could use swap?)
    for (long i = 0; i < G1.NumRows(); ++i)
        kerbas[m1+i] = G1[i];
    return rdegP1;
}

DegVec kernel_basis_intbas(
                           Mat<zz_pX> & kerbas,
                           const Mat<zz_pX> & pmat,
                           const Shift & shift
                          )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    if (pmat.NumRows() == 0)
    {
        kerbas = Mat<zz_pX>();
        return DegVec();
    }

    long rho = 0;
    for (long i = m-n; i < m; i++)
        rho += shift[i];
    long lambda = ceil((rho*1.0)/n);
    long order = 3*lambda;

    //cout << "pmat: " << degree_matrix(pmat) << endl;
    //cout << "shift: ";
    //for (auto i : shift) cout << i << " ";
    //cout << endl;
    //cout << "rho: " << rho << endl;
    //cout << "order: " << order << endl;

    if (order == 0)
    {
        kerbas = Mat<zz_pX>();
        kerbas.SetDims(m-n,m);
        auto res = DegVec();
        res.resize(m-n);
        return res;
    }

    zz_p r;
    random(r);

    Mat<zz_pX> P;
    Vec<zz_p> pts;
    auto dvec = pmbasis_geometric(P, pmat, r, order, shift,pts);

    // find row degrees
    DegVec rdegP;
    rdegP.resize(m);
    //cout << "P: " << degree_matrix(P) << endl;
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

    //cout << "P1: " << degree_matrix(P1) << endl;
    //cout << "P2: " << degree_matrix(P2) << endl;

    if (n == 1)
    {
        kerbas = P1;
        return rdegP1;
    }

    // set up the recursive calls
    for (unsigned long i = 0; i < rdegP2.size(); i++)
        rdegP2[i] -= order; // set rdegP2 = t from paper
    Mat<zz_pX> G;
    if (P2.NumRows() != 0)
        multiply(G,P2,pmat);

    // divide by the polynomial we are working over
    zz_pX poly, x;
    SetCoeff(x,1,1);
    SetCoeff(poly,0,1);
    zz_p pow_r(1);
    for (long i = 0; i < order; i++)
    {
        if (pow_r != pts[i]) cout << "WRONG POINT" << endl;
        poly *= x-pow_r;
        pow_r *= r*r;
    }
    for (long r = 0; r < G.NumRows(); r++)
        for (long c = 0; c < G.NumCols(); c++)
            divide(G[r][c], G[r][c], poly);

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

    DegVec u = kernel_basis_intbas(N1, G1, rdegP2);
    //cout << "N1: " << N1 << endl;
    multiply(G2, N1, G2);

    DegVec v = kernel_basis_intbas(N2, G2, u);

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
