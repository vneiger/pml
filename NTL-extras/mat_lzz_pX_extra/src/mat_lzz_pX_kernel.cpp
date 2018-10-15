#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_partial_linearization.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

// TODO testing generation
// TODO row_wise
// FIXME randomized: probability of not detecting that it is not correct?
bool is_kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const Shift & shift,
                     const PolMatForm & form,
                     const bool row_wise,
                     const bool randomized
                    )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    if (not row_wise)
        throw std::invalid_argument("~~is_kernel_basis~~ column-wise not implemented yet");

    std::cout << "~~is_kernel_basis~~ Warning: checking generation not implemented yet" << std::endl;

    // test whether kerbas has the right dimensions
    if ( (row_wise && kerbas.NumCols() != m)
        || ((not row_wise) && kerbas.NumRows() != n))
        return false;

    // test whether appbas is shift-reduced with form at least 'form'
    if (not is_polmatform(kerbas,form,shift,row_wise))
        return false;

    // verify that the product is zero
    if (not randomized)
    {
        Mat<zz_pX> product;
        if (row_wise)
            multiply(product, kerbas, pmat);
        else
            multiply(product, pmat, kerbas);
        if (not IsZero(product))
            return false;
    }
    else // randomized
    {
        // use left and right random constant projections to speeds up computations
        // --> "projected product" is likely nonzero if product is nonzero
        Mat<zz_p> left_project, right_project;
        Mat<zz_pX> projected_kerbas, projected_pmat;
        zz_pX projected_product;
        if (row_wise)
        {
            random(left_project, 1, kerbas.NumRows());
            mul(projected_kerbas, left_project, kerbas);
            random(right_project, pmat.NumCols(), 1);
            mul(projected_pmat, pmat, right_project);
            for (long i = 0; i < m; ++i)
                projected_product += projected_kerbas[0][i] * projected_pmat[i][0];
        }
        else
        {
            random(left_project, 1, pmat.NumRows());
            mul(projected_pmat, left_project, pmat);
            random(right_project, kerbas.NumCols(), 1);
            mul(projected_kerbas, kerbas, right_project);
            for (long i = 0; i < m; ++i)
                projected_product += projected_pmat[0][i] * projected_kerbas[i][0];
        }
        if (not IsZero(projected_product))
            return false;
    }

    // testing generation in generic case
    // --> compare sum of pivot degrees for input and output
    long ker_dim = (row_wise ? kerbas.NumRows() : kerbas.NumCols());
    std::vector<long> pivind(ker_dim);
    DegVec pivdeg(ker_dim);
    pivot_index(pivind, pivdeg, kerbas, shift, row_wise);
    long kerbas_degdet = std::accumulate(pivdeg.begin(), pivdeg.end(), 0);
    DegVec degs = vector_degree(pmat, Shift(row_wise?m:n,0), not row_wise);
    long input_degdet = std::accumulate(degs.begin(), degs.end(), 0);
    std::cout << "~~is_kernel_basis~~ Generation check assuming generic input of given ";
    std::cout << (row_wise ? "column" : "row") << " degree: ";
    std::cout << ((kerbas_degdet == input_degdet) ? "correct" : "wrong") << std::endl;

    return true;
}


DegVec kernel_basis(
                    Mat<zz_pX> & kerbas,
                    const Mat<zz_pX> & pmat,
                    const Shift & shift
                   )
{
    return kernel_basis_via_approximation(kerbas, pmat, shift);
}


// TODO structure of output: return pivind+pivdeg?
// TODO avoid computing deg(pmat) all the time... give it as input?
DegVec kernel_basis_via_approximation(
                                      Mat<zz_pX> & kerbas,
                                      const Mat<zz_pX> & pmat,
                                      const Shift & shift
                                     )
{
    // parameters
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    const long d = deg(pmat);

    // compute amplitude of the shift
    long amp = amplitude(shift);

    // compute the order for approximation:
    const long order = (n+1)*d + amp + 1;

    // FIXME: improvement: better performance if pmat does not have
    // balanced column degree would be to use a column-degree wise order
    // (however, this is not handled by fast approximant algorithms for now)
    // Warning: code below not up-to-date: order is wrong for non-uniform shifts.
    //Order order(n);
    //column_degree(order, pmat);
    //long sum_cdeg = std::accumulate(order.begin(), order.end(), (long)0);
    //std::transform(order.begin(), order.end(), order.begin(), [&](long ord){return ord+sum_cdeg+1;});

    // compute approximant basis
    Mat<zz_pX> appbas;
    DegVec pivdeg;
    pivdeg = pmbasis(appbas, pmat, order, shift);
    //pivdeg = approximant_basis(appbas, pmat, order, shift);

    // find rows which belong to the kernel
    std::vector<long> pivot_index;
    std::vector<long> pivot_degree;
    for (long i = 0; i < m; ++i)
        if (pivdeg[i]+amp < order-d)
        {
            pivot_index.push_back(i);
            pivot_degree.push_back(pivdeg[i]);
        }
    long ker_dim = pivot_index.size();

    // move these rows to output basis
    kerbas.SetDims(ker_dim, m);
    for (long i = 0; i < ker_dim; ++i)
        kerbas[i].swap(appbas[pivot_index[i]]);

    return pivot_degree;
}

// TODO: use pivdeg/pivind instead of rdeg?
// TODO: currently only supports n < m (for m <= n: directly to divide and conquer on columns)
// and if we know the first approximant basis will not give any kernel vector, should we also directly divide and conquer?
// TODO: doc mentioning requirement: entries of shift should bound row degrees of pmat
// TODO: why is it currently required that shift STRICTLY bounds degrees? (otherwise crashes)
DegVec kernel_basis_zls_via_approximation(
                                          Mat<zz_pX> & kerbas,
                                          const Mat<zz_pX> & pmat,
                                          const Shift & shift
                                         )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    //std::cout << "call: " << m << "," << n << "," << deg(pmat) << std::endl;

    // FIXME assumes m>n
    // find parameter: sum of the n largest entries of shift
    Shift sorted_shift(shift);
    std::sort(sorted_shift.begin(), sorted_shift.end());
    long rho = 0;
    for (long i = m-n; i < m; i++)
        rho += sorted_shift[i];

    // order for call to approximation
    // TODO threshold (factor 2) to be determined
    // currently, choice is such that the approximant basis should capture the
    // whole kernel when n <= m/2 and pmat is generic (FIXME including for
    // strange shifts/degrees?)
    long order = 2 * ceil( (double)rho / n ) + 1;

    // compute approximant basis
    Mat<zz_pX> appbas;
    DegVec rdeg = pmbasis(appbas, pmat, order, shift);

    // rdeg is now the shift-pivot degree of appbas; deduce shift-row degree
    // which is the componentwise addition pivot degree + shift
    std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::plus<long>());

    // identify submatrix of some rows of appbas which are in the kernel
    // (not necessarily all of them; but in most cases yes)
    // note the criterion: since rdeg(pmat) <= shift, we have
    // rdeg(appbas*pmat) <= rdeg and therefore rows with rdeg[i] < order
    // are such that appbas[i] * pmat = 0.
    std::vector<long> ker_rows;
    std::vector<long> other_rows;
    for (long i=0; i<m; ++i)
        if (rdeg[i] < order)
            ker_rows.emplace_back(i);
        else
            other_rows.emplace_back(i);
    long m1 = ker_rows.size();

    // shifted row degree for the kernel part found by the previous call
    DegVec rdeg1(m1);
    for (long i = 0; i < m1; ++i)
        rdeg1[i] = rdeg[ker_rows[i]];

    // if there was just one column or if the kernel is full (matrix was zero),
    // then return
    if (n == 1 || m1 == m )
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; i++)
            kerbas[i].swap(appbas[ker_rows[i]]);
        return rdeg1;
    }
    // FIXME it is annoying that even if full basis was found, we still do
    // multiplications below... the problem being that we cannot assume that
    // the kernel has dimension m-n. To circumvent this, could test whether sum
    // of pivot degree is the expected one (?).

    // dimension and shifted row degree (order removed) of the non-kernel part
    // of the approximant basis
    long m2 = other_rows.size();
    DegVec rdeg2(m2);
    for (long i = 0; i < m2; ++i)
        rdeg2[i] = rdeg[other_rows[i]] - order;

    // retrieve this non-kernel part of the approximant
    Mat<zz_pX> approx;
    approx.SetDims(m2, m);
    for (long i = 0; i < m2; ++i)
        approx[i].swap(appbas[other_rows[i]]);

    // compute residual
    // TODO use middle product?
    Mat<zz_pX> pmat2;
    multiply(pmat2, approx, pmat);
    RightShift(pmat2, pmat2, order);

    // column-split into two submatrices of column dimension ~ n/2
    Mat<zz_pX> pmat2_sub;
    long n1 = n/2;
    long n2 = n-n1;
    pmat2_sub.SetDims(m2, n1);
    for (long r = 0; r < m2; ++r)
        for (long c = 0; c < n1; ++c)
            pmat2_sub[r][c] = pmat2[r][c];

    // recursive calls
    Mat<zz_pX> kerbas1, kerbas2;
    rdeg2 = kernel_basis_zls_via_approximation(kerbas1, pmat2_sub, rdeg2);

    pmat2_sub.SetDims(m2, n2);
    for (long r = 0; r < m2; ++r)
        for (long c = 0; c < n2; ++c)
            pmat2_sub[r][c] = pmat2[r][n1+c];

    multiply(pmat2, kerbas1, pmat2_sub);
    rdeg2 = kernel_basis_zls_via_approximation(kerbas2, pmat2, rdeg2);

    // if kerbas2 is empty, kerbas1 was already the full kernel; return
    if (kerbas2.NumRows() == 0)
    {
        kerbas.SetDims(m1,m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(appbas[ker_rows[i]]);
        return rdeg1;
    }

    // otherwise, compute the remaining kernel as
    // kerbas2 =  kerbas2 * kerbas1 * approx
    // (using pmat2 as a temp)
    // v1:
    //std::cout << "a*b*c, ";
    //std::cout << kerbas2.NumRows() << "," << kerbas2.NumCols() << "," << deg(kerbas2) << " x ";
    //std::cout << kerbas1.NumRows() << "," << kerbas1.NumCols() << "," << deg(kerbas1) << " x ";
    //std::cout << approx.NumRows() << "," << approx.NumCols() << "," << deg(approx) << std::endl;
    multiply(pmat2, kerbas2, kerbas1);
    multiply(kerbas2, pmat2, approx);
    // v2:
    //multiply(pmat2, kerbas1, approx);
    //multiply(kerbas2, kerbas2, pmat2);

    // update the shifted row degree
    rdeg1.reserve(m1+rdeg2.size());
    for (auto &i: rdeg2)
        rdeg1.emplace_back(i);

    // copy the kernel basis, joining the kernel rows of the approximant basis
    // and the newly found rows of kerbas2
    kerbas.SetDims(m1+kerbas2.NumRows(), m);
    for (long i = 0; i < m1; ++i)
        kerbas[i].swap(appbas[ker_rows[i]]);
    for (long i = 0; i < kerbas2.NumRows(); ++i)
        kerbas[m1+i].swap(kerbas2[i]);

    return rdeg1;
}

DegVec kernel_basis_zls_via_interpolation(
                                          Mat<zz_pX> & kerbas,
                                          const Mat<zz_pX> & pmat,
                                          const Shift & shift
                                         )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    if (pmat.NumRows() == pmat.NumCols())
    {
        kerbas = Mat<zz_pX>();
        return DegVec();
    }

    // find parameter: sum of the m-n largest entries of shift
    // TODO assumes m<n ?
    // all below assumes pmat nonzero?
    Shift sorted_shift(shift);
    std::sort(sorted_shift.begin(), sorted_shift.end());
    long rho = 0;
    for (long i = m-n; i < m; i++)
        rho += sorted_shift[i]+1;

    // order for call to approximation
    // TODO threshold ( 3* ?) to determine
    long order = 3 * ceil( (double)rho / n);

    cout << "order: " << order << endl;
    cout << "shift: ";
    for (auto i : shift)
        cout << i << " ";
    cout << endl;

    zz_p r;
    random(r);

    Mat<zz_pX> P;
    Vec<zz_p> pts;
    Vec<Mat<zz_p>> evals;
    auto dvec = pmbasis_geometric(P,pmat,r,order,shift,evals,pts);

    // find row degrees
    DegVec rdegP;
    rdegP.resize(m);
    row_degree(rdegP,P,shift);
    cout << "rdegP: ";
    for (auto i: rdegP)
        cout << i << " ";
    cout << endl;

    // partition
    Mat<zz_pX> &P1 = kerbas;
    Mat<zz_pX> P2;
    DegVec rdegP1, rdegP2;

    long m1 = 0;
    for (auto &i : rdegP)
        if (i < order) m1++;
    P1.SetDims(m1, P.NumCols());
    rdegP1.resize(m1);
    P2.SetDims(m-m1, P.NumCols());
    rdegP2.resize(m-m1);

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

    cout << "pmat: " << degree_matrix(pmat) << endl;
    cout << "P1: " << degree_matrix(P1) << endl;
    cout << "P2: " << degree_matrix(P2) << endl;
    if (P1.NumRows() != 0)
    {
        Mat<zz_pX> tmp;
        multiply(tmp,P1,pmat);
        cout << "P1 prod: " << degree_matrix(tmp) << endl;
    }

    if (n == 1 || r1 == m)
    {
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
        poly *= x-pts[i];
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

    cout << "\n\ncall 1" << endl;
    DegVec u = kernel_basis_zls_via_interpolation(N1, G1, rdegP2);
    cout << "u: ";
    for (auto i : u)
        cout << i << " ";
    cout << endl;

    multiply(G2, N1, G2);

    cout << "\n\ncall 2" << endl;
    DegVec v = kernel_basis_zls_via_interpolation(N2, G2, u);
    cout << "v: ";
    for (auto i : v)
        cout << i << " ";
    cout << endl;

    if (N2.NumRows() == 0)
    {
        return rdegP1;
    }

    // collect output
    multiply(G1,N2,N1);
    multiply(G1,G1,P2);
    kerbas.SetDims(P1.NumRows()+G1.NumRows(), P1.NumCols());
    for (long r = 0; r < G1.NumRows(); r++)
    {
        kerbas[r+m1] = G1[r];
    }
    for (auto &i: v)
        rdegP1.emplace_back(i);
    cout << "return degree: ";
    for (auto i : rdegP1)
        cout << i << " ";
    cout << endl;
    return rdegP1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
