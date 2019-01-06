#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"
#include "mat_lzz_pX_kernel.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*            DEGREE BOUNDS FOR KERNEL BASES                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// Assuming pmat has full column rank n, an unpublished result by Vu Thi Xuan
// (Master's research report, Lemma 10) shows that the sum of the pivot degrees
// of the shifted Popov kernel basis (for any shift) is at most the degree of
// the determinant of the complement part in 'pmat' (complement: rows not in
// the shifted pivot support of the kernel).  This can probably be derived from
// knowledge about irreducible fractions.
//
// As a result, we have the following non-strict upper bounds on the sum of
// pivot degree of the shifted Popov kernel basis of pmat:
//   * pmat.NumCols() * degree(pmat)
//   * sum of column degrees of pmat
//   * sum of degrees of the rows in the complement of pmat
//          (useful if pivot support is known)
//   * if pivot support is unknown, since the complement has <= n rows, this
//   can be relaxed as: sum of degrees of the n largest-degree rows of pmat.
//
// The max of the pivot degrees, which is also the degree of the pivot part of
// the shifted Popov kernel, can be equal to the sum of pivot degrees (although
// this is not expected generically). Hence we have the same bounds on degree
// of pivot part. This directly gives bounds on the non-pivot part: pick one of
// the bounds above, and add max(shift) - min(shift). Note that this might be
// pessimistic for shifts of large amplitude, since we also have bounds on the
// non-pivot part involving the adjugate of the pivot support complement part
// of pmat, ensuring in particular that the non-pivot part of the kernel has
// degree at most n * deg(pmat)   (cf Lemma 12 of Vu Thi Xuan's research
// report).
//
// The max of pivot degrees of the kernel also gives a bound on
// max(rdeg_s(kernel)): pick one of the bounds above for sum of pivot degrees,
// and add max(s).
//
// Case where pmat does not have full column rank: the left kernel of pmat is
// equal to the left kernel of any column basis of pmat; taking a column
// reduced form of pmat allows us to reduce to the full column rank case. This
// shows that the bounds above still hold without full column rank assumption.




// TODO testing generation
// TODO row_wise
// FIXME randomized: probability of not detecting that it is not correct?
bool is_kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift,
                     const PolMatForm & form,
                     const bool randomized
                    )
{
    const long m = pmat.NumRows();

    std::cout << "~~is_kernel_basis~~ Warning: checking generation not implemented yet" << std::endl;

    // test whether kerbas has the right dimensions
    if (kerbas.NumCols() != m)
        return false;

    // test whether appbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,kerbas,shift))
        return false;

    // verify that the product is zero
    if (not randomized)
    {
        Mat<zz_pX> product;
        multiply(product, kerbas, pmat);
        if (not IsZero(product))
            return false;
    }
    else // randomized
    {
        // Freivalds like:
        // use left and right random constant projections to speeds up computations
        // --> "projected product" is likely nonzero if product is nonzero
        Mat<zz_p> left_project, right_project;
        Mat<zz_pX> projected_kerbas, projected_pmat;
        zz_pX projected_product;
        random(left_project, 1, kerbas.NumRows());
        mul(projected_kerbas, left_project, kerbas);
        random(right_project, pmat.NumCols(), 1);
        mul(projected_pmat, pmat, right_project);
        for (long i = 0; i < m; ++i)
            projected_product += projected_kerbas[0][i] * projected_pmat[i][0];
        if (not IsZero(projected_product))
            return false;
    }

    //// testing generation in generic case
    //// --> compare sum of pivot degrees for input and output
    //long ker_dim = kerbas.NumRows();
    //VecLong pivind(ker_dim);
    //VecLong pivdeg(ker_dim);
    //row_pivots(pivind, pivdeg, kerbas, shift);
    //long kerbas_degdet = std::accumulate(pivdeg.begin(), pivdeg.end(), 0);
    //VecLong degs;
    //row_degree(degs, pmat);
    //long input_degdet = std::accumulate(degs.begin(), degs.end(), 0);
    //std::cout << "~~is_kernel_basis~~ Generation check assuming generic input of given ";
    //std::cout << "column degree: ";
    //std::cout << ((kerbas_degdet == input_degdet) ? "correct" : "wrong") << std::endl;

    return true;
}


VecLong kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift
                    )
{
    return kernel_basis_via_approximation(kerbas, pmat, shift);
}


// TODO improvement: better performance if pmat has very unbalanced column
// degree would be to use a column-degree wise order (however, this is not
// handled by fast approximant algorithms for now)
VecLong kernel_basis_via_approximation(
                                       Mat<zz_pX> & kerbas,
                                       const Mat<zz_pX> & pmat,
                                       const VecLong & shift
                                      )
{
    // parameters
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    const long d = deg(pmat);

    if (d==-1)
    {
        ident(kerbas, m);
        return VecLong(m);
    }

    // compute amplitude of the shift
    const long amp = amplitude(shift);

    // compute the order for approximation:
    const long order = (n+1)*d + amp + 1;

    // compute approximant basis
    Mat<zz_pX> appbas;
    VecLong pivdeg = pmbasis(appbas, pmat, order, shift);

    // find rows which belong to the kernel
    VecLong pivot_index;
    VecLong pivot_degree;
    for (long i = 0; i < m; ++i)
        if (pivdeg[i]+amp < order-d)
        {
            pivot_index.push_back(i);
            pivot_degree.push_back(pivdeg[i]);
        }
    const long ker_dim = pivot_index.size();

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
VecLong kernel_basis_zls_via_approximation(
                                           Mat<zz_pX> & kerbas,
                                           Mat<zz_pX> & pmat,
                                           const VecLong & shift
                                          )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    //std::cout << "call: " << m << "," << n << "," << deg(pmat) << std::endl;

    // FIXME here, assumes m>n
    // find parameter rho: sum of the n largest entries of shift
    VecLong sorted_shift(shift);
    std::sort(sorted_shift.begin(), sorted_shift.end());
    const long rho = std::accumulate(sorted_shift.begin()+m-n, sorted_shift.end(), 0);

    // order for call to approximation
    // TODO threshold (factor 2) to be determined
    // choosing 2 ensures that the approximant basis captures the whole kernel
    // when n <= m/2, pmat is sufficiently generic, and shift is uniform
    const long order = 2 * ceil((double)rho / n) + 1;

    // compute approximant basis
    // --> it does not necessarily capture the whole kernel; but does in two
    // notable cases: if n=1, or in the generic situation mentioned above
    Mat<zz_pX> appbas;
    VecLong rdeg = pmbasis(appbas, pmat, order, shift); // does not change pmat or shift

    // rdeg is now the shift-pivot degree of appbas; deduce shift-row degree
    // which is the componentwise addition pivot degree + shift
    std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::plus<long>());

    // Identify submatrix of some rows of appbas which are in the kernel
    // Note the criterion: since rdeg(pmat) <= shift, we have
    // rdeg(appbas*pmat) <= rdeg and therefore rows with rdeg[i] < order
    // are such that appbas[i] * pmat = 0.
    // Note that this may miss rows in the kernel but with rdeg >= order
    // (this is usually not the case).
    VecLong ker_rows;
    VecLong other_rows;
    for (long i=0; i<m; ++i)
        if (rdeg[i] < order)
            ker_rows.emplace_back(i);
        else
            other_rows.emplace_back(i);
    const long m1 = ker_rows.size();

    // shifted row degree for the kernel part found by the previous call
    VecLong rdeg1(m1);
    for (long i = 0; i < m1; ++i)
        rdeg1[i] = rdeg[ker_rows[i]];

    // if there was just one column or if the kernel is full (matrix was zero),
    // then return
    if (n == 1 || m1 == m )
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(appbas[ker_rows[i]]);
        return rdeg1;
    }
    // FIXME it is annoying that even if full basis was found, we still do
    // multiplications below... the problem being that we cannot assume that
    // the kernel has dimension m-n. To circumvent this, could test whether sum
    // of pivot degree is the expected one (?).

    // dimension and shifted row degree (order removed) of the non-kernel part
    // of the approximant basis
    const long m2 = other_rows.size();
    VecLong rdeg2(m2);
    for (long i = 0; i < m2; ++i)
        rdeg2[i] = rdeg[other_rows[i]] - order;

    // retrieve this non-kernel part of the approximant
    Mat<zz_pX> approx(INIT_SIZE, m2, m);
    for (long i = 0; i < m2; ++i)
        approx[i].swap(appbas[other_rows[i]]);

    // compute residual
    // TODO use middle product?
    multiply(pmat, approx, pmat);
    RightShift(pmat, pmat, order);

    // pmat will be column-splitted into two submatrices of column dimension ~ n/2
    const long n1 = n/2;
    const long n2 = n-n1;

    // recursive call 1, with left submatrix of the residual pmat
    Mat<zz_pX> pmat_sub(INIT_SIZE, m2, n1);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n1; ++j)
            pmat_sub[i][j] = pmat[i][j];

    Mat<zz_pX> kerbas1;
    rdeg2 = kernel_basis_zls_via_approximation(kerbas1, pmat_sub, rdeg2);

    // recursive call 2, with right submatrix of the residual pmat
    pmat_sub.SetDims(m2, n2);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n2; ++j)
            pmat_sub[i][j] = pmat[i][n1+j];
    multiply(pmat, kerbas1, pmat_sub);
    pmat_sub.kill();

    // recursive call 2
    rdeg2 = kernel_basis_zls_via_approximation(kerbas, pmat, rdeg2);

    // if kerbas is non-empty, it corresponds to new kernel rows: complete
    // the computation of these rows via kerbas =  kerbas * kerbas1 * approx
    // (otherwise, the approximant basis already captured the whole kernel)
    if (kerbas.NumRows() != 0)
    {
        // (using pmat as a temp)
        // v1:
        //multiply(pmat, kerbas, kerbas1);
        //multiply(kerbas, pmat, approx);
        // v2:
        multiply(pmat, kerbas1, approx);
        multiply(kerbas, kerbas, pmat);
    }

    // update the shifted row degree
    rdeg2.reserve(m1+rdeg2.size());
    rdeg2.insert(rdeg2.end(), rdeg1.begin(), rdeg1.end());

    // insert the kernel rows of the approximant basis
    kerbas.SetDims(m1+kerbas.NumRows(), m);
    for (long i = 0; i < m1; ++i)
        kerbas[kerbas.NumRows()-m1+i].swap(appbas[ker_rows[i]]);

    return rdeg2;
}

VecLong kernel_basis_zls_via_interpolation(
                                           Mat<zz_pX> & kerbas,
                                           const Mat<zz_pX> & pmat,
                                           const VecLong & shift
                                          )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    if (pmat.NumRows() == pmat.NumCols())
    {
        kerbas = Mat<zz_pX>();
        return VecLong();
    }

    // find parameter: sum of the m-n largest entries of shift
    // TODO assumes m<n ?
    // all below assumes pmat nonzero?
    VecLong sorted_shift(shift);
    std::sort(sorted_shift.begin(), sorted_shift.end());
    long rho = 0;
    for (long i = m-n; i < m; i++)
        rho += sorted_shift[i]+1;

    // order for call to approximation
    // TODO threshold ( 3* ?) to determine
    long order = 3 * ceil( (double)rho / n);

    //cout << "order: " << order << endl;
    //cout << "shift: ";
    //for (auto i : shift)
    //    cout << i << " ";
    //cout << endl;

    zz_p r;
    random(r);

    Mat<zz_pX> P;
    Vec<zz_p> pts;
    auto dvec = pmbasis_geometric(P,pmat,r,order,shift,pts);

    // find row degrees
    VecLong rdegP;
    rdegP.resize(m);
    row_degree(rdegP,P,shift);
    //cout << "rdegP: ";
    //for (auto i: rdegP)
    //    cout << i << " ";
    //cout << endl;

    // partition
    Mat<zz_pX> &P1 = kerbas;
    Mat<zz_pX> P2;
    VecLong rdegP1, rdegP2;

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

    //cout << "pmat: " << degree_matrix(pmat) << endl;
    //cout << "P1: " << degree_matrix(P1) << endl;
    //cout << "P2: " << degree_matrix(P2) << endl;
    if (P1.NumRows() != 0)
    {
        Mat<zz_pX> tmp;
        multiply(tmp,P1,pmat);
        //cout << "P1 prod: " << degree_matrix(tmp) << endl;
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

    //cout << "\n\ncall 1" << endl;
    VecLong u = kernel_basis_zls_via_interpolation(N1, G1, rdegP2);
    //cout << "u: ";
    //for (auto i : u)
    //    cout << i << " ";
    //cout << endl;

    multiply(G2, N1, G2);

    //cout << "\n\ncall 2" << endl;
    VecLong v = kernel_basis_zls_via_interpolation(N2, G2, u);
    //cout << "v: ";
    //for (auto i : v)
    //    cout << i << " ";
    //cout << endl;

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
    //cout << "return degree: ";
    //for (auto i : rdegP1)
    //    cout << i << " ";
    //cout << endl;
    return rdegP1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
