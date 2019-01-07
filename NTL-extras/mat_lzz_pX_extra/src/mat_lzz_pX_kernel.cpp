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
    {
        std::cout << "Wrong dimension" << std::endl;
        return false;
    }

    // test whether appbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,kerbas,shift))
    {
        std::cout << "Wrong polmatform" << std::endl;
        return false;
    }

    // verify that the product is zero
    if (not randomized)
    {
        Mat<zz_pX> product;
        multiply(product, kerbas, pmat);
        if (not IsZero(product))
        {
            std::cout << "Not in kernel" << std::endl;
            return false;
        }
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
        {
            std::cout << "Not in kernel" << std::endl;
            return false;
        }
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



void kernel_basis_zls_via_approximation_new(
                                        Mat<zz_pX> & kerbas,
                                        Mat<zz_pX> & pmat,
                                        VecLong & shift,
                                        VecLong & pivind,
                                        VecLong & pivdeg
                                       )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    //std::cout << m << "\t" << n << "\t" << deg(pmat) << std::endl;

    // if m==1, just check whether pmat is zero
    if (m==1)
    {
        if (not IsZero(pmat))
        {
            kerbas.SetDims(0,1);
            pivind.clear();
            pivdeg.clear();
            return;
        }
        // pmat is the zero 1 x n matrix
        ident(kerbas, 1);
        pivind.resize(1); pivind[0] = 0;
        pivdeg.resize(1); pivdeg[0] = 0;
        return;
    }

    // if n > m/2 (with m>=2), just split the columns into two submatrices or
    // roughly equal dimensions, and call the algorithm recursively
    if (2*n > m)
    {
        // pmat will be column-splitted into two submatrices of column dimension ~ n/2
        const long n1 = n/2;
        const long n2 = n-n1;

        // recursive call 1, with left submatrix of pmat
        Mat<zz_pX> pmat_sub(INIT_SIZE, m, n1);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < n1; ++j)
                pmat_sub[i][j] = pmat[i][j];

        Mat<zz_pX> kerbas1;
        VecLong pivind1, pivdeg1;
        kernel_basis_zls_via_approximation_new(kerbas1, pmat_sub, shift, pivind1, pivdeg1);

        // recursive call 2, with residual (kerbas * right submatrix of pmat)
        pmat_sub.SetDims(m, n2);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < n2; ++j)
                pmat_sub[i][j] = pmat[i][n1+j];
        multiply(pmat, kerbas1, pmat_sub);
        pmat_sub.kill();

        // recursive call 2
        kernel_basis_zls_via_approximation_new(kerbas, pmat, shift, pivind, pivdeg);

        // multiply bases and combine pivots
        multiply(kerbas, kerbas, kerbas1);
        for (size_t i = 0; i < pivind.size(); ++i)
        {
            pivdeg[i] += pivdeg1[pivind[i]];
            pivind[i] = pivind1[pivind[i]];
        }

        return;
    }

    // Here, we are in the case m>=2 and n <= m/2 (in particular n < m)

    VecLong rdeg; // buffer, used to store shifts/row degrees

    // add a constant to all the shift entries to ensure that the difference
    // diff_shift = min(shift - rdeg(pmat))  is zero
    row_degree(rdeg, pmat);
    long diff_shift = shift[0] - rdeg[0];
    for (long i = 1; i < m; ++i)
    {
        const long diff = shift[i]-rdeg[i];
        if (diff_shift > diff)
            diff_shift = diff;
    }
    std::transform(shift.begin(), shift.end(), shift.begin(),
                   [diff_shift](long x) -> long { return x - diff_shift;});

    // find parameter rho: sum of the n largest entries of (reduced) shift
    rdeg = shift;
    std::sort(rdeg.begin(), rdeg.end());
    const long rho = std::accumulate(rdeg.begin()+m-n, rdeg.end(), 0);

    // order for call to approximation
    // choosing this order (with factor 2) is sufficient to ensure that the
    // approximant basis captures the whole kernel when n <= m/2 and pmat is
    // sufficiently generic
    const long order = 2 * ceil((double)rho / n) + 1;

    // compute approximant basis, along with shift-row degree and pivot degree
    // --> it does not necessarily capture the whole kernel; but does in two
    // notable cases: if n=1, or in the generic situation mentioned above
    Mat<zz_pX> appbas;
    pmbasis_new(appbas, pmat, order, shift, pivdeg); // does not change pmat

    // Identify submatrix of some rows of appbas which are in the kernel
    // Note the criterion: since before the call we have rdeg(pmat) <= shift,
    // the after the call we have rdeg(appbas*pmat) <= shift and therefore rows
    // with shift[i] < order are such that appbas[i] * pmat = 0.
    // Note that this may miss rows in the kernel, if there are some with shift >= order
    // which are in this appbas (this is usually not the case)
    VecLong rdeg1, pivind_app, pivdeg_app;
    VecLong rdeg2, pivind0, pivdeg0;
    for (long i=0; i<m; ++i)
    {
        if (shift[i] < order)
        {
            rdeg1.emplace_back(shift[i]);
            pivdeg_app.emplace_back(pivdeg[i]);
            pivind_app.emplace_back(i);
        }
        else
        {
            rdeg2.emplace_back(shift[i]);
            pivdeg0.emplace_back(pivdeg[i]);
            pivind0.emplace_back(i);
        }
    }

    // number of found kernel rows
    const long m1 = rdeg1.size();
    // number of non-kernel rows in the approximant basis
    const long m2 = pivind0.size();

    // if the sum of pivot degrees in the part of kernel basis computed is
    // equal to the sum of column degrees of pmat, then we know that we have
    // captured the whole kernel
    // --> this is often the case (e.g. in the generic situation mentioned
    // above), and testing this avoids going through a few multiplications and
    // the two recursive calls (which only lead to the conclusion that there is
    // no new kernel row to be found)
    // TODO use row degree? valuation stuff?
    VecLong cdeg; col_degree(cdeg, pmat);
    const bool early_exit =
        (std::accumulate(cdeg.begin(), cdeg.end(),0)
         == std::accumulate(pivdeg_app.begin(), pivdeg_app.end(), 0));
    //std::cout << "\tsum_pivdeg = " << sum_pivdeg << " ; sum cdeg = " << std::accumulate(cdeg.begin(), cdeg.end(),0) << std::endl;

    // if the whole kernel was captured according to the above test, or if
    // there was just one column, or if the kernel is full (matrix was zero),
    // then we have the whole kernel: just copy and return
    if (early_exit || n == 1 || m1 == m)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(appbas[pivind_app[i]]);
        shift = rdeg1;
        pivind = pivind_app;
        pivdeg = pivdeg_app;
        return;
    }

    // retrieve the non-kernel part of the approximant
    Mat<zz_pX> approx(INIT_SIZE, m2, m);
    for (long i = 0; i < m2; ++i)
        approx[i].swap(appbas[pivind0[i]]);

    // compute residual:
    // pmat = trunc(trunc(approx, dA+1)*pmat div x^order, deg(pmat)+deg(approx)-order+1)
    middle_product(pmat,approx,pmat,order,deg(pmat)+deg(approx)-order);

    // pmat will be column-splitted into two submatrices of column dimension ~ n/2
    const long n1 = n/2;
    const long n2 = n-n1;

    // recursive call 1, with left submatrix of the residual pmat
    Mat<zz_pX> pmat_sub(INIT_SIZE, m2, n1);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n1; ++j)
            pmat_sub[i][j] = pmat[i][j];

    Mat<zz_pX> kerbas1;
    VecLong pivind1, pivdeg1;
    kernel_basis_zls_via_approximation_new(kerbas1, pmat_sub, rdeg2, pivind1, pivdeg1);

    // recursive call 2, with right submatrix of the residual pmat
    pmat_sub.SetDims(m2, n2);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n2; ++j)
            pmat_sub[i][j] = pmat[i][n1+j];
    multiply(pmat, kerbas1, pmat_sub);
    pmat_sub.kill();

    // recursive call 2
    VecLong pivind2, pivdeg2;
    kernel_basis_zls_via_approximation_new(kerbas, pmat, rdeg2, pivind2, pivdeg2);

    // if kerbas is empty: (i.e. the approximant basis already captured the
    // whole kernel, although we had not guessed it with early_exit)
    // --> just copy and return
    if (kerbas.NumRows() == 0)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(appbas[pivind_app[i]]);
        shift = rdeg1;
        pivind = pivind_app;
        pivdeg = pivdeg_app;
        return;
    }

    // kerbas is non-empty: we have found new kernel rows
    // I/ complete the computation of these rows via
    // kerbas =  kerbas * kerbas1 * approx
    // II/ update also the pivot index / pivot degree
    // III/ merge this with rows / pivots from the kernel part of appbas

    // I/ complete the computation of new kernel rows
    // We use pmat as a temp, and store the result in approx
    // v1:
    //multiply(pmat, kerbas, kerbas1);
    //multiply(approx, pmat, approx);
    // v2:
    multiply(pmat, kerbas1, approx);
    multiply(approx, kerbas, pmat);
    // TODO when FFT / eval is used, this kind of product may be faster by
    // avoiding interpolation in the middle (?)

    // II/ complete the computation of the pivots of the new kernel rows
    for (size_t i = 0; i < pivind.size(); ++i)
    {
        pivdeg2[i] += (pivdeg0[pivind1[pivind2[i]]] + pivdeg1[pivind2[i]]);
        pivind2[i] = pivind0[pivind1[pivind2[i]]];
    }

    // III/ merge with previously obtained kernel rows, ensuring that
    // the resulting basis is in shifted ordered weak Popov form
    // Note concerning shift: to ensure that shift contains the correct shifted
    // row degree of kerbas, we add the constant diff_shift that had been
    // removed from the input shift
    const long ker_dim = m1 + approx.NumRows();
    kerbas.SetDims(ker_dim, m);
    shift.resize(ker_dim); pivind.resize(ker_dim); pivdeg.resize(ker_dim);
    long i=0, i_appbas=0, i_approx=0;
    while (i_appbas < m1 && i_approx < approx.NumRows())
    {
        if (pivind_app[i_appbas] < pivind2[i_approx])
        {
            // row already captured in preliminary appbas computation
            kerbas[i].swap(appbas[pivind_app[i_appbas]]);
            shift[i] = rdeg1[i_appbas]+diff_shift;
            pivind[i] = pivind_app[i_appbas];
            pivdeg[i] = pivdeg_app[i_appbas];
            ++i_appbas; ++i;
        }
        else
        {
            // row computed via recursive calls and stored in `approx`
            kerbas[i].swap(approx[i_approx]);
            shift[i] = rdeg2[i_approx]+diff_shift;
            pivind[i] = pivind2[i_approx];
            pivdeg[i] = pivdeg2[i_approx];
            ++i_approx; ++i;
        }
    }
    while (i_appbas < m1)
    {
        // row already captured in preliminary appbas computation
        kerbas[i].swap(appbas[pivind_app[i_appbas]]);
        shift[i] = rdeg1[i_appbas]+diff_shift;
        pivind[i] = pivind_app[i_appbas];
        pivdeg[i] = pivdeg_app[i_appbas];
        ++i_appbas; ++i;
    }
    while (i_approx < approx.NumRows())
    {
        // row computed via recursive calls and stored in `approx`
        kerbas[i].swap(approx[i_approx]);
        shift[i] = rdeg2[i_approx]+diff_shift;
        pivind[i] = pivind2[i_approx];
        pivdeg[i] = pivdeg2[i_approx];
        ++i_approx; ++i;
    }

    return;
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
