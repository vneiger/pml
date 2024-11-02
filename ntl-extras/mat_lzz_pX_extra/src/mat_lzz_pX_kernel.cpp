#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_forms.h"
#include "numeric" // for std::iota

#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"
#include <NTL/mat_lzz_p.h>
#include "mat_lzz_pX_kernel.h"

//#define GENERIC_KER_PROFILE

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




bool is_kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift,
                     const PolMatForm & form,
                     const bool randomized
                    )
{
    const long m = pmat.NumRows();

    //std::cout << "~~is_kernel_basis~~ Warning: checking generation not implemented yet" << std::endl;

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


void kernel_basis(
                  Mat<zz_pX> & kerbas,
                  const Mat<zz_pX> & pmat,
                  const VecLong & shift
                 )
{
    VecLong rdeg(shift);
    Mat<zz_pX> copy_mat(pmat);
    kernel_basis_zls_via_approximation(kerbas, copy_mat, rdeg);
}


void kernel_basis_via_approximation(
                                    Mat<zz_pX> & kerbas,
                                    VecLong & pivind,
                                    const Mat<zz_pX> & pmat,
                                    VecLong & shift
                                   )
{
    // parameters
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    const long d = deg(pmat);

    if (d==-1) // matrix is zero
    {
        ident(kerbas, m);
        pivind.resize(m);
        std::iota(pivind.begin(), pivind.end(), 0);
        return;
    }

    // if m==1, the kernel basis is empty (since pmat is nonzero)
    if (m==1)
    {
        kerbas.SetDims(0,1);
        shift.clear();
        pivind.clear();
        return;
    }

    // more generally, if m<=n, check whether coefficient of degree 0 has full
    // rank, if yes then kernel is empty
    if (m <= n)
    {
        Mat<zz_p> coeff0 = coeff(pmat, 0);
        const long r = gauss(coeff0);
        if (r == m)
        {
            kerbas.SetDims(0,m);
            shift.clear();
            pivind.clear();
            return;
        }
    }

    // compute amplitude of the shift
    const long amp = amplitude(shift);

    // compute the order for approximation:
    const long order = (n+1)*d + amp + 1;

    // compute approximant basis and shifted row degree
    Mat<zz_pX> appbas;
    VecLong rdeg(shift);
    pmbasis(appbas, pmat, order, rdeg);

    // find rows which belong to the kernel and record their pivot index
    // (since appbas is in ordered weak Popov form, the pivot index is also
    // the index of the row in appbas)
    pivind.clear();
    for (long i = 0; i < m; ++i)
        if (rdeg[i]-shift[i]+amp < order-d)
            pivind.emplace_back(i);

    // copy the kernel rows and the corresponding shifted row degrees
    const size_t ker_dim = pivind.size();
    shift.resize(ker_dim);
    kerbas.SetDims(ker_dim, m);
    for (size_t i = 0; i < ker_dim; ++i)
    {
        shift[i] = rdeg[pivind[i]];
        kerbas[i].swap(appbas[pivind[i]]);
    }
}

void kernel_basis_via_interpolation(
                                    Mat<zz_pX> & kerbas,
                                    VecLong & pivind,
                                    const Mat<zz_pX> & pmat,
                                    VecLong & shift
                                   )
{
    // parameters
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    const long d = deg(pmat);

    if (d==-1) // matrix is zero
    {
        ident(kerbas, m);
        pivind.resize(m);
        std::iota(pivind.begin(), pivind.end(), 0);
        return;
    }

    // if m==1, the kernel basis is empty (since pmat is nonzero)
    if (m==1)
    {
        kerbas.SetDims(0,1);
        shift.clear();
        pivind.clear();
        return;
    }

    // more generally, if m<=n, check whether coefficient of degree 0 has full
    // rank, if yes then kernel is empty
    if (m <= n)
    {
        Mat<zz_p> coeff0 = coeff(pmat, 0);
        const long r = gauss(coeff0);
        if (r == m)
        {
            kerbas.SetDims(0,m);
            shift.clear();
            pivind.clear();
            return;
        }
    }

    // compute amplitude of the shift
    const long amp = amplitude(shift);

    // compute the order for approximation:
    const long order = (n+1)*d + amp + 1;

    // define geometric points r^2, r^4, r^6...
    zz_p r;
    if (2*order >= zz_p::modulus()) // small field --> approximation
    {
        kernel_basis_via_approximation(kerbas, pivind, pmat, shift);
        return;
    }

    element_of_order(r, 2*order);
    if (IsZero(r))
    {
        std::cout << "Could not find element of sufficiently large order; field is probably too small for kernel via interpolation with such degrees --> calling the approximation version" << std::endl;
        kernel_basis_via_approximation(kerbas, pivind, pmat, shift);
        return;
    }

    // compute interpolant basis and shifted row degree
    VecLong rdeg(shift);
    Mat<zz_pX> intbas;
    Vec<zz_p> pts;
    pmbasis_geometric(intbas,pmat,r,order,rdeg,pts);

    // find rows which belong to the kernel and record their pivot index
    // (since appbas is in ordered weak Popov form, the pivot index is also
    // the index of the row in appbas)
    pivind.clear();
    for (long i = 0; i < m; ++i)
        if (rdeg[i]-shift[i]+amp < order-d)
            pivind.emplace_back(i);

    // copy the kernel rows and the corresponding shifted row degrees
    const size_t ker_dim = pivind.size();
    shift.resize(ker_dim);
    kerbas.SetDims(ker_dim, m);
    for (size_t i = 0; i < ker_dim; ++i)
    {
        shift[i] = rdeg[pivind[i]];
        kerbas[i].swap(intbas[pivind[i]]);
    }
}


void kernel_basis_zls_via_approximation(
                                        Mat<zz_pX> & kerbas,
                                        Mat<zz_pX> & pmat,
                                        VecLong & shift
                                       )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // if m==1, just check whether pmat is zero
    if (m==1)
    {
        if (not IsZero(pmat))
        {
            kerbas.SetDims(0,1);
            shift.clear();
            return;
        }
        // pmat is the zero 1 x n matrix
        // --> kerbas is [1], shift is unchanged
        ident(kerbas, 1);
        return;
    }

    // if m<=n, check whether coefficient of degree 0 has full rank, if yes
    // then kernel is empty
    if (m <= n)
    {
        Mat<zz_p> coeff0 = coeff(pmat, 0);
        const long r = gauss(coeff0);
        if (r == m)
        {
            kerbas.SetDims(0,m);
            shift.clear();
            return;
        }
    }

    // compute row degree of pmat (will be used afterwards for making the
    // shift such that shift >= rdeg)
    const VecLong rdeg = row_degree(pmat);

    // if matrix is zero, kernel is identity (shift unchanged)
    if (*max_element(rdeg.begin(), rdeg.end()) == -1)
    {
        ident(kerbas, m);
        return;
    }
    // TODO handle constant case with NTL's kernel? (involves permutations to
    // respect the shift... not so easy... is the current algo much slower than
    // NTL's kernel in that case (it does rely on it via pmbasis)?

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
        kernel_basis_zls_via_approximation(kerbas1, pmat_sub, shift);

        // if kerbas1 is empty, just return empty kernel
        if (kerbas1.NumRows() == 0)
        {
            kerbas.SetDims(0, m);
            shift.clear();
            return;
        }

        // recursive call 2, with residual (kerbas * right submatrix of pmat)
        pmat_sub.SetDims(m, n2);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < n2; ++j)
                pmat_sub[i][j] = pmat[i][n1+j];
        multiply(pmat, kerbas1, pmat_sub);
        pmat_sub.kill();

        // recursive call 2
        kernel_basis_zls_via_approximation(kerbas, pmat, shift);

        // multiply bases
        multiply(kerbas, kerbas, kerbas1);
        return;
    }

    // Here, we are in the case m>=2 and n <= m/2 (in particular n < m)

    // add a constant to all the shift entries to ensure that the difference
    // diff_shift = min(shift - rdeg(pmat))  is zero
    // Recall that currently rdeg is rdeg(pmat)
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
    VecLong rdeg1(shift); // rdeg1 will be re-used later for another purpose, hence the name
    std::sort(rdeg1.begin(), rdeg1.end());
    const long rho = std::accumulate(rdeg1.begin()+m-n, rdeg1.end(), 0);

    // order for call to approximation
    // choosing this order (with factor 2) is sufficient to ensure that the
    // approximant basis captures the whole kernel when n <= m/2 and pmat is
    // sufficiently generic
    const long order = 2 * ceil((double)rho / n) + 1;

    // compute approximant basis, along with shift-row degree and pivot degree
    // --> it does not necessarily capture the whole kernel; but does in two
    // notable cases: if n=1, or in the generic situation mentioned above

    // copy of input shift, will be needed at the end to compute shifted pivot index
    VecLong copy_shift(shift);
    Mat<zz_pX> appbas;
    pmbasis(appbas, pmat, order, shift); // does not modify pmat

    // Identify submatrix of some rows of appbas which are in the kernel
    // Note the criterion: since before the call we have rdeg(pmat) <= shift,
    // the after the call we have rdeg(appbas*pmat) <= shift and therefore rows
    // with shift[i] < order are such that appbas[i] * pmat = 0.
    // Note that this may miss rows in the kernel, if there are some with shift >= order
    // which are in this appbas (this is usually not the case)
    rdeg1.clear();
    VecLong rows_app_ker;
    VecLong rdeg2, rows_app_notker;
    for (long i=0; i<m; ++i)
    {
        if (shift[i] < order)
        {
            rdeg1.emplace_back(shift[i]);
            rows_app_ker.emplace_back(i);
        }
        else
        {
            rdeg2.emplace_back(shift[i]);
            rows_app_notker.emplace_back(i);
        }
    }

    // number of found kernel rows
    const long m1 = rdeg1.size();
    // number of non-kernel rows in the approximant basis
    const long m2 = rdeg2.size();

    // if the sum of pivot degrees in the part of kernel basis computed is
    // equal to one of the upper bounds (here, sum of row degrees of pmat for
    // non-pivot rows), then we know that we have captured the whole kernel
    // --> this is often the case (e.g. in the generic situation mentioned
    // above), and testing this avoids going through a few multiplications and
    // the two recursive calls (which only lead to the conclusion that there is
    // no new kernel row to be found)
    long sum_matdeg = 0;
    for (long i = 0; i < m2; ++i)
        sum_matdeg += rdeg[rows_app_notker[i]];
    long sum_kerdeg = 0;
    for (long i = 0; i < m1; ++i)
        sum_kerdeg += deg(appbas[rows_app_ker[i]][rows_app_ker[i]]);

    const bool early_exit = (sum_kerdeg == sum_matdeg);

    // if the whole kernel was captured according to the above test, or if
    // there was just one column, or if the kernel is full (matrix was zero),
    // then we have the whole kernel: just copy and return
    if (early_exit || n == 1 || m1 == m)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(appbas[rows_app_ker[i]]);
        shift.swap(rdeg1);
        return;
    }

    // Note: for most input matrices/shifts, the following code will not be
    // executed because we will have taken early exit

    // retrieve the non-kernel part of the approximant
    Mat<zz_pX> approx(INIT_SIZE, m2, m);
    for (long i = 0; i < m2; ++i)
        approx[i].swap(appbas[rows_app_notker[i]]);

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
    kernel_basis_zls_via_approximation(kerbas1, pmat_sub, rdeg2);

    // recursive call 2, with right submatrix of the residual pmat
    pmat_sub.SetDims(m2, n2);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n2; ++j)
            pmat_sub[i][j] = pmat[i][n1+j];
    multiply(pmat, kerbas1, pmat_sub);
    pmat_sub.kill();

    kernel_basis_zls_via_approximation(kerbas, pmat, rdeg2);

    // if kerbas is empty: (i.e. the approximant basis already captured the
    // whole kernel, although we had not guessed it with early_exit)
    // --> just copy and return
    if (kerbas.NumRows() == 0)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(appbas[rows_app_ker[i]]);
        shift.swap(rdeg1);
        return;
    }

    // kerbas is non-empty: we have found new kernel rows
    // I/ complete the computation of these rows via
    //         kerbas =  kerbas * kerbas1 * approx
    // II/ merge this with rows from the kernel part of appbas

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

    // II/ merge with previously obtained kernel rows, ensuring that
    // the resulting basis is in shifted ordered weak Popov form

    // Note concerning shift: to ensure that shift contains the correct shifted
    // row degree of kerbas, we add the constant diff_shift that had been
    // removed from the input shift

    // Note that `rows_app_ker` is the shifted pivot index of the kernel part
    // of `appbas`
    // --> we need to compute the shifted pivot index of the newly found rows
    // of the kernel, as follows (note that copy_shift is a copy of the input
    // shift, while here `shift` is used as a temporary variable which will
    // store pivot degrees, which we discard):
    VecLong pivots_approx;
    row_pivots(pivots_approx, shift, approx, copy_shift);

    const long ker_dim = m1 + approx.NumRows();
    kerbas.SetDims(ker_dim, m);
    shift.resize(ker_dim);
    long i=0, i_appbas=0, i_approx=0;
    while (i_appbas < m1 && i_approx < approx.NumRows())
    {
        if (rows_app_ker[i_appbas] < pivots_approx[i_approx])
        {
            // row already captured in preliminary appbas computation
            kerbas[i].swap(appbas[rows_app_ker[i_appbas]]);
            shift[i] = rdeg1[i_appbas]+diff_shift;
            ++i_appbas; ++i;
        }
        else
        {
            // row computed via recursive calls and stored in `approx`
            kerbas[i].swap(approx[i_approx]);
            shift[i] = rdeg2[i_approx]+diff_shift;
            ++i_approx; ++i;
        }
    }
    while (i_appbas < m1)
    {
        // row already captured in preliminary appbas computation
        kerbas[i].swap(appbas[rows_app_ker[i_appbas]]);
        shift[i] = rdeg1[i_appbas]+diff_shift;
        ++i_appbas; ++i;
    }
    while (i_approx < approx.NumRows())
    {
        // row computed via recursive calls and stored in `approx`
        kerbas[i].swap(approx[i_approx]);
        shift[i] = rdeg2[i_approx]+diff_shift;
        ++i_approx; ++i;
    }
}

void kernel_basis_zls_via_interpolation(
                                        Mat<zz_pX> & kerbas,
                                        Mat<zz_pX> & pmat,
                                        VecLong & shift
                                       )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // if m==1, just check whether pmat is zero
    if (m==1)
    {
        if (not IsZero(pmat))
        {
            kerbas.SetDims(0,1);
            shift.clear();
            return;
        }
        // pmat is the zero 1 x n matrix
        // --> kerbas is [1], shift is unchanged
        ident(kerbas, 1);
        return;
    }

    // if m<=n, check whether coefficient of degree 0 has full rank, if yes
    // then kernel is empty
    if (m <= n)
    {
        Mat<zz_p> coeff0 = coeff(pmat, 0);
        const long r = gauss(coeff0);
        if (r == m)
        {
            kerbas.SetDims(0,m);
            shift.clear();
            return;
        }
    }

    // compute row degree of pmat (will be used afterwards for making the
    // shift such that shift >= rdeg)
    const VecLong rdeg = row_degree(pmat);

    // if matrix is zero, kernel is identity (shift unchanged)
    if (*max_element(rdeg.begin(), rdeg.end()) == -1)
    {
        ident(kerbas, m);
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
        kernel_basis_zls_via_interpolation(kerbas1, pmat_sub, shift);

        // if kerbas1 is empty, just return empty kernel
        if (kerbas1.NumRows() == 0)
        {
            kerbas.SetDims(0, m);
            shift.clear();
            return;
        }

        // recursive call 2, with residual (kerbas * right submatrix of pmat)
        pmat_sub.SetDims(m, n2);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < n2; ++j)
                pmat_sub[i][j] = pmat[i][n1+j];
        multiply(pmat, kerbas1, pmat_sub);
        pmat_sub.kill();

        // recursive call 2
        kernel_basis_zls_via_interpolation(kerbas, pmat, shift);

        // multiply bases
        multiply(kerbas, kerbas, kerbas1);
        return;
    }

    // Here, we are in the case m>=2 and n <= m/2 (in particular n < m)

    // add a constant to all the shift entries to ensure that the difference
    // diff_shift = min(shift - rdeg(pmat))  is zero
    // Recall that currently rdeg is rdeg(pmat)
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
    VecLong rdeg1(shift); // rdeg1 will be re-used later for another purpose, hence the name
    std::sort(rdeg1.begin(), rdeg1.end());
    const long rho = std::accumulate(rdeg1.begin()+m-n, rdeg1.end(), 0);

    // order for call to interpolation
    // choosing this order (with factor 2) is sufficient to ensure that the
    // interpolant basis captures the whole kernel when n <= m/2 and pmat is
    // sufficiently generic
    const long order = 2 * ceil((double)rho / n) + 1;

    // compute interpolant basis, along with shift-row degree and pivot degree
    // --> it does not necessarily capture the whole kernel; but does in two
    // notable cases: if n=1, or in the generic situation mentioned above
    zz_p r; // defining geometric points r^2, r^4, r^6...
    if (2*order >= zz_p::modulus()) // small field --> approximation
    {
        kernel_basis_zls_via_approximation(kerbas, pmat, shift);
        return;
    }

    element_of_order(r, 2*order);
    if (IsZero(r))
    {
        std::cout << "Could not find element of sufficiently large order; field is probably too small for kernel via interpolation with such degrees --> calling the approximation version" << std::endl;
        kernel_basis_zls_via_approximation(kerbas, pmat, shift);
        return;
    }

    // copy of input shift, will be needed at the end to compute shifted pivot index
    VecLong copy_shift(shift);
    Mat<zz_pX> intbas;
    Vec<zz_p> pts;
    pmbasis_geometric(intbas,pmat,r,order,shift,pts); // does not modify pmat

    // Identify submatrix of some rows of intbas which are in the kernel
    // Note the criterion: since before the call we have rdeg(pmat) <= shift,
    // the after the call we have rdeg(intbas*pmat) <= shift and therefore rows
    // with shift[i] < order are such that intbas[i] * pmat = 0.
    // Note that this may miss rows in the kernel, if there are some with shift >= order
    // which are in this intbas (this is usually not the case)
    rdeg1.clear();
    VecLong rows_int_ker;
    VecLong rdeg2, rows_int_notker;
    for (long i=0; i<m; ++i)
    {
        if (shift[i] < order)
        {
            rdeg1.emplace_back(shift[i]);
            rows_int_ker.emplace_back(i);
        }
        else
        {
            rdeg2.emplace_back(shift[i]);
            rows_int_notker.emplace_back(i);
        }
    }

    // number of found kernel rows
    const long m1 = rdeg1.size();
    // number of non-kernel rows in the interpolant basis
    const long m2 = rows_int_notker.size();

    // if the sum of pivot degrees in the part of kernel basis computed is
    // equal to one of the upper bounds (here, sum of row degrees of pmat for
    // non-pivot rows), then we know that we have captured the whole kernel
    // --> this is often the case (e.g. in the generic situation mentioned
    // above), and testing this avoids going through a few multiplications and
    // the two recursive calls (which only lead to the conclusion that there is
    // no new kernel row to be found)
    long sum_matdeg = 0;
    for (long i = 0; i < m2; ++i)
        sum_matdeg += rdeg[rows_int_notker[i]];
    long sum_kerdeg = 0;
    for (long i = 0; i < m1; ++i)
        sum_kerdeg += deg(intbas[rows_int_ker[i]][rows_int_ker[i]]);

    const bool early_exit = (sum_kerdeg == sum_matdeg);

    // if the whole kernel was captured according to the above test, or if
    // there was just one column, or if the kernel is full (matrix was zero),
    // then we have the whole kernel: just copy and return
    if (early_exit || n == 1 || m1 == m)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(intbas[rows_int_ker[i]]);
        shift.swap(rdeg1);
        return;
    }

    // Note: for most input matrices/shifts, the following code will not be
    // executed because we will have taken early exit

    // retrieve the non-kernel part of the interpolant
    Mat<zz_pX> interp(INIT_SIZE, m2, m);
    for (long i = 0; i < m2; ++i)
        interp[i].swap(intbas[rows_int_notker[i]]);

    // compute residual:
    // pmat = (pmat*interp) div (x-pts[0])...(x-pts[order-1])

    // construct the modulus polynomial (x-pts[0]) ... (x-pts[order-1])
    const zz_pX mod = BuildFromRoots(pts);
    zz_pXModulus Mod(mod);
    multiply(pmat,interp,pmat);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n; ++j)
            div(pmat[i][j], pmat[i][j], Mod);

    // pmat will be column-splitted into two submatrices of column dimension ~ n/2
    const long n1 = n/2;
    const long n2 = n-n1;

    // recursive call 1, with left submatrix of the residual pmat
    Mat<zz_pX> pmat_sub(INIT_SIZE, m2, n1);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n1; ++j)
            pmat_sub[i][j] = pmat[i][j];

    Mat<zz_pX> kerbas1;
    kernel_basis_zls_via_interpolation(kerbas1, pmat_sub, rdeg2);

    // recursive call 2, with right submatrix of the residual pmat
    pmat_sub.SetDims(m2, n2);
    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n2; ++j)
            pmat_sub[i][j] = pmat[i][n1+j];
    multiply(pmat, kerbas1, pmat_sub);
    pmat_sub.kill();

    kernel_basis_zls_via_interpolation(kerbas, pmat, rdeg2);

    // if kerbas is empty: (i.e. the interpolant basis already captured the
    // whole kernel, although we had not guessed it with early_exit)
    // --> just copy and return
    if (kerbas.NumRows() == 0)
    {
        kerbas.SetDims(m1, m);
        for (long i = 0; i < m1; ++i)
            kerbas[i].swap(intbas[rows_int_ker[i]]);
        shift.swap(rdeg1);
        return;
    }

    // kerbas is non-empty: we have found new kernel rows
    // I/ complete the computation of these rows via
    // kerbas =  kerbas * kerbas1 * interp
    // II/ merge this with rows from the kernel part of intbas

    // I/ complete the computation of new kernel rows
    // We use pmat as a temp, and store the result in interp
    // v1:
    //multiply(pmat, kerbas, kerbas1);
    //multiply(interp, pmat, interp);
    // v2:
    multiply(pmat, kerbas1, interp);
    multiply(interp, kerbas, pmat);
    // TODO when FFT / eval is used, this kind of product may be faster by
    // avoiding interpolation in the middle (?)

    // II/ merge with previously obtained kernel rows, ensuring that
    // the resulting basis is in shifted ordered weak Popov form

    // Note concerning shift: to ensure that shift contains the correct shifted
    // row degree of kerbas, we add the constant diff_shift that had been
    // removed from the input shift

    // Note that `rows_int_ker` is the shifted pivot index of the kernel part
    // of `intbas`
    // --> we need to compute the shifted pivot index of the newly found rows
    // of the kernel, as follows (note that copy_shift is a copy of the input
    // shift, while here `shift` is used as a temporary variable which will
    // store pivot degrees, which we discard):
    VecLong pivots_interp;
    row_pivots(pivots_interp, shift, interp, copy_shift);

    const long ker_dim = m1 + interp.NumRows();
    kerbas.SetDims(ker_dim, m);
    shift.resize(ker_dim);
    long i=0, i_intbas=0, i_interp=0;
    while (i_intbas < m1 && i_interp < interp.NumRows())
    {
        if (rows_int_ker[i_intbas] < pivots_interp[i_interp])
        {
            // row already captured in preliminary intbas computation
            kerbas[i].swap(intbas[rows_int_ker[i_intbas]]);
            shift[i] = rdeg1[i_intbas]+diff_shift;
            ++i_intbas; ++i;
        }
        else
        {
            // row computed via recursive calls and stored in `interp`
            kerbas[i].swap(interp[i_interp]);
            shift[i] = rdeg2[i_interp]+diff_shift;
            ++i_interp; ++i;
        }
    }
    while (i_intbas < m1)
    {
        // row already captured in preliminary intbas computation
        kerbas[i].swap(intbas[rows_int_ker[i_intbas]]);
        shift[i] = rdeg1[i_intbas]+diff_shift;
        ++i_intbas; ++i;
    }
    while (i_interp < interp.NumRows())
    {
        // row computed via recursive calls and stored in `interp`
        kerbas[i].swap(interp[i_interp]);
        shift[i] = rdeg2[i_interp]+diff_shift;
        ++i_interp; ++i;
    }
}


// input: (n+m) x n polynomial matrix of the form
//     [ X**d I + A(X) ]
//     [      V(X)     ]
// where A is n x n and V is m x n, both of degree < d
// assumption: generic left kernel basis degrees
// output: 0-weak Popov kernel basis of pmat
// If d <= 1, kerbas is replaced with the kernel basis
// If d > 1, kerbas is multiplied by the kernel basis (so kerbas should
// already be allocated with appropriate dimensions)
void kernel_basis_generic(Mat<zz_pX> & kerbas, Mat<zz_pX> & pmat)
{
    const long n = pmat.NumCols();
    const long m = pmat.NumRows() - n;
    const long d = deg(pmat[0][0]);
    //std::cout << degree_matrix(pmat) << std::endl;

#ifdef GENERIC_KER_PROFILE
    double t;
    std::cout << "enter with dims = " << m+n << " x " << n << ", deg = " << d << std::endl;
#endif // GENERIC_KER_PROFILE

    if (d == 0)
    {
        Mat<zz_p> K;
        kernel(K, coeff(pmat, 0));
        conv(kerbas, K);
    }
    else if (m >= n)
    {
        // base case: one call to pmbasis
        VecLong shift = row_degree(pmat);

        // compute the kernel via approximant basis at high order
        Mat<zz_pX> appbas;
        // degree of kernel basis will be d (generically)
        // --> compute approximants at order 2d is enough
        // (one might need 2d+1, but thanks to the weak Popov form, 2d suffices)
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime();
#endif // GENERIC_KER_PROFILE
        popov_pmbasis_generic(appbas, pmat, 2*d);
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime()-t;
        std::cout << "\tpmbasis order = " << 2*d << " || time " << t << std::endl;
#endif // GENERIC_KER_PROFILE

        // minimal left kernel basis of pmat: last m rows of app
        Mat<zz_pX> kerbas2;
        kerbas2.SetDims(m, m+n);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < m+n; ++j)
                kerbas2[i][j] = appbas[n+i][j];

        // compute product of kernel bases
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime();
#endif // GENERIC_KER_PROFILE
        multiply(kerbas, kerbas2, kerbas);
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime()-t;
        std::cout << "\tmultiply ker, degrees " << deg(kerbas2) << "," << deg(kerbas) << " || time " << t << std::endl;
#endif // GENERIC_KER_PROFILE
    }
    else if (d == 1)
    {
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime();
#endif // GENERIC_KER_PROFILE
        // case of characteristic matrix:
        //     [ X I + A ]
        //     [    V    ]
        // where A and V are constant, and m < n
        const long nn1 = (m+n)/2;
        const long nn2 = n - nn1;
        const long eps = m + nn2 - nn1;

        // Let F be the submatrix formed by the left nn1 columns of pmat and
        // compute its left kernel; these columns are
        //             nn1
        //         [ X I + A1 ]  nn1
        //     F = [    A2    ]  nn2
        //         [     V    ]   m 
        //
        // The choice nn1 == floor((m+n)/2) <= (m+n)/2 == NumRows/2
        // guarantees a degree one kernel basis (under the genericity assumption)
        // of the form [ K1   K2 + X*K3 ], dimensions (m+nn2) x (m+n),
        // but we can be more precise.
        //
        // Also, note m + nn2 == m + n - nn1 == ceil((m+n)/2)
        // is either nn1 or nn1+1. Let eps = m + nn2 - nn1.
        // Note also nn1+nn1+eps = nn1+m+nn2 = m+n.
        //
        // More precisely, let
        //           nn1  nn1  eps
        //         [ K1    I   0 ]  nn1
        //         [ K2    0   1 ]  eps
        // be the generic kernel of the constant part [[A1], [A2], [V]].
        // Then the matrix
        //      K =   [ -B   -B*K1i + X I  0]
        //            [ K2         0       1]
        // where B = F[nn1:2*nn1, :] == A[nn1:2*nn1, :nn1]
        // and K1i == inverse(K1)
        // is the Popov left kernel basis for F
        // -> it has the right rank
        // -> it is in Popov form
        // -> it has the right sum(rdeg) == nn1
        // -> K*F == 0:
        //    --> clear for the row [K2  0  1]
        //    --> the other one yields
        //           -B * (X I + A1) + (-B*K1i + X I) * B
        //           == X * (-B + B) + (-B*A1 - B*K1i*B)
        //           ==                -B*K1i * (K1 * A1 + B)
        //           == 0.

        Mat<zz_p> F0;
        F0.SetDims(m+n, nn1);
        for (long i = 0; i < m+n; i++)
            for (long j = 0; j < nn1; j++)
                F0[i][j].LoopHole() = ConstTerm(pmat[i][j])._zz_p__rep;

        Mat<zz_p> K;
        kernel(K, F0);

        Mat<zz_p> mB;  // this is -B
        mB.SetDims(nn1, nn1);
        for (long i = 0; i < nn1; i++)
            NTL::negate(mB[i], F0[nn1+i]);
        F0.kill();

        Mat<zz_p> K1, K2;
        K1.SetDims(nn1, nn1);
        K2.SetDims(eps, nn1);
        for (long i = 0; i < nn1; i++)
            VectorCopy(K1[i], K[i], nn1);
        if (eps > 0)
            VectorCopy(K2[0], K[nn1], nn1);

        inv(K1, K1);
        mul(K1, mB, K1);

        // FIXME the first multiplication (kerbas2 in rec call by the below
        // kerbas) could be accelerated
        // currently we just build kerbas explicitly now, and multiply later in recursive call
        kerbas.SetDims(nn1+eps, m+n);
        for (long i = 0; i < nn1; i++)
        {
            for (long j = 0; j < nn1; j++)
                SetCoeff(kerbas[i][j], 0, mB[i][j]);
            for (long j = nn1; j < 2*nn1; j++)
            {
                if (i+nn1 == j)
                    SetX(kerbas[i][j]);
                SetCoeff(kerbas[i][j], 0, K1[i][j-nn1]);
            }
        }
        if (eps > 0)
        {
            for (long j = 0; j < nn1; j++)
                SetCoeff(kerbas[nn1][j], 0, K2[0][j]);
            set(kerbas[nn1][m+n-1]); 
        }

        // multiply kernel basis with right-hand columns
        //  
        //      K =   [ -B   -B*K1i + X I  0]
        //            [ K2         0       1]
        //      G =   [  A3  ]     [ 0 ]
        //            [  A4  ] + X [ L ]
        //            [  v2  ]     [ 0 ]
        //     with L = [[Id nn2 x nn2], [Zero (nn1-nn2) x nn2]], nn1-nn2 == m-eps
        // - the possible last row is K2 * A3 + v2
        // - the first nn1 rows are mB*A3 + mB*K1i*A4 + X*(A4 + mB*K1i*L) + X**2 * L
        // recall mB*K1i is currently stored in K1
        Mat<zz_pX> F;
        F.SetDims(nn1+eps, nn2);

        // retrieve A4 and do degree 1 and degree 2 terms
        Mat<zz_p> A4;
        A4.SetDims(nn1, nn2);
        for (long i = 0; i < nn1; i++)
        {
            if (i < nn2)
                SetCoeff(F[i][i], 2);
            for (long j = 0; j < nn2; j++)
            {
                zz_p c = ConstTerm(pmat[nn1+i][nn1+j]);
                A4[i][j] = c;
                SetCoeff(F[i][j], 1, c + K1[i][j]);
            }
        }

        // retrieve A3 and do degree 0 term
        Mat<zz_p> A3;
        A3.SetDims(nn1, nn2);
        for (long i = 0; i < nn1; i++)
            for (long j = 0; j < nn2; j++)
                A3[i][j].LoopHole() = ConstTerm(pmat[i][nn1+j])._zz_p__rep;

        // handle possible last row
        if (eps > 0)
        {
            mul(K2, K2, A3);
            for (long j = 0; j < nn2; j++)
                SetCoeff(F[nn1][j], 0, K2[0][j] + ConstTerm(pmat[m+n-1][j]));
        }

        mul(A3, mB, A3);

        mul(A4, K1, A4);
        add(A3, A4, A3);
        for (long i = 0; i < nn1; i++)
            for (long j = 0; j < nn2; j++)
                SetCoeff(F[i][j], 0, A3[i][j]);
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime()-t;
        std::cout << "\tdegree 1, main time " << t << std::endl;
#endif // GENERIC_KER_PROFILE

        F0.kill();
        A3.kill();
        A4.kill();
        K.kill();
        K1.kill();
        K2.kill();
        mB.kill();

        //std::cout << kerbas*pmat << std::endl;
        //std::cout << F << std::endl;
        kernel_basis_generic(kerbas, F);
    }
    else
    {
        // assuming we have gone through degree 1 case,
        // which brought us to almost square case;
        // otherwise for rectangular input this might be accelerated
        const long nn1 = (m+n)/2;
        const long nn2 = n-nn1;

        Mat<zz_pX> pmat_l;
        Mat<zz_pX> pmat_r;
        pmat_l.SetDims(m+n,nn1);
        pmat_r.SetDims(m+n,nn2);

        for (long i = 0; i < m+n; ++i)
        {
            for (long j = 0; j < nn1; ++j)
                pmat_l[i][j] = pmat[i][j];
            for (long j = 0; j < nn2; ++j)
                pmat_r[i][j] = pmat[i][j+nn1];
        }

        // compute the kernel of pmat_l via approximant basis at high order
        Mat<zz_pX> appbas;
        // by choice of nn1, degree of kernel basis will be d (generically)
        // --> compute approximants at order 2d is enough
        // (one might need 2d+1, but thanks to the weak Popov form, 2d suffices)
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime();
#endif // GENERIC_KER_PROFILE
        popov_pmbasis_generic(appbas, pmat_l, 2*d);
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime()-t;
        std::cout << "\tpmbasis order = " << 2*d << " || time " << t << std::endl;
#endif // GENERIC_KER_PROFILE

        // minimal left kernel basis of pmat_l : last rows of app
        Mat<zz_pX> kerbas2;
        kerbas2.SetDims(m+nn2, m+n);
        for (long i = 0; i < m+nn2; ++i)
            for (long j = 0; j < m+n; ++j)
                kerbas2[i][j] = appbas[nn1+i][j];

        // then compute product with right-hand columns
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime();
#endif // GENERIC_KER_PROFILE
        multiply(pmat, kerbas2, pmat_r);
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime()-t;
        std::cout << "\tmultiply mat, degrees " << deg(kerbas2) << "," << deg(pmat_r) << " || time " << t << std::endl;
#endif // GENERIC_KER_PROFILE

        // and product of kernel bases
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime();
#endif // GENERIC_KER_PROFILE
        multiply(kerbas, kerbas2, kerbas);
#ifdef GENERIC_KER_PROFILE
        t = GetWallTime()-t;
        std::cout << "\tmultiply ker, degrees " << deg(kerbas2) << "," << deg(kerbas) << " || time " << t << std::endl;
#endif // GENERIC_KER_PROFILE

        pmat_l.kill();
        pmat_r.kill();
        appbas.kill();
        kerbas2.kill();

        //cout << kerbas2 * pmat_l << endl;
        //cout << degree_matrix(kerbas2) << endl;
        //cout << degree_matrix(appbas) << endl;

        // finally, the recursive call
        kernel_basis_generic(kerbas, pmat);
    }

    return;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
