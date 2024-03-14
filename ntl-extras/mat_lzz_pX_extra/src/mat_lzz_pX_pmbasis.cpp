#include <numeric> // for std::iota

#include "lzz_pX_middle_product.h"
#include "mat_lzz_pX_approximant.h"

//#define MBASIS1_PROFILE
//#define MBASIS_PROFILE
//#define PMBASIS_PROFILE

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MBASIS: ITERATIVE ALGORITHMS FOR UNIFORM ORDER             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// These algorithms are defined for arbitrary shifts, but work best for shifts
// close to uniform
// (except popov_mbasis1/mbasis1, where the shift has roughly no influence)

/*------------------------------------------------------------*/
/* base case: order 1, i.e. working modulo X                  */
/* returning Popov form                                       */
/*------------------------------------------------------------*/

VecLong popov_mbasis1(
                      Mat<zz_p> & kerbas,
                      const Mat<zz_p> & pmat,
                      const VecLong & shift
                     )
{
    long m = pmat.NumRows();
    long n = pmat.NumCols();
    // compute permutation which realizes stable sort of the shift
    // (i.e. sorts (shift[0],0)....(shift[len],len) lexicographically increasingly)
#ifdef MBASIS1_PROFILE
    double t_perm1,t_perm2,t_pivind,t_ker,t_now;
    t_now = GetWallTime();
#endif
    VecLong perm_shift(m);
    std::iota(perm_shift.begin(), perm_shift.end(), 0);
    stable_sort(perm_shift.begin(), perm_shift.end(),
                [&](const long& a, const long& b)->bool
                {
                return (shift[a] < shift[b]);
                } );

    // permute rows of pmat accordingly
    Mat<zz_p> mat;
    mat.SetDims(m,n);
    for (long i = 0; i < m; ++i)
        mat[i] = pmat[perm_shift[i]];
#ifdef MBASIS1_PROFILE
    t_perm1 = GetWallTime() - t_now;
#endif

    // find the permuted kernel basis in row echelon form
#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    Mat<zz_p> p_kerbas;
    kernel(p_kerbas,mat);
#ifdef MBASIS1_PROFILE
    t_ker = GetWallTime() - t_now;
#endif
    long k = p_kerbas.NumRows();
    if (k==0)
    {
        // kerbas is empty, same as p_kerbas
        kerbas.move(p_kerbas);
        return VecLong(m,1);
    }
    if (k==m)
    {
        // kerbas should be the identity, same as p_kerbas
        kerbas.move(p_kerbas);
        return VecLong(m,0);
    }

#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    // Check whether we have the expected pivot indices in permuted kernel
    // (pivot = rightmost nonzero entry)
    // NTL ensures that (according to experiments):
    //   * kernel is expected to be of the form [ K | Id ]
    //   * in general it is a column-permutation of such a matrix
    // However note that a column-permutation is not sufficient for our needs.
    // Another property: if pivots are in the expected location (diagonal of
    // rightmost square submatrix), then the corresponding column is the identity column.
    bool expected_pivots = true;
    VecLong p_pivind(k,m-1);
    for (long i = 0; i<k; ++i)
    {
        while (IsZero(p_kerbas[i][p_pivind[i]]))
            --p_pivind[i];
        if (p_pivind[i] != m-k+i)
        { expected_pivots=false; break; }
    }

    if (not expected_pivots)
    {
        // the kernel is not in the wanted shape
        // --> let's compute its lower triangular reduced row echelon form
        Mat<zz_p> column_permuted_ker;
        column_permuted_ker.SetDims(k,m);
        for (long i = 0; i < k; ++i)
            for (long j = 0; j < m; ++j)
                column_permuted_ker[i][j] = p_kerbas[i][m-1-j];
        image(column_permuted_ker, column_permuted_ker);
        // now column_permuted_ker is in upper triangular row echelon form
        for (long i = 0; i < k; ++i)
            for (long j = 0; j < m; ++j)
                p_kerbas[i][j] = column_permuted_ker[k-i-1][m-1-j];
        // and now p_kerbas is in lower triangular row echelon form
        // it remains to make it reduced
        // --> we use a rather naive way, which should not impact timings
        // too much for most use cases, since this case (non-generic kernel) is
        // rare over large fields with not-too-special input:

        // compute the actual pivot indices
        for (long i = 0; i<k; ++i)
        {
            p_pivind[i] = m-1;
            while (IsZero(p_kerbas[i][p_pivind[i]]))
                --p_pivind[i];
        }

        Mat<zz_p> tmat;
        tmat.SetDims(k,k);
        for (long i = 0; i < k; ++i)
            for (long j = 0; j < k; ++j)
                tmat[i][j] = p_kerbas[i][p_pivind[j]];
        inv(tmat, tmat);
        mul(p_kerbas, tmat, p_kerbas);
    }

    // now p_kerbas is in lower triangular reduced row echelon form
    // permute everything back to original order:
    // prepare kernel permutation by permuting kernel pivot indices
    // pivot degrees corresponding to kernel pivot indices are 0, others are 1
    VecLong pivind(k);
    VecLong pivdeg(m,1);
    for (long i = 0; i < k; ++i)
    {
        pivind[i] = perm_shift[p_pivind[i]];
        pivdeg[pivind[i]] = 0;
    }

    VecLong perm_rows_ker(k);
    std::iota(perm_rows_ker.begin(), perm_rows_ker.end(), 0);
    std::sort(perm_rows_ker.begin(), perm_rows_ker.end(),
              [&](const long& a, const long& b)->bool
              {
              return (pivind[a] < pivind[b]);
              } );

    kerbas.SetDims(k,m);
    for (long i = 0; i < k; ++i)
        for (long j = 0; j < m; ++j)
            kerbas[i][perm_shift[j]] = p_kerbas[perm_rows_ker[i]][j];
#ifdef MBASIS1_PROFILE
    t_perm2 = GetWallTime() - t_now;
#endif

#ifdef MBASIS1_PROFILE
    std::cout << "~~popov_mbasis1~~ input dimensions " << m << " x " << n << std::endl;
    std::cout << "  Initial permutations: " << t_perm1 << std::endl;
    std::cout << "  Computing kernel: " << t_ker << std::endl;
    std::cout << "  Computing pivot indices: " << t_pivind << std::endl;
    std::cout << "  Final permutations: " << t_perm2 << std::endl;
#endif

    return pivdeg;
}

/*------------------------------------------------------------*/
/* base case: order 1, i.e. working modulo X                  */
/* returning ordered weak Popov form                          */
/*------------------------------------------------------------*/

VecLong mbasis1(
                Mat<zz_p> & kerbas,
                const Mat<zz_p> & pmat,
                const VecLong & shift
               )
{
    long m = pmat.NumRows();
    long n = pmat.NumCols();
    // compute permutation which realizes stable sort of the shift
    // (i.e. sorts (shift[0],0)....(shift[len],len) lexicographically increasingly)
#ifdef MBASIS1_PROFILE
    double t_perm1,t_perm2,t_pivind,t_ker,t_now;
    t_now = GetWallTime();
#endif
    VecLong perm_shift(m);
    std::iota(perm_shift.begin(), perm_shift.end(), 0);
    stable_sort(perm_shift.begin(), perm_shift.end(),
                [&](const long& a, const long& b)->bool
                {
                return (shift[a] < shift[b]);
                } );

    // permute rows of pmat accordingly
    Mat<zz_p> mat;
    mat.SetDims(m,n);
    for (long i = 0; i < m; ++i)
        mat[i] = pmat[perm_shift[i]];
#ifdef MBASIS1_PROFILE
    t_perm1 = GetWallTime() - t_now;
#endif

    // find the permuted kernel basis in row echelon form
#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    Mat<zz_p> p_kerbas;
    kernel(p_kerbas,mat);
#ifdef MBASIS1_PROFILE
    t_ker = GetWallTime() - t_now;
#endif
    long k = p_kerbas.NumRows();
    if (k==0)
    {
        // kerbas is empty, same as p_kerbas
        kerbas.move(p_kerbas);
        return VecLong(m,1);
    }
    if (k==m)
    {
        // kerbas should be the identity, same as p_kerbas
        kerbas.move(p_kerbas);
        return VecLong(m,0);
    }

#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    // Compute pivot indices (pivot = rightmost nonzero entry)
    // Experiments show that:
    //   * kernel is expected to be of the form [ K | Id ]
    //   * in general it is a column-permutation of such a matrix
    // However note that a column-permutation is not sufficient for our needs
    // Another property: if pivots are in the expected location (diagonal of
    // rightmost square submatrix), then the corresponding column is the identity column.
    bool expected_pivots = true;
    VecLong p_pivind(k,m-1);
    for (long i = 0; i<k; ++i)
    {
        while (IsZero(p_kerbas[i][p_pivind[i]]))
            --p_pivind[i];
        if (p_pivind[i] != m-k+i)
            expected_pivots = false;
    }

    if (not expected_pivots)
    {
        // find whether p_pivind has pairwise distinct entries
        VecLong duplicate(p_pivind);
        std::sort(duplicate.begin(), duplicate.end());
        // if it has, then fine, the basis will not be Popov
        // but will be ordered weak Popov
        // (the goal of expected_pivots above was just to
        // avoid this call to sort in the most usual case)

        if (std::adjacent_find(duplicate.begin(),duplicate.end()) != duplicate.end())
        {
            // the kernel is not in a shape we can deduce the appbas from (some pivots collide)
            // --> let's compute its lower triangular row echelon form
            // (use kerbas as tmp space)
            kerbas.SetDims(k,m);
            kerbas.SetDims(k,m);
            for (long i = 0; i < k; ++i)
                for (long j = 0; j < m; ++j)
                    kerbas[i][j] = p_kerbas[i][m-1-j];
            image(kerbas, kerbas);
            // now column_permuted_ker is in upper triangular row echelon form
            for (long i = 0; i < k; ++i)
                for (long j = 0; j < m; ++j)
                    p_kerbas[i][j] = kerbas[k-i-1][m-1-j];
            // and now p_kerbas is the sought lower triangular row echelon kernel

            // compute the actual pivot indices
            for (long i = 0; i<k; ++i)
            {
                p_pivind[i] = m-1;
                while (IsZero(p_kerbas[i][p_pivind[i]]))
                    --p_pivind[i];
            }
        }
    }

    // now, up to row permutation, the kernel is in "lower triangular" row
    // echelon form (almost there: we want the non-permuted one)
    // permute everything back to original order:
    // prepare kernel permutation by permuting kernel pivot indices
    // pivot degrees corresponding to kernel pivot indices are 0, others are 1
    VecLong pivind(k);
    VecLong pivdeg(m,1);
    for (long i = 0; i < k; ++i)
    {
        pivind[i] = perm_shift[p_pivind[i]];
        pivdeg[pivind[i]] = 0;
    }

    VecLong perm_rows_ker(k);
    std::iota(perm_rows_ker.begin(), perm_rows_ker.end(), 0);
    std::sort(perm_rows_ker.begin(), perm_rows_ker.end(),
              [&](const long& a, const long& b)->bool
              {
              return (pivind[a] < pivind[b]);
              } );

    kerbas.SetDims(k,m);
    for (long i = 0; i < k; ++i)
        for (long j = 0; j < m; ++j)
            kerbas[i][perm_shift[j]] = p_kerbas[perm_rows_ker[i]][j];

#ifdef MBASIS1_PROFILE
    t_perm2 = GetWallTime() - t_now;
#endif

#ifdef MBASIS1_PROFILE
    std::cout << "~~mbasis1~~ input dimensions " << m << " x " << n << std::endl;
    std::cout << "  Initial permutations: " << t_perm1 << std::endl;
    std::cout << "  Computing kernel: " << t_ker << std::endl;
    std::cout << "  Computing pivot indices: " << t_pivind << std::endl;
    std::cout << "  Final permutations: " << t_perm2 << std::endl;
#endif

    return pivdeg;
}

/*------------------------------------------------------------*/
/* plain mbasis with polynomial matrices                      */
/*------------------------------------------------------------*/
void mbasis_plain(
                  Mat<zz_pX> & appbas,
                  const Mat<zz_pX> & pmat,
                  const long order,
                  VecLong & shift
                 )
{
    // initially, appbas is the identity matrix
    ident(appbas,pmat.NumRows());

    // `shift` will be updated along the algorithm,
    // to always be the shifted row degree of appbas
    // (initially this is exactly the input shift)

    long deg_pmat = deg(pmat);

    // will store the pivot degree at each call of mbasis1
    VecLong diff_pivdeg;

    // matrix to store the kernels in mbasis1 calls
    Mat<zz_p> kerbas;

    // matrix to store residuals, initially constant coeff of pmat
    Mat<zz_p> residual(coeff(pmat,0));

    // declare matrices
    Mat<zz_p> res_coeff,res_coeff1,res_coeff2; // will store coefficient matrices used to compute the residual
    Mat<zz_pX> kerapp; // will store constant-kernel * appbas

    for (long ord = 1; ord <= order; ++ord)
    {
        diff_pivdeg = mbasis1(kerbas,residual,shift);

        if (kerbas.NumRows()==0)
        {
            // computation is already finished: the final basis is X^(order-ord+1)*appbas
            appbas <<= (order-ord+1);
            // compute pivot degree, and return
            for (long i = 0; i < pmat.NumRows(); ++i)
                shift[i] += order-ord+1;
            return;
        }

        if (kerbas.NumRows()<residual.NumRows())
        {
            // I/ Update degrees:
            // new shifted row degree = old one + diff_pivdeg
            std::transform(shift.begin(), shift.end(), diff_pivdeg.begin(), shift.begin(), std::plus<long>());

            // II/ update approximant basis

            //--> submatrix of rows with diff_pivdeg==0 is replaced by kerbas*appbas
            //--> rows with diff_pivdeg=1 are simply multiplied by X
            mul(kerapp,kerbas,appbas);
            long row=0;
            for (long i = 0; i < appbas.NumRows(); ++i)
            {
                if (diff_pivdeg[i]==0)
                {
                    appbas[i].swap(kerapp[row]);
                    ++row;
                }
                else
                    LeftShiftRow(appbas,appbas,i,1);
            }

            // III/ compute next residual, if needed
            // this is coefficient of degree ord in appbas * pmat
            if (ord<order)
            {
                long deg_appbas = deg(appbas);
                clear(residual);
                for (long d = std::max<long>(0,ord-deg_pmat); d <= deg_appbas; ++d) // note that deg_appbas <= ord holds
                {
                    res_coeff1 = coeff(appbas,d);
                    res_coeff2 = coeff(pmat,ord-d);
                    mul(res_coeff, res_coeff1, res_coeff2);
                    add(residual, residual, res_coeff);
                }
            }
        }
        else
        {
            // Exceptional case:  kerbas.NumRows()==residual.NumRows())
            // --> residual is zero, kerbas is identity
            // --> approximant basis is already correct for this order, no need
            // to change it or to change pivdeg

            // Find degree of appbas
            long deg_appbas = deg(appbas);

            // compute next residual, if needed (if ord<order)
            // this is coefficient of degree ord in appbas * pmat
            if (ord<order)
            {
                // at this point, before the following loop, residual is zero
                for (long d = std::max<long>(0,ord-deg_pmat); d <= deg_appbas; ++d) // note that deg_appbas <= ord holds
                {
                    res_coeff1 = coeff(appbas,d);
                    res_coeff2 = coeff(pmat,ord-d);
                    mul(res_coeff, res_coeff1, res_coeff2);
                    add(residual, residual, res_coeff);
                }
            }
        }
    }
}

/*------------------------------------------------------------*/
/* mbasis, using vectors of matrices                          */
/*------------------------------------------------------------*/

// version with residual constant matrix computed at each iteration
void mbasis_rescomp(
                    Mat<zz_pX> & appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    VecLong & shift
                   )
{
#ifdef MBASIS_PROFILE
    double t_others=0.0,t_residual=0.0,t_appbas=0.0,t_kernel=0.0,t_now;
    t_now = GetWallTime();
#endif
    // A. General

    // A.1 dimensions of input matrix
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // A.2 store iota since it will be used at each iteration
    VecLong iota(m);
    std::iota(iota.begin(), iota.end(), 0);

    // B. Input representation; initialize output

    // B.1 convert input into vector of constant matrices (its "coefficients")
    const Vec<Mat<zz_p>> coeffs_pmat = conv(pmat,order);

    // B.2 vector of coefficients of output approximant basis
    Vec<Mat<zz_p>> coeffs_appbas;

    // B.3 initially, appbas is the identity matrix
    coeffs_appbas.SetLength(1);
    ident(coeffs_appbas[0], m);
    // degree of approximant basis, initially zero
    long deg_appbas = 0;
    // `shift` will always hold the shifted row degree of appbas
    // (initially, this is the input shift)

    // C. Residual matrix (m x n constant matrix, next coefficient
    // of appbas * pmat which we want to annihilate)

    // C.1 stores the residual, initially coeffs_pmat[0]
    Mat<zz_p> residuals(coeffs_pmat[0]);

    // C.2 temporary matrices used during the computation of residuals
    Mat<zz_p> res_coeff;

    // C.3 permuted residual, used as input to the kernel at the "base case"
    Mat<zz_p> p_residual(INIT_SIZE, m, n);

    // D. Base case (working modulo X, essentially amounts to finding the left
    // kernel of the permuted residual p_residual)

    // D.1 pivot indices in kernel basis (which is in row echelon form)
    // Note: length is probably overestimated (usually kernel has m-n rows),
    // but this avoids reallocating the right length at each iteration
    VecLong pivind(m-1);
    // Vector indicating if a given column index appears in this pivot index
    // i.e. is_pivind[pivind[i]] = true and others are false
    std::vector<bool> is_pivind(m, false);

    // D.2 permutation for the rows of the constant kernel
    VecLong perm_rows_ker;
    // pivot indices of row echelon form before permutation
    VecLong p_pivind(m-1);

    // D.3 permutation which stable-sorts the shift, used at the base case
    VecLong p_shift;

    // D.4 the constant kernel, and its permuted version
    Mat<zz_p> kerbas;
    Mat<zz_p> p_kerbas;

    // E. Updating appbas
    // stores the product "constant-kernel * coeffs_appbas[d]"
    Mat<zz_p> kerapp; 

#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif

    for (long ord = 1; ord <= order; ++ord)
    {
#ifdef MBASIS_PROFILE
        t_now = GetWallTime();
#endif
        // compute permutation which realizes stable sort of the row degree `shift`
        // --> we need to permute things, to take into account the "priority"
        // (i.e. "weights") indicated by the shift
        p_shift = iota;
        stable_sort(p_shift.begin(), p_shift.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (shift[a] < shift[b]);
                    } );

        // permute rows of the residual accordingly
        for (long i = 0; i < m; ++i)
            p_residual[i].swap(residuals[p_shift[i]]);
        // content of residual has been changed --> let's make it zero
        // (already the case if ord==1, since residual is then the former p_residual which was zero)
        if (ord>1)
            clear(residuals);
#ifdef MBASIS_PROFILE
        t_others += GetWallTime()-t_now;
        t_now = GetWallTime();
#endif
        // find the (permuted) left kernel basis, hopefully in row echelon form;
        kernel(p_kerbas,p_residual);
#ifdef MBASIS_PROFILE
        t_kernel += GetWallTime()-t_now;
#endif
        const long ker_dim = p_kerbas.NumRows();

        if (ker_dim==0)
        {
            // Exceptional case: the residual matrix has empty left kernel
            // --> no need to compute more: the final basis is X^(order-ord+1)*coeffs_appbas
            // TODO improve: simultaneous convert+shift!!
            appbas = conv(coeffs_appbas);
            appbas <<= (order-ord+1);
            for (long i = 0; i < m; ++i)
                shift[i] += order-ord+1;
            return;
        }

        else if (ker_dim==m)
        {
            // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
            // --> approximant basis is already correct for this order, no need to
            // change it or to change shift
            // --> we just need to compute the next residual
            // (unless ord == order, in which case the algorithm returns)
            // this "residual" is the coefficient of degree ord in appbas * pmat:
            // Note: at this point, residuals==0
            if (ord<order)
            {
#ifdef MBASIS_PROFILE
                t_now = GetWallTime();
#endif
                for (long d = std::max<long>(0,ord-coeffs_pmat.length()+1); d <= deg_appbas; ++d)
                {
                    mul(res_coeff, coeffs_appbas[d], coeffs_pmat[ord-d]);
                    add(residuals, residuals, res_coeff);
                }
#ifdef MBASIS_PROFILE
                t_residual += GetWallTime()-t_now;
#endif
            }
        }

        else
        {
            // here, we are in the "usual" case, where the left kernel of the
            // residual has no special shape

            // first, we permute everything back to original order

#ifdef MBASIS_PROFILE
            t_now = GetWallTime();
#endif
            // Compute pivots indices (pivot = rightmost nonzero entry)
            // Experiments show that:
            //   * kernel is expected to be of the form [ K | Id ]
            //   * in general it is a column-permutation of such a matrix
            // However note that a column-permutation is not sufficient for our needs
            // Another property: if pivots are in the expected location (diagonal of
            // rightmost square submatrix), then the corresponding column is the identity column.
            bool expected_pivots = true;
            for (long i = 0; i<ker_dim; ++i)
            {
                p_pivind[i] = m-1;
                while (IsZero(p_kerbas[i][p_pivind[i]]))
                    --p_pivind[i];
                if (p_pivind[i] != m-ker_dim+i)
                    expected_pivots = false;
            }

            if (not expected_pivots)
            {
                // find whether p_pivind has pairwise distinct entries
                // (use pivind as temp space)
                pivind = p_pivind;
                std::sort(pivind.begin(), pivind.end());
                // if pairwise distinct, then fine, the basis will not
                // be Popov but will be ordered weak Popov (the goal of
                // expected_pivots above was just to avoid this call to
                // sort in the most usual case)

                if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
                {
                    // the kernel is not in a shape we can deduce the appbas from (some pivots collide)
                    // --> let's compute its lower triangular row echelon form
                    // (use kerbas as temporary space)
                    kerbas.SetDims(ker_dim,m);
                    for (long i = 0; i < ker_dim; ++i)
                        for (long j = 0; j < m; ++j)
                            kerbas[i][j] = p_kerbas[i][m-1-j];
                    image(kerbas, kerbas);
                    // now column_permuted_ker is in upper triangular row echelon form
                    for (long i = 0; i < ker_dim; ++i)
                        for (long j = 0; j < m; ++j)
                            p_kerbas[i][j] = kerbas[ker_dim-i-1][m-1-j];
                    // and now p_kerbas is the sought lower triangular row echelon kernel

                    // compute the actual pivot indices
                    for (long i = 0; i<ker_dim; ++i)
                    {
                        p_pivind[i] = m-1;
                        while (IsZero(p_kerbas[i][p_pivind[i]]))
                            --p_pivind[i];
                    }
                }
            }


            // up to row permutation, the kernel is in "lower triangular" row
            // echelon form (almost there: we want the non-permuted one)
            // prepare kernel permutation by permuting kernel pivot indices;
            // also record which rows are pivot index in this kernel
            // (note that before this loop, is_pivind is filled with 'false')
            for (long i = 0; i < ker_dim; ++i)
            {
                pivind[i] = p_shift[p_pivind[i]];
                is_pivind[pivind[i]] = true;
            }

            // perm_rows_ker = [0 1 2 ... ker_dim-1]
            perm_rows_ker.resize(ker_dim);
            std::copy_n(iota.begin(), ker_dim, perm_rows_ker.begin());
            // permutation putting the pivot indices pivind in increasing order
            sort(perm_rows_ker.begin(), perm_rows_ker.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (pivind[a] < pivind[b]);
                    } );

            // permute rows and columns of kernel back to original order
            kerbas.SetDims(ker_dim,m);
            for (long i = 0; i < ker_dim; ++i)
                for (long j = 0; j < m; ++j)
                    kerbas[i][p_shift[j]] = p_kerbas[perm_rows_ker[i]][j];
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Now, update shifted row degree:
            // entries corresponding to kernel pivot indices are kept, others are +1
            // Also, deduce the degree of appbas
            bool deg_updated=false;
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                {
                    ++shift[i];
                    if (not deg_updated && not IsZero(coeffs_appbas[deg_appbas][i]))
                    { ++deg_appbas; deg_updated=true; }
                }

            // this new degree is either unchanged (== coeffs_appbas.length()-1),
            // or is the old one + 1 (== coeffs_appbas.length())
            if (deg_appbas==coeffs_appbas.length())
            {
                coeffs_appbas.SetLength(deg_appbas+1);
                coeffs_appbas[deg_appbas].SetDims(m, m);
            }
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Update approximant basis

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*coeffs_appbas (note: these rows currently have degree
            // at most deg_appbas)
            // TODO possible small improvement for uniform shift: these rows
            // have degree less than deg_appbas, in this case (and deg_appbas
            // is reached on the diagonal, among the pivot degrees)
            for (long d = 0; d <= deg_appbas; ++d)
            {
                mul(kerapp, kerbas, coeffs_appbas[d]);
                for (long i = 0; i < ker_dim; ++i)
                    coeffs_appbas[d][pivind[perm_rows_ker[i]]].swap(kerapp[i]);
            }

            // rows with !is_pivind are multiplied by X (note: these rows
            // currently have degree less than deg_appbas)
            for (long d = deg_appbas-1; d >= 0; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                        coeffs_appbas[d+1][i].swap(coeffs_appbas[d][i]);
            // Note: after this, the row coeffs_appbas[0][i] is zero
#ifdef MBASIS_PROFILE
            t_appbas += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Find next residual: coefficient of degree ord in appbas*pmat
            // (this is not necessary if ord==order, since in this case
            // we have finished: appbas*pmat is zero mod X^order)
            // Note: at this point, residuals==0
            if (ord<order)
            {
                long dmin=std::max<long>(0,ord-coeffs_pmat.length()+1);
                for (long d = dmin; d < deg_appbas+1; ++d) // we have deg_appbas <= ord
                {
                    mul(res_coeff, coeffs_appbas[d], coeffs_pmat[ord-d]);
                    add(residuals, residuals, res_coeff);
                }
#ifdef MBASIS_PROFILE
                t_residual += GetWallTime()-t_now;
                t_now = GetWallTime();
#endif
                // Restore is_pivind to all false, as it should be at the beginning of
                // the iteration
                for (long i = 0; i < ker_dim; ++i)
                    is_pivind[pivind[i]] = false;
#ifdef MBASIS_PROFILE
                t_others += GetWallTime()-t_now;
#endif
            }
        }
    }

#ifdef MBASIS_PROFILE
    t_now = GetWallTime();
#endif
    // Convert approximant basis to polynomial matrix representation
    appbas = conv(coeffs_appbas);
#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif
#ifdef MBASIS_PROFILE
    double t_total = t_residual + t_appbas + t_kernel + t_others;
    std::cout << "~~mbasis_rescomp~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_residual/t_total << "," << t_appbas/t_total << "," <<
    t_kernel/t_total << "," << t_others/t_total << std::endl;
#endif
}

// version with residual constant matrix computed at each iteration, multi-threaded
// TODO work in progress
// before more work, should be merged with updated rescomp above
// (among others, it is probably wrong for some non-generic-rank-profile inputs)
//VecLong mbasis_rescomp_multithread(
//                                   Mat<zz_pX> & appbas,
//                                   const Mat<zz_pX> & pmat,
//                                   const long order,
//                                   const VecLong & shift
//                                  )
//{
//#ifdef MBASIS_PROFILE
//    double t_others=0.0,t_residual=0.0,t_appbas=0.0,t_kernel=0.0,t_now;
//    t_now = GetWallTime();
//#endif
//    // A. General
//
//    // A.1 store context, useful if multiple threads
//    zz_pContext context;
//    context.save();
//    const long nthreads = AvailableThreads();
//
//    // A.2 dimensions of input matrix
//    const long m = pmat.NumRows();
//    const long n = pmat.NumCols();
//
//    // A.3 store iota since it will be used at each iteration
//    VecLong iota(m);
//    std::iota(iota.begin(), iota.end(), 0);
//
//    // B. Input representation; initialize output
//
//    // B.1 convert input into vector of constant matrices (its "coefficients")
//    const Vec<Mat<zz_p>> coeffs_pmat = conv(pmat,order);
//
//    // B.2 vector of coefficients of output approximant basis
//    Vec<Mat<zz_p>> coeffs_appbas;
//
//    // B.3 initially, appbas is the identity matrix
//    coeffs_appbas.SetLength(1);
//    ident(coeffs_appbas[0], m);
//    // degree of approximant basis, initially zero
//    long deg_appbas = 0;
//    // shifted row degree of appbas, initially equal to shift
//    VecLong rdeg(shift);
//    // pivot degree of appbas, initially zero
//    // (note that along this algorithm we have pivdeg+shift = rdeg, entrywise,
//    // since we will compute appbas in ordered weak Popov form)
//    VecLong pivdeg(m);
//
//    // C. Residual matrix (m x n constant matrix, next coefficient
//    // of appbas * pmat which we want to annihilate)
//
//    // C.1 stores the residual (residual is in residuals[0], yet
//    // during some parallel computations it is decomposed into some parts that
//    // are computed in residuals[1], ... residuals[nthreads] and only afterwards
//    // combined into residuals[0]
//    Vec<Mat<zz_p>> residuals;
//    residuals.SetLength(nthreads);
//    residuals[0] = coeffs_pmat[0];
//    for (long i = 0; i < nthreads; ++i)
//        residuals[i].SetDims(m, n);
//
//    // C.2 temporary matrices used during the computation of residuals[i]
//    Vec<Mat<zz_p>> res_coeff;
//    res_coeff.SetLength(nthreads);
//
//    // C.3 permuted residual, used as input to the kernel at the "base case"
//    Mat<zz_p> p_residual;
//    p_residual.SetDims(m, n);
//
//    // D. Base case (working modulo X, essentially amounts to finding the left
//    // kernel of the permuted residual p_residual)
//
//    // D.1 pivot indices in kernel basis (which is in row echelon form)
//    // Note: length is probably overestimated (usually kernel has m-n rows),
//    // but this avoids reallocating the right length at each iteration
//    VecLong pivind(m-1);
//    // Vector indicating if a given column index appears in this pivot index
//    // i.e. is_pivind[pivind[i]] = true and others are false
//    std::vector<bool> is_pivind(m, false);
//
//    // D.2 permutation for the rows of the constant kernel
//    VecLong perm_rows_ker;
//    // pivot indices of row echelon form before permutation
//    VecLong p_pivind(m-1);
//
//    // D.3 permutation which stable-sorts the shift, used at the base case
//    VecLong p_rdeg;
//
//    // D.4 the constant kernel, and its permuted version
//    Mat<zz_p> kerbas;
//    Mat<zz_p> p_kerbas;
//
//    // E. Updating appbas
//    // stores the product "constant-kernel * coeffs_appbas[d]"
//    Vec<Mat<zz_p>> kerapp; 
//    kerapp.SetLength(nthreads);
//
//#ifdef MBASIS_PROFILE
//    t_others += GetWallTime()-t_now;
//#endif
//
//    for (long ord = 1; ord <= order; ++ord)
//    {
//#ifdef MBASIS_PROFILE
//        t_now = GetWallTime();
//#endif
//        // compute permutation which realizes stable sort of rdeg
//        // --> we need to permute things, to take into account the "priority"
//        // (i.e. "weights") indicated by the shift; at this stage, the input
//        // shift is the shift-row degree 'rdeg' of appbas
//        p_rdeg = iota;
//        stable_sort(p_rdeg.begin(), p_rdeg.end(),
//                    [&](const long& a, const long& b)->bool
//                    {
//                    return (rdeg[a] < rdeg[b]);
//                    } );
//
//        // permute rows of the residual accordingly
//        for (long i = 0; i < m; ++i)
//            p_residual[i].swap(residuals[0][p_rdeg[i]]);
//        // content of residual has been changed --> let's make it zero
//        // (already the case if ord==1, since residual is then the former p_residual which was zero)
//        if (ord>1)
//            clear(residuals[0]);
//#ifdef MBASIS_PROFILE
//        t_others += GetWallTime()-t_now;
//        t_now = GetWallTime();
//#endif
//        // find the (permuted) left kernel basis in row echelon form;
//        // note it might not exactly be the usual row echelon form but a
//        // row-permuted version --> this is taken into account below in the
//        // computation of the pivot indices
//        kernel(p_kerbas,p_residual);
//#ifdef MBASIS_PROFILE
//        t_kernel += GetWallTime()-t_now;
//#endif
//        long ker_dim = p_kerbas.NumRows();
//
//        if (ker_dim==0)
//        {
//            // Exceptional case: the residual matrix has empty left kernel
//            // --> no need to compute more: the final basis is X^(order-ord+1)*coeffs_appbas
//            appbas = conv(coeffs_appbas);
//            appbas <<= (order-ord+1);
//            std::for_each(pivdeg.begin(), pivdeg.end(), [&order,&ord](long& a) { a+=order-ord+1; });
//            return pivdeg;
//        }
//
//        else if (ker_dim==m)
//        {
//            // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
//            // --> approximant basis is already correct for this order, no need to
//            // change it or to change pivdeg
//            // --> we just need to compute the next residual
//            // (unless ord == order, in which case the algorithm returns)
//            // this "residual" is the coefficient of degree ord in appbas * pmat:
//            if (ord<order)
//            {
//#ifdef MBASIS_PROFILE
//                t_now = GetWallTime();
//#endif
//                // TODO make this parallel as below
//                for (long d = std::max<long>(0,ord-coeffs_pmat.length()+1); d <= deg_appbas; ++d)
//                {
//                    mul(res_coeff[0], coeffs_appbas[d], coeffs_pmat[ord-d]);
//                    add(residuals[0], residuals[0], res_coeff[0]);
//                }
//#ifdef MBASIS_PROFILE
//                t_residual += GetWallTime()-t_now;
//#endif
//            }
//        }
//
//        else
//        {
//            // here, we are in the "usual" case, where the left kernel of the
//            // residual has no special shape
//
//            // first, we permute everything back to original order
//
//            // compute the (permuted) pivot indices
//#ifdef MBASIS_PROFILE
//            t_now = GetWallTime();
//#endif
//            for (long i = 0; i<ker_dim; ++i)
//            {
//                p_pivind[i] = m-1;
//                while (IsZero(p_kerbas[i][p_pivind[i]]))
//                    --p_pivind[i];
//            }
//
//            // prepare kernel permutation by permuting kernel pivot indices;
//            // also record which rows are pivot index in this kernel
//            // (note that before this loop, is_pivind is filled with 'false')
//            for (long i = 0; i < ker_dim; ++i)
//            {
//                pivind[i] = p_rdeg[p_pivind[i]];
//                is_pivind[pivind[i]] = true;
//            }
//
//            // perm_rows_ker = [0 1 2 ... ker_dim-1]
//            perm_rows_ker.resize(ker_dim);
//            std::copy_n(iota.begin(), ker_dim, perm_rows_ker.begin());
//            // permutation putting the pivot indices pivind in increasing order
//            sort(perm_rows_ker.begin(), perm_rows_ker.end(),
//                 [&](const long& a, const long& b)->bool
//                 {
//                 return (pivind[a] < pivind[b]);
//                 } );
//
//            // permute rows and columns of kernel back to original order
//            kerbas.SetDims(ker_dim,m);
//            for (long i = 0; i < ker_dim; ++i)
//                for (long j = 0; j < m; ++j)
//                    kerbas[i][p_rdeg[j]] = p_kerbas[perm_rows_ker[i]][j];
//#ifdef MBASIS_PROFILE
//            t_others += GetWallTime()-t_now;
//            t_now = GetWallTime();
//#endif
//
//            // Now, update shifted pivot degree and shifted row degree:
//            // entries corresponding to kernel pivot indices are kept, others are +1
//            for (long i = 0; i < m; ++i)
//                if (not is_pivind[i])
//                {
//                    ++rdeg[i];
//                    ++pivdeg[i];
//                }
//
//            // Deduce the degree of appbas; it is a property of this algorithm
//            // that deg(appbas) = max(pivot degree) (i.e. max(degree of
//            // diagonal entries); this does not hold in general for shifted
//            // ordered weak Popov approximant bases
//            deg_appbas = *std::max_element(pivdeg.begin(), pivdeg.end());
//            // this new degree is either unchanged (== coeffs_appbas.length()-1),
//            // or is the old one + 1 (== coeffs_appbas.length())
//            if (deg_appbas==coeffs_appbas.length())
//            {
//                coeffs_appbas.SetLength(deg_appbas+1);
//                coeffs_appbas[deg_appbas].SetDims(m, m);
//            }
//#ifdef MBASIS_PROFILE
//            t_others += GetWallTime()-t_now;
//            t_now = GetWallTime();
//#endif
//
//            // Update approximant basis
//
//            // Submatrix of rows corresponding to pivind are replaced by
//            // kerbas*coeffs_appbas
//            PartitionInfo pinfo(deg_appbas);
//            NTL_EXEC_INDEX(pinfo.NumIntervals(), index)
//            context.restore();
//            long first, last;
//            pinfo.interval(first, last, index);
//            for (long d = first; d < last; ++d)
//            {
//                kerapp[index] = kerbas * coeffs_appbas[d];
//                for (long i = 0; i < ker_dim; ++i)
//                    coeffs_appbas[d][pivind[perm_rows_ker[i]]].swap(kerapp[index][i]);
//            }
//            NTL_EXEC_INDEX_END
//
//            // rows with !is_pivind are multiplied by X (note: these rows have
//            // degree less than deg_appbas)
//            for (long d = deg_appbas-1; d >= 0; --d)
//                for (long i = 0; i < m; ++i)
//                    if (not is_pivind[i])
//                        coeffs_appbas[d+1][i].swap(coeffs_appbas[d][i]);
//            // Note: after this, the row coeffs_appbas[0][i] is zero
//#ifdef MBASIS_PROFILE
//            t_appbas += GetWallTime()-t_now;
//            t_now = GetWallTime();
//#endif
//            // Find next residual: coefficient of degree ord in appbas*pmat
//            // (this is not necessary if ord==order, since in this case
//            // we have finished: appbas*pmat is zero mod X^order)
//            if (ord<order)
//            {
//                long dmin=std::max<long>(0,ord-coeffs_pmat.length()+1);
//                PartitionInfo pinfo(deg_appbas-dmin+1);
//                for (long i = 0; i < residuals.length(); ++i)
//                    clear(residuals[i]);
//                NTL_EXEC_INDEX(pinfo.NumIntervals(), index)
//                context.restore();
//                long first, last;
//                pinfo.interval(first, last, index);
//                for (long d = first; d < last; ++d) // we have deg_appbas <= ord
//                {
//                    mul(res_coeff[index], coeffs_appbas[d+dmin], coeffs_pmat[ord-d-dmin]);
//                    add(residuals[index], residuals[index], res_coeff[index]);
//                }
//                NTL_EXEC_INDEX_END
//                for (long i = 1; i < residuals.length(); ++i)
//                    residuals[0] += residuals[i];
//#ifdef MBASIS_PROFILE
//                t_residual += GetWallTime()-t_now;
//                t_now = GetWallTime();
//#endif
//                // Restore is_pivind to all false, as it should be at the beginning of
//                // the iteration
//                for (long i = 0; i < ker_dim; ++i)
//                    is_pivind[pivind[i]] = false;
//#ifdef MBASIS_PROFILE
//                t_others += GetWallTime()-t_now;
//#endif
//            }
//        }
//    }
//
//#ifdef MBASIS_PROFILE
//    t_now = GetWallTime();
//#endif
//    // Convert approximant basis to polynomial matrix representation
//    appbas = conv(coeffs_appbas);
//#ifdef MBASIS_PROFILE
//    t_others += GetWallTime()-t_now;
//#endif
//#ifdef MBASIS_PROFILE
//    double t_total = t_residual + t_appbas + t_kernel + t_others;
//    std::cout << "~~mbasis_rescomp~~\t (residuals,basis,kernel,others): \t ";
//    std::cout << t_residual/t_total << "," << t_appbas/t_total << "," <<
//    t_kernel/t_total << "," << t_others/t_total << std::endl;
//#endif
//
//    return pivdeg;
//}

// version with full residual matrix continuously updated along the iterations
void mbasis_resupdate(
                      Mat<zz_pX> & appbas,
                      const Mat<zz_pX> & pmat,
                      const long order,
                      VecLong & shift
                     )
{
#ifdef MBASIS_PROFILE
    double t_others=0.0,t_residual=0.0,t_appbas=0.0,t_kernel=0.0,t_now;
    t_now = GetWallTime();
#endif
    // A. General

    // A.1 dimensions of input matrix
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // A.2 store iota since it will be used at each iteration
    VecLong iota(m);
    std::iota(iota.begin(), iota.end(), 0);

    // B. Initialize output

    // B.1 vector of coefficients of output approximant basis
    Vec<Mat<zz_p>> coeffs_appbas;

    // B.2 initially, appbas is the identity matrix
    coeffs_appbas.SetLength(1);
    ident(coeffs_appbas[0], m);
    // degree of approximant basis, initially zero
    long deg_appbas = 0;

    // C. Residual matrix (holds the next constant matrix we want to annihilate; in this version, also holds further matrices that are continuously updated)

    // C.1 convert input into vector of constant matrices (its "coefficients")
    // in this version "resupdate", this will be our residual
    Vec<Mat<zz_p>> residuals = conv(pmat,order);

    // C.2 temporary matrices used during the update of residuals: kernel * residual
    Mat<zz_p> kerres;

    // C.3 permuted residual, used as input to the kernel at the "base case"
    Mat<zz_p> p_residual;
    p_residual.SetDims(m, n);

    // D. Base case (working modulo X, essentially amounts to finding the left
    // kernel of the permuted residual p_residual)

    // D.1 pivot indices in kernel basis (which is in row echelon form)
    // Note: length is probably overestimated (usually kernel has m-n rows),
    // but this avoids reallocating the right length at each iteration
    VecLong pivind(m-1);
    // Vector indicating if a given column index appears in this pivot index
    // i.e. is_pivind[pivind[i]] = true and others are false
    std::vector<bool> is_pivind(m, false);

    // D.2 permutation for the rows of the constant kernel
    VecLong perm_rows_ker;
    // pivot indices of row echelon form before permutation
    VecLong p_pivind(m-1);

    // D.3 permutation which stable-sorts the shift, used at the base case
    VecLong p_shift;

    // D.4 the constant kernel, and its permuted version
    Mat<zz_p> kerbas;
    Mat<zz_p> p_kerbas;

    // E. Updating appbas
    // stores the product "constant-kernel * coeffs_appbas[d]"
    Mat<zz_p> kerapp; 

#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif

    for (long ord = 1; ord <= order; ++ord)
    {
#ifdef MBASIS_PROFILE
        t_now = GetWallTime();
#endif
        // compute permutation which realizes stable sort of shift
        // --> we need to permute things, to take into account the "priority"
        // (i.e. "weights") indicated by the shift
        p_shift = iota;
        stable_sort(p_shift.begin(), p_shift.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (shift[a] < shift[b]);
                    } );

        // permute rows of the residual accordingly
        for (long i = 0; i < m; ++i)
            p_residual[i] = residuals[ord-1][p_shift[i]];
        // Note: not using swap here, since we will need residuals[ord-1] for updating the residual
#ifdef MBASIS_PROFILE
        t_others += GetWallTime()-t_now;
        t_now = GetWallTime();
#endif
        // find the (permuted) left kernel basis in row echelon form;
        // note it might not exactly be the usual row echelon form but a
        // row-permuted version --> this is taken into account below in the
        // computation of the pivot indices
        kernel(p_kerbas,p_residual);
#ifdef MBASIS_PROFILE
        t_kernel += GetWallTime()-t_now;
#endif
        long ker_dim = p_kerbas.NumRows();

        if (ker_dim==0)
        {
            // Exceptional case: the residual matrix has empty left kernel
            // --> no need to compute more: the final basis is X^(order-ord+1)*coeffs_appbas
            // TODO improve: simultaneous convert+shift!!
            appbas = conv(coeffs_appbas);
            appbas <<= (order-ord+1);
            for (long i = 0; i < m; ++i)
                shift[i] += order-ord+1;
            return;
        }

        if (ker_dim < m)
        {
#ifdef MBASIS_PROFILE
            t_now = GetWallTime();
#endif
            // Compute pivots indices (pivot = rightmost nonzero entry)
            // Experiments show that:
            //   * kernel is expected to be of the form [ K | Id ]
            //   * in general it is a column-permutation of such a matrix
            // However note that a column-permutation is not sufficient for our needs
            // Another property: if pivots are in the expected location (diagonal of
            // rightmost square submatrix), then the corresponding column is the identity column.
            bool expected_pivots = true;
            for (long i = 0; i<ker_dim; ++i)
            {
                p_pivind[i] = m-1;
                while (IsZero(p_kerbas[i][p_pivind[i]]))
                    --p_pivind[i];
                if (p_pivind[i] != m-ker_dim+i)
                    expected_pivots = false;
            }

            if (not expected_pivots)
            {
                // find whether p_pivind has pairwise distinct entries
                // (use pivind as temp space)
                pivind = p_pivind;
                std::sort(pivind.begin(), pivind.end());
                // if pairwise distinct, then fine, the basis will not
                // be Popov but will be ordered weak Popov (the goal of
                // expected_pivots above was just to avoid this call to
                // sort in the most usual case)

                if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
                {
                    // the kernel is not in a shape we can deduce the appbas from (some pivots collide)
                    // --> let's compute its lower triangular row echelon form
                    // (use kerbas as temporary space)
                    kerbas.SetDims(ker_dim,m);
                    for (long i = 0; i < ker_dim; ++i)
                        for (long j = 0; j < m; ++j)
                            kerbas[i][j] = p_kerbas[i][m-1-j];
                    image(kerbas, kerbas);
                    // now column_permuted_ker is in upper triangular row echelon form
                    for (long i = 0; i < ker_dim; ++i)
                        for (long j = 0; j < m; ++j)
                            p_kerbas[i][j] = kerbas[ker_dim-i-1][m-1-j];
                    // and now p_kerbas is the sought lower triangular row echelon kernel

                    // compute the actual pivot indices
                    for (long i = 0; i<ker_dim; ++i)
                    {
                        p_pivind[i] = m-1;
                        while (IsZero(p_kerbas[i][p_pivind[i]]))
                            --p_pivind[i];
                    }
                }
            }


            // up to row permutation, the kernel is in "lower triangular" row
            // echelon form (almost there: we want the non-permuted one)
            // prepare kernel permutation by permuting kernel pivot indices;
            // also record which rows are pivot index in this kernel
            // (note that before this loop, is_pivind is filled with 'false')
            for (long i = 0; i < ker_dim; ++i)
            {
                pivind[i] = p_shift[p_pivind[i]];
                is_pivind[pivind[i]] = true;
            }

            // perm_rows_ker = [0 1 2 ... ker_dim-1]
            perm_rows_ker.resize(ker_dim);
            std::copy_n(iota.begin(), ker_dim, perm_rows_ker.begin());
            // permutation putting the pivot indices pivind in increasing order
            sort(perm_rows_ker.begin(), perm_rows_ker.end(),
                 [&](const long& a, const long& b)->bool
                 {
                 return (pivind[a] < pivind[b]);
                 } );

            // permute rows and columns of kernel back to original order
            kerbas.SetDims(ker_dim,m);
            for (long i = 0; i < ker_dim; ++i)
                for (long j = 0; j < m; ++j)
                    kerbas[i][p_shift[j]] = p_kerbas[perm_rows_ker[i]][j];
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Now, update shifted row degree:
            // entries corresponding to kernel pivot indices are kept, others are +1
            // Also, deduce the new degree of appbas
            bool deg_updated = false;
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                {
                    ++shift[i];
                    if (not deg_updated && not IsZero(coeffs_appbas[deg_appbas][i]))
                    { ++deg_appbas; deg_updated=true; }
                }

            // this new degree is either unchanged (== coeffs_appbas.length()-1),
            // or is the old one + 1 (== coeffs_appbas.length())
            if (deg_appbas==coeffs_appbas.length())
            {
                coeffs_appbas.SetLength(deg_appbas+1);
                coeffs_appbas[deg_appbas].SetDims(m, m);
            }
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Update approximant basis

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*coeffs_appbas (note: these rows have degree at most
            // deg_appbas)
            // TODO possible small improvement for uniform shift: these rows
            // have degree less than deg_appbas, in this case (and deg_appbas
            // is reached on the diagonal, among the pivot degrees)
            for (long d = 0; d <= deg_appbas; ++d)
            {
                mul(kerapp, kerbas, coeffs_appbas[d]);
                for (long i = 0; i < ker_dim; ++i)
                    coeffs_appbas[d][pivind[perm_rows_ker[i]]].swap(kerapp[i]);
            }

            // rows with !is_pivind are multiplied by X (note: these rows have
            // degree less than deg_appbas)
            for (long d = deg_appbas-1; d >= 0; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                        coeffs_appbas[d+1][i].swap(coeffs_appbas[d][i]);
            // Note: after this, the row coeffs_appbas[0][i] is zero
#ifdef MBASIS_PROFILE
            t_appbas += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Update residual[ord...order-1] similarly
            // We do not consider residual[ord-1], since it is known
            // that its update will be zero (and it will not be used
            // in further iterations on ord)

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*residuals
            for (long d = ord; d < order; ++d)
            {
                mul(kerres, kerbas, residuals[d]);
                for (long i = 0; i < ker_dim; ++i)
                    residuals[d][pivind[perm_rows_ker[i]]].swap(kerres[i]);
            }

            // rows with !is_pivind are multiplied by X
            for (long d = order-2; d >= ord-1; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                        residuals[d+1][i].swap(residuals[d][i]);
#ifdef MBASIS_PROFILE
            t_residual += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Restore is_pivind to all false, as it should be at the beginning of
            // the iteration
            for (long i = 0; i < ker_dim; ++i)
                is_pivind[pivind[i]] = false;
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
#endif
        }
        // if (ker_dim==m)
        // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
        // --> approximant basis is already correct for this order, no need to
        // change it or to change rdeg
        // --> we just advance to the next residual matrix
    }

#ifdef MBASIS_PROFILE
    t_now = GetWallTime();
#endif
    // Convert approximant basis to polynomial matrix representation
    conv(appbas, coeffs_appbas);
#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif
#ifdef MBASIS_PROFILE
    double t_total = t_residual + t_appbas + t_kernel + t_others;
    std::cout << "~~mbasis_resupdate~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_residual/t_total << "," << t_appbas/t_total << "," <<
    t_kernel/t_total << "," << t_others/t_total << std::endl;
#endif
}


/*------------------------------------------------------------*/
/* M-Basis returning Popov basis                              */
/*------------------------------------------------------------*/
void popov_mbasis(
                  Mat<zz_pX> &appbas,
                  const Mat<zz_pX> & pmat,
                  const long order,
                  VecLong & shift
                 )
{
    // compute first basis, copy the shift because we need it afterwards
    VecLong rdeg(shift);
    mbasis(appbas,pmat,order,rdeg);

    // shift for second call: negated pivot degree
    VecLong popov_shift(pmat.NumRows());
    std::transform(shift.begin(), shift.end(), rdeg.begin(),
                   popov_shift.begin(), std::minus<long>());

    // output shifted row degree
    shift=rdeg;

    // save `popov_shift` using `rdeg` as a buffer
    rdeg=popov_shift;

    // second call, basis shifted Popov up to constant transformation
    mbasis(appbas,pmat,order,popov_shift);

    // perform the constant transformation
    Mat<zz_p> lmat;
    row_leading_matrix(lmat, appbas, rdeg);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas);
}




/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis                                */
/*------------------------------------------------------------*/
void pmbasis(
             Mat<zz_pX> & appbas,
             const Mat<zz_pX> & pmat,
             const long order,
             VecLong & shift
            )
{
    if (order <= 32) // TODO thresholds to be determined
    {
        mbasis(appbas,pmat,order,shift);
        return;
    }

    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call

    // first recursive call, with 'pmat' and 'shift'
    Mat<zz_pX> trunc_pmat; // truncated pmat for first call
    // TODO is this trunc necessary here? it seems to be just for avoiding
    // big copies at the base mbasis (right?) (probably just for resupdate?) so
    // it could be done directly inside mbasis when defining the residual?
    trunc(trunc_pmat,pmat,order1);
    pmbasis(appbas,trunc_pmat,order1,shift);

    // shift is now the shifted row degree of appbas,
    // which is the shift for second call

    // TODO remove once middle_product has been fixed? (currently that's why
    // it's here)
    // TODO require this for the input?
    trunc(trunc_pmat, pmat, order);

    // residual = (appbas * pmat * X^-order1) mod X^order2
    Mat<zz_pX> residual; // for the residual
    middle_product(residual, appbas, trunc_pmat, order1, order2-1);

    // second recursive call, with 'residual' and 'rdeg'
    Mat<zz_pX> appbas2; // basis for second call
    pmbasis(appbas2,residual,order2,shift);

    // final basis = appbas2 * appbas
    multiply(appbas,appbas2,appbas);
}

/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis returning Popov                */
/*------------------------------------------------------------*/
void popov_pmbasis(
                   Mat<zz_pX> &appbas,
                   const Mat<zz_pX> & pmat,
                   const long order,
                   VecLong & shift
                  )
{
    // compute first basis, copy the shift because we need it afterwards
    VecLong rdeg(shift);
    pmbasis(appbas,pmat,order,rdeg);

    // shift for second call: negated pivot degree
    VecLong popov_shift(pmat.NumRows());
    std::transform(shift.begin(), shift.end(), rdeg.begin(),
                   popov_shift.begin(), std::minus<long>());

    // output shifted row degree
    shift=rdeg;

    // save `popov_shift` using `rdeg` as a buffer
    rdeg=popov_shift;

    // second call, basis shifted Popov up to constant transformation
    pmbasis(appbas,pmat,order,popov_shift);

    // perform the constant transformation
    Mat<zz_p> lmat;
    row_leading_matrix(lmat, appbas, rdeg);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas);
}

void pmbasis_2x1(
                 zz_pX & p00,
                 zz_pX & p01,
                 zz_pX & p10,
                 zz_pX & p11,
                 const zz_pX & f0,
                 const zz_pX & f1,
                 long order,
                 long & s0,
                 long & s1,
                 long threshold
                )
{
    if (order <= threshold) // TODO thresholds to be determined
    {
        appbas_iterative_2x1(p00,p01,p10,p11,f0,f1,order,s0,s1);
        return;
    }

    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call

    double t=0;
    if (order>200000)
        t=GetWallTime();
    // first recursive call, with 'pmat' and 'shift'
    pmbasis_2x1(p00,p01,p10,p11,trunc(f0,order1),trunc(f1,order1),order1,s0,s1,threshold);
    if (order>200000)
    {
        t=GetWallTime()-t;
        std::cout << order << " -- 1st rec call: " << t << std::endl;
    }

    // (s0,s1) is now the shifted row degree of appbas,
    // which is the shift for second call

    // residual = (appbas * pmat * X^-order1) mod X^order2
    ////version 1:
    //zz_pX tmp,r0,r1;
    //middle_product(r0,p00,f0,order1,order2-1);
    //middle_product(r1,p01,f1,order1,order2-1); // using r1 as temporary
    //r0 += r1;
    //middle_product(r1,p10,f0,order1,order2-1);
    //middle_product(tmp,p11,f1,order1,order2-1); // using tf0 as temporary
    //r1 += tmp;
    ////version 2:
    //zz_pX r0,r1;
    //RightShift(r0, MulTrunc(p00, f0, order) + MulTrunc(p01, f1, order), order1);
    //RightShift(r1, MulTrunc(p10, f0, order) + MulTrunc(p11, f1, order), order1);
    ////version3:
    if (order>200000)
        t=GetWallTime();
    Mat<zz_pX> res, tmp;
    Mat<zz_pX> appbas,appbas2;
    appbas.SetDims(2,2);
    appbas[0][0] = p00;
    appbas[0][1] = p01;
    appbas[1][0] = p10;
    appbas[1][1] = p11;
    tmp.SetDims(2,1);
    tmp[0][0] = trunc(f0,order);
    tmp[1][0] = trunc(f1,order);
    RightShift(tmp, tmp, order1-deg(appbas));
    middle_product_evaluate_FFT_new(res, appbas, tmp, deg(appbas), order2-1);
    if (order>200000)
    {
        t=GetWallTime()-t;
        std::cout << order << " -- residual: " << t << std::endl;
    }

    // second recursive call, with the 'residual' (r1,r2) and updated (s0,s1)
    if (order>200000)
        t=GetWallTime();
    zz_pX q00,q01,q10,q11; // basis for second call
    pmbasis_2x1(q00,q01,q10,q11,res[0][0],res[1][0],order2,s0,s1,threshold);
    if (order>200000)
    {
        t=GetWallTime()-t;
        std::cout << order << " -- 2nd rec call: " << t << std::endl;
    }

    // second recursive call, with the 'residual' (r1,r2) and updated (s0,s1)
    //zz_pX q00,q01,q10,q11; // basis for second call
    //pmbasis_2x1(q00,q01,q10,q11,r0,r1,order2,s0,s1,threshold);

    // final basis = product of the two
    //Mat<zz_pX> appbas,appbas2;
    //appbas.SetDims(2,2);
    //appbas[0][0] = p00;
    //appbas[0][1] = p01;
    //appbas[1][0] = p10;
    //appbas[1][1] = p11;
    if (order>200000)
        t=GetWallTime();
    appbas2.SetDims(2,2);
    appbas2[0][0] = q00;
    appbas2[0][1] = q01;
    appbas2[1][0] = q10;
    appbas2[1][1] = q11;
    multiply(appbas,appbas2,appbas);
    p00 = appbas[0][0];
    p01 = appbas[0][1];
    p10 = appbas[1][0];
    p11 = appbas[1][1];
    if (order>200000)
    {
        t=GetWallTime()-t;
        std::cout << order << " -- multiply bases: " << t << std::endl;
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
