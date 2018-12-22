#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota
#include <NTL/BasicThreadPool.h>

#include "mat_lzz_pX_interpolant.h"
#include "mat_lzz_pX_approximant.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* interpolant basis verification                             */
/* TODO incomplete                                            */
/*------------------------------------------------------------*/
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & evals,
                          const Vec<zz_p> & pts,
                          const VecLong & shift,
                          const PolMatForm & form,
                          const bool randomized
                         )
{
    if (randomized)
        throw std::logic_error("==is_interpolant_basis== Fast randomized interpolant basis verification not implemented yet");

    //std::cout << "==is_interpolant_basis== WARNING: not fully implemented: not checking generation" << std::endl;

    // test that appbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,intbas,shift))
        return false;

    // test that the matrix consists of interpolants
    Vec<Mat<zz_p>> intbas_vec;
    conv(intbas_vec,intbas);
    long deg = intbas_vec.length()-1;
    for (long pt=0; pt<pts.length(); ++pt)
    {
        // evaluate with Horner
        Mat<zz_p> ev = intbas_vec[deg];
        for (long d = deg-1; d >= 0; --d)
        {
            ev = pts[pt] * ev;
            add(ev, ev, intbas_vec[d]);
        }
        mul(ev, ev, evals[pt]);
        if (not IsZero(ev))
            return false;
    }
    return true;
}

/*------------------------------------------------------------*/
/* M-Basis for interpolation points                           */
/*------------------------------------------------------------*/

// plain algorithm
VecLong mbasis(
               Mat<zz_pX> & intbas,
               const Vec<Mat<zz_p>> & evals,
               const Vec<zz_p> & pts,
               const VecLong & shift
              )
{
    const long nrows = evals[0].NumRows();

    // interpolation basis represented as a vector (polynomial) of constant matrices
    Vec<Mat<zz_p>> coeffs_intbas;

    // initially, intbas is the identity matrix
    coeffs_intbas.SetLength(1);
    ident(coeffs_intbas[0], nrows);

    // holds the current shifted row degree of coeffs_intbas
    // initially, this is exactly shift
    VecLong rdeg(shift);

    // holds the current pivot degree of coeffs_intbas
    // initially tuple of zeroes
    // (note that at all times pivdeg+shift = rdeg entrywise)
    VecLong pivdeg(nrows);

    // will store the pivot degree at each call of mbasis1
    VecLong diff_pivdeg;

    // matrix to store the kernels in mbasis1 calls
    Mat<zz_p> kerbas;
    // matrix to store residuals, initially constant coeff of evals
    Mat<zz_p> residual( evals[0] );
    Mat<zz_p> intbas_eval;
    intbas_eval.SetDims(nrows,nrows);

    // declare matrices
    Mat<zz_p> res_coeff; // will store coefficient matrices used to compute the residual
    Mat<zz_p> kerapp; // will store constant-kernel * coeffs_intbas[d]

    for (long pt = 1; pt <= pts.length(); ++pt)
    {
        // call MBasis1 to retrieve kernel and pivdeg
        diff_pivdeg = popov_mbasis1(kerbas,residual,rdeg);

        if (kerbas.NumRows()==0)
        {
            // the basis is simply multiplied by X-pt
            // new degree is old degree + 1, that is, the length
            long d = coeffs_intbas.length();
            coeffs_intbas.SetLength(d+1);
            // note that d>0
            coeffs_intbas[d] = coeffs_intbas[d-1];
            --d;
            for (; d > 0; --d)
            {
                mul(coeffs_intbas[d], coeffs_intbas[d], -pts[pt-1]);
                add(coeffs_intbas[d], coeffs_intbas[d], coeffs_intbas[d-1]);
            }
            mul(coeffs_intbas[0], coeffs_intbas[0], -pts[pt-1]);
            // update pivdeg accordingly
            std::for_each(pivdeg.begin(), pivdeg.end(), [](long& a) { ++a; });

            // TODO update residual??
        }

        // kerbas.NumRows()==residual.NumRows() --> interpolant basis is already
        // correct for this point, just go to the next

        else if (kerbas.NumRows()<residual.NumRows())
        {
            // I/ Update degrees:
            // new shifted row degree = old rdeg + diff_pivdeg
            std::transform(rdeg.begin(), rdeg.end(), diff_pivdeg.begin(), rdeg.begin(), std::plus<long>());
            // new pivot degree = old pivot_degree + diff_pivdeg
            std::transform(pivdeg.begin(), pivdeg.end(), diff_pivdeg.begin(), pivdeg.begin(), std::plus<long>());
            // deduce degree of coeffs_intbas; note that it is a property of this algorithm
            // that deg(coeffs_intbas) = max(pivot degree) (i.e. max(degree of diagonal
            // entries); this does not hold in general for ordered weak Popov forms
            long deg_coeffs_intbas = *std::max_element(pivdeg.begin(), pivdeg.end());

            // II/ update approximant basis

            // submatrix of rows with diff_pivdeg==0 is replaced by kerbas*coeffs_intbas
            // while rows with diff_pivdeg=1 are multiplied by X-pt
            // --> the loop goes downwards, so that we can do both in the same iteration
            long row;
            long d = deg_coeffs_intbas;
            // separate treatment of highest degree terms, if degree has increased
            if (d>=coeffs_intbas.length()) // then, it is equal
            {
                coeffs_intbas.SetLength(d+1); // increases length by 1
                coeffs_intbas[d].SetDims(nrows,nrows);
                for (long i = 0; i < nrows; ++i)
                    if (diff_pivdeg[i]==1)
                        coeffs_intbas[d][i] = coeffs_intbas[d-1][i];
                --d;
            }

            // normal treatment of other terms
            for (; d >= 0; --d)
            {
                mul(kerapp,kerbas,coeffs_intbas[d]);
                row=0;
                for (long i = 0; i < nrows; ++i)
                {
                    if (diff_pivdeg[i]==0)
                    {
                        coeffs_intbas[d][i] = kerapp[row];
                        ++row;
                    }
                    else  // diff_pivdeg[i]==1 --> multiply by X-pt
                    {
                        coeffs_intbas[d][i] *= -pts[pt-1];
                        if (d>0)
                            coeffs_intbas[d][i] += coeffs_intbas[d-1][i];
                    }
                }
            }

            // III/ compute next residual, if needed
            // it is coeffs_intbas(pts[pt]) * evals[pt]
            // evaluate with Horner method
            if (pt<pts.length())
            {
                intbas_eval = coeffs_intbas[deg_coeffs_intbas];
                d = deg_coeffs_intbas-1;
                for (; d >= 0; --d)
                {
                    intbas_eval = pts[pt] * intbas_eval;
                    add(intbas_eval, intbas_eval, coeffs_intbas[d]);
                }
                mul(residual,intbas_eval,evals[pt]);
            }
        }
    }

    conv(intbas,coeffs_intbas);
    return pivdeg;
}


// version with residual constant matrix computed at each iteration
VecLong mbasis_rescomp(
                       Mat<zz_pX> & intbas,
                       const Vec<Mat<zz_p>> & evals,
                       const Vec<zz_p> & pts,
                       const VecLong & shift
                      )
{
#ifdef MBASIS_PROFILE
    double t_others=0.0,t_residual=0.0,t_intbas=0.0,t_kernel=0.0,t_now;
    t_now = GetWallTime();
#endif
    // A. General

    // A.1 dimensions of input matrix
    const long m = evals[0].NumRows();
    const long n = evals[0].NumCols();

    // A.2 order of interpolation (number of points)
    const long order = evals.length();

    // A.3 store iota since it will be used at each iteration
    VecLong iota(m);
    std::iota(iota.begin(), iota.end(), 0);

    // B. Input representation; initialize output

    // B.1 vector of coefficients of output approximant basis
    Vec<Mat<zz_p>> coeffs_intbas;

    // B.2 initially, intbas is the identity matrix
    coeffs_intbas.SetLength(1);
    ident(coeffs_intbas[0], m);
    // degree of approximant basis, initially zero
    long deg_intbas = 0;
    // shifted row degree of intbas, initially equal to shift
    VecLong rdeg(shift);

    // C. Residual matrix (m x n constant matrix, next coefficient
    // of intbas * pmat which we want to annihilate)

    // C.1 stores the residual, initially evals[0]
    Mat<zz_p> residuals(evals[0]);

    // C.2 temporary matrices used during the computation of residuals
    Mat<zz_p> res_coeff;

    // C.3 permuted residual, used as input to the kernel at the "base case"
    Mat<zz_p> p_residual(INIT_SIZE, m, n);

    // D. Base case (essentially amounts to finding the left kernel of the
    // permuted residual p_residual)

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
    VecLong p_rdeg;

    // D.4 the constant kernel, and its permuted version
    Mat<zz_p> kerbas;
    Mat<zz_p> p_kerbas;

    // E. Updating intbas
    // stores the product "constant-kernel * coeffs_intbas[d]"
    Mat<zz_p> kerint; 

#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif

    for (long ord = 0; ord < order; ++ord)
    {
#ifdef MBASIS_PROFILE
        t_now = GetWallTime();
#endif
        // compute permutation which realizes stable sort of rdeg
        // --> we need to permute things, to take into account the "priority"
        // (i.e. "weights") indicated by the shift; at this stage, the input
        // shift is the shift-row degree 'rdeg' of intbas
        p_rdeg = iota;
        stable_sort(p_rdeg.begin(), p_rdeg.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (rdeg[a] < rdeg[b]);
                    } );

        // permute rows of the residual accordingly
        for (long i = 0; i < m; ++i)
            p_residual[i].swap(residuals[p_rdeg[i]]);
        // content of residual has been changed --> let's make it zero
        // (already the case if ord==0, since residual is then the former p_residual which was zero)
        if (ord>0)
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
            // --> left-multiply intbas by (x-pt[ord])
            // --> if the next point is different from pt[ord], the residual
            // becomes the next matrix in evals (no need to multiply all
            // entries by the same nonzero constant pt[ord]-pt[ord+1]); while
            // the next point is equal, advance until it is not equal

            // update intbas
            ++deg_intbas; // note that this degree is now > 0
            coeffs_intbas.SetLength(deg_intbas+1);
            coeffs_intbas[deg_intbas] = coeffs_intbas[deg_intbas-1];
            const zz_p pt(-pts[ord]);
            for (long d=deg_intbas-1; d > 0; --d)
            {
                mul(coeffs_intbas[d], coeffs_intbas[d], pt);
                add(coeffs_intbas[d], coeffs_intbas[d], coeffs_intbas[d-1]);
            }
            mul(coeffs_intbas[0], coeffs_intbas[0], pt);

            // update rdeg accordingly
            std::for_each(rdeg.begin(), rdeg.end(), [](long& a) { ++a; });

            // update residual, if not the last iteration
            while (ord<order-1 && pts[ord]==pts[ord+1])
                ++ord;

            // if we have found a different point, take the next residual
            if (ord<order-1)
                residuals = evals[ord+1];
        }

        else if (ker_dim==m)
        {
            // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
            // --> approximant basis is already correct for this point, no need to
            // change it or to change rdeg
            // --> we just need to compute the next residual
            // (unless we are at the last iteration, in which case the algorithm returns)
            // this "residual" is the coefficient of degree ord in intbas * pmat:
            if (ord<order)
            {
#ifdef MBASIS_PROFILE
                t_now = GetWallTime();
#endif
                for (long d = std::max<long>(0,ord-coeffs_pmat.length()+1); d <= deg_intbas; ++d)
                {
                    mul(res_coeff, coeffs_intbas[d], coeffs_pmat[ord-d]);
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
                    // the kernel is not in a shape we can deduce the intbas from (some pivots collide)
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
                pivind[i] = p_rdeg[p_pivind[i]];
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
                    kerbas[i][p_rdeg[j]] = p_kerbas[perm_rows_ker[i]][j];
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Now, update shifted row degree:
            // entries corresponding to kernel pivot indices are kept, others are +1
            // Also, deduce the degree of intbas
            bool deg_updated=false;
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                {
                    ++rdeg[i];
                    if (not deg_updated && not IsZero(coeffs_intbas[deg_intbas][i]))
                    { ++deg_intbas; deg_updated=true; }
                }

            // this new degree is either unchanged (== coeffs_intbas.length()-1),
            // or is the old one + 1 (== coeffs_intbas.length())
            if (deg_intbas==coeffs_intbas.length())
            {
                coeffs_intbas.SetLength(deg_intbas+1);
                coeffs_intbas[deg_intbas].SetDims(m, m);
            }
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Update approximant basis

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*coeffs_intbas (note: these rows currently have degree
            // at most deg_intbas)
            // TODO possible small improvement for uniform shift: these rows
            // have degree less than deg_intbas, in this case (and deg_intbas
            // is reached on the diagonal, among the pivot degrees)
            for (long d = 0; d <= deg_intbas; ++d)
            {
                kerint = kerbas * coeffs_intbas[d];
                for (long i = 0; i < ker_dim; ++i)
                    coeffs_intbas[d][pivind[perm_rows_ker[i]]].swap(kerint[i]);
            }

            // rows with !is_pivind are multiplied by X (note: these rows
            // currently have degree less than deg_intbas)
            for (long d = deg_intbas-1; d >= 0; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                        coeffs_intbas[d+1][i].swap(coeffs_intbas[d][i]);
            // Note: after this, the row coeffs_intbas[0][i] is zero
#ifdef MBASIS_PROFILE
            t_intbas += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Find next residual: coefficient of degree ord in intbas*pmat
            // (this is not necessary if ord==order, since in this case
            // we have finished: intbas*pmat is zero mod X^order)
            if (ord<order)
            {
                long dmin=std::max<long>(0,ord-coeffs_pmat.length()+1);
                clear(residuals);
                for (long d = dmin; d < deg_intbas+1; ++d) // we have deg_intbas <= ord
                {
                    mul(res_coeff, coeffs_intbas[d], coeffs_pmat[ord-d]);
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
    intbas = conv(coeffs_intbas);
#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif
#ifdef MBASIS_PROFILE
    double t_total = t_residual + t_intbas + t_kernel + t_others;
    std::cout << "~~mbasis_rescomp~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_residual/t_total << "," << t_intbas/t_total << "," <<
    t_kernel/t_total << "," << t_others/t_total << std::endl;
#endif

    // deduce pivot degree
    for (long i = 0; i < m; ++i)
        rdeg[i] -= shift[i];
    return rdeg;
}


VecLong pmbasis_geometric(
                          Mat<zz_pX> & intbas,
                          const Mat<zz_pX> & pmat,
                          const zz_p & r,
                          const long order,
                          const VecLong & shift,
                          Vec<Mat<zz_p>> &evals,
                          Vec<zz_p> &pts
                         )
{
    VecLong rdeg;
    rdeg.resize(pmat.NumRows());
    row_degree(rdeg, pmat);
    long max_rowdeg = rdeg[0];
    for (auto i : rdeg)
        if (i > max_rowdeg) max_rowdeg = i;
    max_rowdeg = max(order, max_rowdeg);
    zz_pX_Multipoint_Geometric eval(r, max_rowdeg+1);

    // set up pts
    zz_pX x;
    SetCoeff(x, 1, 1);
    Vec<zz_p> pts2;
    eval.evaluate(pts2, x); // just gets powers of r

    // set up evaluations of pmat
    evals.SetLength(order);
    pts.SetLength(order);
    for (long d = 0; d < order; d++)
    {
        pts[d] = pts2[d];
        evals[d].SetDims(pmat.NumRows(), pmat.NumCols());
    }
    for (long r = 0; r < pmat.NumRows(); r++)
    {
        for (long c = 0; c < pmat.NumCols(); c++)
        {
            Vec<zz_p> vals;
            eval.evaluate(vals, pmat[r][c]);
            for (long d = 0; d < order; d++)
                evals[d][r][c] = vals[d];
        }
    }

    return pmbasis_geometric(intbas, evals, pts, r, shift);
}

VecLong pmbasis_geometric(
                          Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & evals,
                          const Vec<zz_p> & pts,
                          const zz_p & r,
                          const VecLong & shift
                         )
{
    const long order = pts.length();
    zz_pContext context;

    if (order <= 32)
        return mbasis(intbas, evals, pts, shift);

    VecLong pivdeg; // pivot degree, first call
    VecLong pivdeg2; // pivot degree, second call
    VecLong rdeg(evals[0].NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> intbas2; // basis for second call

    // first recursive call
    Vec<zz_p>  pts1;
    pts1.SetLength(order1);
    for (long i = 0; i < order1; i++)
        pts1[i] = pts[i];
    pivdeg = pmbasis_geometric(intbas, evals, pts1, r, shift);

    long max_pivdeg = pivdeg[0];
    for (auto i : pivdeg)
        if (max_pivdeg < i) max_pivdeg = i;

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

    // get the product of evaluations intbas(x_i) * pmat(x_i)
    // for the second half of the points
    Vec<zz_p> pts2;
    pts2.SetLength(order2);
    for (long i=0; i<order2; ++i)
        pts2[i] = pts[order1+i];

    // geometric progression of r, starting at pts2[0]
    // first ensure we have enough points for the degree
    max_pivdeg = max(order2, max_pivdeg); 
    zz_pX_Multipoint_Geometric ev(r,pts2[0], max_pivdeg+1);

    Vec<Mat<zz_p>> evals2;
    evals2.SetLength(order2);
    for (long i = 0; i < order2; i++)
        evals2[i].SetDims(intbas.NumRows(),intbas.NumCols());

    context.save();  

    // evaluate and store
    NTL_EXEC_RANGE(intbas.NumRows(),first,last)

    context.restore();    
    for (long r = first; r < last; r++)
        for (long c = 0; c < intbas.NumCols(); c++)
        {
            Vec<zz_p> val;
            ev.evaluate(val, intbas[r][c]);

            for (long i = 0; i < order2; i++)
                evals2[i][r][c] = val[i];

        }
    NTL_EXEC_RANGE_END

    context.save();  

    // multiply and store    
    NTL_EXEC_RANGE(order2,first,last)
    context.restore();    

    for (long i = first; i < last; i++)
    {
        evals2[i] = evals2[i] * evals[order1+i];
    }
    NTL_EXEC_RANGE_END

    // second recursive call
    pivdeg2 = pmbasis_geometric(intbas2, evals2, pts2, r, rdeg);

    multiply(intbas,intbas2,intbas);

    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;    
}

VecLong pmbasis(
                Mat<zz_pX> & intbas,
                const Vec<Mat<zz_p>> & evals,
                const Vec<zz_p> & pts,
                const VecLong & shift
               )
{

    const long order = pts.length();
    zz_pContext context;

    if (order <= 32)
        return mbasis(intbas, evals, pts, shift);

    VecLong pivdeg; // pivot degree, first call
    VecLong pivdeg2; // pivot degree, second call
    VecLong rdeg(evals[0].NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> intbas2; // basis for second call

    // first recursive call
    Vec<zz_p>  pts1;
    pts1.SetLength(order1);
    for (long i = 0; i < order1; i++)
        pts1[i] = pts[i];
    pivdeg = pmbasis(intbas, evals, pts1, shift);

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

    // get the product of evaluations intbas(x_i) * pmat(x_i)
    // for the second half of the points
    Vec<zz_p> pts2;
    pts2.SetLength(order2);
    for (long i=0; i<order2; ++i)
        pts2[i] = pts[order1+i];

    zz_pX_Multipoint_General ev(pts2);

    Vec<Mat<zz_p>> evals2;
    evals2.SetLength(order2);
    for (long i = 0; i < order2; i++)
        evals2[i].SetDims(intbas.NumRows(),intbas.NumCols());

    context.save();  

    // evaluate and store
    NTL_EXEC_RANGE(intbas.NumRows(),first,last)

    context.restore();    
    for (long r = first; r < last; r++)
        for (long c = 0; c < intbas.NumCols(); c++)
        {
            Vec<zz_p> val;
            ev.evaluate(val, intbas[r][c]);
            for (long i = 0; i < order2; i++)
                evals2[i][r][c] = val[i];

        }
    NTL_EXEC_RANGE_END

    context.save();  

    // multiply and store    
    NTL_EXEC_RANGE(order2,first,last)
    context.restore();    

    for (long i = first; i < last; i++)
    {
        evals2[i] = evals2[i] * evals[order1+i];
    }
    NTL_EXEC_RANGE_END

    // second recursive call
    pivdeg2 = pmbasis(intbas2, evals2, pts2, rdeg);

    multiply(intbas,intbas2,intbas);

    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
