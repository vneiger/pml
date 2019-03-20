#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota
//#include <NTL/BasicThreadPool.h>

#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_interpolant.h"

//#define MBASIS_PROFILE

NTL_CLIENT

/*------------------------------------------------------------*/
/* interpolant basis verification                             */
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

    // test that intbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,intbas,shift))
        return false;

    // test that the matrix consists of interpolants
    zz_pX_Multipoint_General ev(pts);
    Vec<Mat<zz_p>> intbas_evals;
    ev.evaluate_matrix(intbas_evals, intbas);
    Mat<zz_p> res;
    for (long pt=0; pt<pts.length(); ++pt)
    {
        mul(res, intbas_evals[pt], evals[pt]);
        if (not IsZero(res))
            return false;
    }

    // test that intbas generates all interpolant, or equivalently,
    // that intbas has the right degree of determinant
    // --> since intbas is reduced, we can find the degree of determinant as
    // the sum of shifted row degrees minus the sum of shifts

    // compute degree of determinant of intbas
    VecLong rdeg;
    row_degree(rdeg, intbas, shift);
    long actual_degdet = std::accumulate(rdeg.begin(), rdeg.end(), 0) - std::accumulate(shift.begin(), shift.end(), 0);

    // compute sum of ranks of input matrices: this should be the
    // degree of determinant of intbas
    long target_degdet=0;
    for (long pt = 0; pt < pts.length(); ++pt)
    {
        Mat<zz_p> evals_pt(evals[pt]);
        long rank = gauss(evals_pt);
        target_degdet += rank;
    }

    if (actual_degdet != target_degdet)
        return false;

    return true;
}

bool is_interpolant_basis_geometric(
                                    const Mat<zz_pX> & intbas,
                                    const Vec<Mat<zz_p>> & evals,
                                    const zz_p & pt, // geometric case
                                    const long order,
                                    const VecLong & shift,
                                    const PolMatForm & form,
                                    const bool randomized
                                   )
{
    if (randomized)
        throw std::logic_error("==is_interpolant_basis== Fast randomized interpolant basis verification not implemented yet");

    // test that intbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,intbas,shift))
        return false;

    // test that the matrix consists of interpolants
    zz_pX_Multipoint_Geometric ev(pt, order);
    Vec<Mat<zz_p>> intbas_evals;
    ev.evaluate_matrix(intbas_evals, intbas);
    Mat<zz_p> res;
    for (long pt=0; pt<order; ++pt)
    {
        mul(res, intbas_evals[pt], evals[pt]);
        if (not IsZero(res))
            return false;
    }

    // test that intbas generates all interpolant, or equivalently,
    // that intbas has the right degree of determinant
    // --> since intbas is reduced, we can find the degree of determinant as
    // the sum of shifted row degrees minus the sum of shifts

    // compute degree of determinant of intbas
    VecLong rdeg;
    row_degree(rdeg, intbas, shift);
    long actual_degdet = std::accumulate(rdeg.begin(), rdeg.end(), 0) - std::accumulate(shift.begin(), shift.end(), 0);

    // compute sum of ranks of input matrices: this should be the
    // degree of determinant of intbas
    long target_degdet=0;
    for (long pt = 0; pt < order; ++pt)
    {
        Mat<zz_p> evals_pt(evals[pt]);
        long rank = gauss(evals_pt);
        target_degdet += rank;
    }

    if (actual_degdet != target_degdet)
        return false;

    return true;
}

bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Mat<zz_pX> & pmat,
                          const Vec<zz_p> & pts,
                          const VecLong & shift,
                          const PolMatForm & form,
                          const bool randomized
                         )
{
    zz_pX_Multipoint_General ev(pts);
    Vec<Mat<zz_p>> evals;
    ev.evaluate_matrix(evals, pmat);
    return is_interpolant_basis(intbas,evals,pts,shift,form,randomized);
}

bool is_interpolant_basis_geometric(
                                    const Mat<zz_pX> & intbas,
                                    const Mat<zz_pX> & pmat,
                                    const zz_p & pt, // geometric case
                                    const long order,
                                    const VecLong & shift,
                                    const PolMatForm & form,
                                    const bool randomized
                                   )
{
    zz_pX_Multipoint_Geometric ev(pt, order);
    Vec<Mat<zz_p>> evals;
    ev.evaluate_matrix(evals, pmat);
    return is_interpolant_basis_geometric(intbas,evals,pt,order,shift,form,randomized);
}


/*------------------------------------------------------------*/
/* M-Basis for interpolation points                           */
/*------------------------------------------------------------*/

// version with residual constant matrix computed at each iteration
void mbasis_rescomp(
                    Mat<zz_pX> & intbas,
                    const Vec<Mat<zz_p>> & evals,
                    const Vec<zz_p> & pts,
                    VecLong & shift,
                    long offset,
                    long order
                   )
{
#ifdef MBASIS_PROFILE
    double t_others=0.0,t_residual=0.0,t_intbas=0.0,t_kernel=0.0,t_now;
    t_now = GetWallTime();
#endif
    // A. General

    // A.1 dimensions of input matrix
    const long m = evals[offset].NumRows();
    const long n = evals[offset].NumCols();

    // A.3 store iota since it will be used at each iteration
    VecLong iota(m);
    std::iota(iota.begin(), iota.end(), 0);

    // B. Input representation; initialize output

    // B.1 vector of coefficients of output interpolant basis
    Vec<Mat<zz_p>> coeffs_intbas;

    // B.2 initially, intbas is the identity matrix
    coeffs_intbas.SetLength(1);
    ident(coeffs_intbas[0], m);
    // degree of interpolant basis, initially zero
    long deg_intbas = 0;
    // shifted row degree of intbas: initially equal to shift,
    // and shift will be continuously updated to always be the
    // "input-shift"-row degree of intbas

    // C. Residual matrix (m x n constant matrix, next coefficient
    // of intbas * pmat which we want to annihilate)

    // C.1 stores the residual, initially the first evaluation
    Mat<zz_p> residuals(evals[offset]);

    // C.2 temporary matrices used during the computation of residuals
    Mat<zz_p> intbas_eval;

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
    VecLong p_shift;

    // D.4 the constant kernel, and its permuted version
    Mat<zz_p> kerbas;
    Mat<zz_p> p_kerbas;

    // E. Updating intbas
    // stores the product "constant-kernel * coeffs_intbas[d]"
    Mat<zz_p> kerint; 

#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif

    for (long ord = offset; ord < offset+order; ++ord)
    {
#ifdef MBASIS_PROFILE
        t_now = GetWallTime();
#endif
        // compute permutation which realizes stable sort of shift
        // --> we need to permute things, to take into account the "priority"
        // (i.e. "weights") indicated by the shift; at this stage, the input
        // shift is the "input-shift"-row degree of intbas
        p_shift = iota;
        stable_sort(p_shift.begin(), p_shift.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (shift[a] < shift[b]);
                    } );

        // permute rows of the residual accordingly
        for (long i = 0; i < m; ++i)
            p_residual[i].swap(residuals[p_shift[i]]);
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
            // --> compute next residual

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

            // update residual, if not the last iteration
            if (ord<offset+order-1)
            {
                eval(intbas_eval, coeffs_intbas, pts[ord+1]);
                mul(residuals, intbas_eval, evals[ord+1]);
            }

            // update shift
            std::for_each(shift.begin(), shift.end(), [](long& a) { ++a; });
        }

        else if (ker_dim==m)
        {
            // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
            // --> interpolant basis is already correct for this point, no need to
            // change it or to change shift
            // --> we just need to compute the next residual
            // (unless we are at the last iteration, in which case the algorithm returns)
            if (ord<offset+order-1)
            {
#ifdef MBASIS_PROFILE
                t_now = GetWallTime();
#endif
                eval(intbas_eval, coeffs_intbas, pts[ord+1]);
                mul(residuals, intbas_eval, evals[ord+1]);
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
            // Also, deduce the degree of intbas
            bool deg_updated=false;
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                {
                    ++shift[i];
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

            // Update interpolant basis

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*coeffs_intbas (note: these rows currently have degree
            // at most deg_intbas)
            // TODO possible small improvement for uniform shift: these rows
            // have degree less than deg_intbas, in this case (and deg_intbas
            // is reached on the diagonal, among the pivot degrees)
            for (long d = 0; d <= deg_intbas; ++d)
            {
                mul(kerint, kerbas, coeffs_intbas[d]);
                for (long i = 0; i < ker_dim; ++i)
                    coeffs_intbas[d][pivind[perm_rows_ker[i]]].swap(kerint[i]);
            }

            // rows with !is_pivind are multiplied by X-pts[ord] (note: these rows
            // currently have degree less than deg_intbas)
            const zz_p pt(-pts[ord]);
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                    coeffs_intbas[deg_intbas][i] = coeffs_intbas[deg_intbas-1][i];
            for (long d = deg_intbas-1; d > 0; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                    {
                        mul(coeffs_intbas[d][i], coeffs_intbas[d][i], pt);
                        add(coeffs_intbas[d][i], coeffs_intbas[d][i], coeffs_intbas[d-1][i]);
                    }
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                    mul(coeffs_intbas[0][i], coeffs_intbas[0][i], pt);
#ifdef MBASIS_PROFILE
            t_intbas += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Find next residual: evaluation at pts[ord+1] of intbas*pmat
            // (this is not necessary if we are at the last iteration, i.e. ord==order-1)
            // we have finished: intbas*pmat is zero mod (X-pts[offset])...(X-pts[offset+order-1])
            if (ord<offset+order-1)
            {
                eval(intbas_eval, coeffs_intbas, pts[ord+1]);
                mul(residuals, intbas_eval, evals[ord+1]);
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
    // Convert interpolant basis to polynomial matrix representation
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
}

// version with vector residuals continuously updated
void mbasis_resupdate(
                      Mat<zz_pX> & intbas,
                      Vec<Mat<zz_p>> & evals,
                      const Vec<zz_p> & pts,
                      VecLong & shift,
                      long offset,
                      long order
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

    // A.2 store iota since it will be used at each iteration
    VecLong iota(m);
    std::iota(iota.begin(), iota.end(), 0);

    // B. Initialize output

    // B.1 vector of coefficients of output interpolant basis
    Vec<Mat<zz_p>> coeffs_intbas;

    // B.2 initially, intbas is the identity matrix
    coeffs_intbas.SetLength(1);
    ident(coeffs_intbas[0], m);
    // degree of interpolant basis, initially zero
    long deg_intbas = 0;
    // shifted row degree of intbas: initially equal to shift,
    // and shift will be continuously updated to always be the
    // "input-shift"-row degree of intbas

    // C. Residual matrix (holds the next constant matrix we want to annihilate; in this version, also holds further matrices that are continuously updated)

    // in this version "resupdate", we will continuously update
    // evals[offset],...,evals[offset+order-1] to get our residuals

    // C.1 temporary matrices used during the update of residuals: kernel * residual
    Mat<zz_p> kerres;

    // C.2 permuted residual, used as input to the kernel at the "base case"
    Mat<zz_p> p_residual;
    p_residual.SetDims(m, n);

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
    VecLong p_shift;

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
        // compute permutation which realizes stable sort of shift
        // --> we need to permute things, to take into account the "priority"
        // (i.e. "weights") indicated by the shift; at this stage, the input
        // shift is the "input-shift"-row degree of intbas
        p_shift = iota;
        stable_sort(p_shift.begin(), p_shift.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (shift[a] < shift[b]);
                    } );

        // permute rows of the residual accordingly
        for (long i = 0; i < m; ++i)
            p_residual[i] = evals[offset+ord][p_shift[i]];
        // Note: not using swap here, since we will need evals[offset] for
        // updating the residual
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
            // --> left-multiply intbas by (x-pt[offset+ord])
            // --> no need to update the residual:
            //    its j-th entry should be left-multiplied by (pt[j]-pt[offset+ord])
            //    (==> if pt[j] == pt[offset+ord], just set to zero --> now assuming
            //    points pairwise distinct)
            //    ==> but changing residual[j] by multiplying all its entries
            //    by the same nonzero constant does not change the interpolants

            // update intbas
            ++deg_intbas; // note that this degree is now > 0
            coeffs_intbas.SetLength(deg_intbas+1);
            coeffs_intbas[deg_intbas] = coeffs_intbas[deg_intbas-1];
            const zz_p pt(-pts[offset+ord]);
            for (long d=deg_intbas-1; d > 0; --d)
            {
                mul(coeffs_intbas[d], coeffs_intbas[d], pt);
                add(coeffs_intbas[d], coeffs_intbas[d], coeffs_intbas[d-1]);
            }
            mul(coeffs_intbas[0], coeffs_intbas[0], pt);

            // update shift accordingly
            std::for_each(shift.begin(), shift.end(), [](long& a) { ++a; });
        }

        else if (ker_dim < m)
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
            // Also, deduce the new degree of intbas
            bool deg_updated = false;
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                {
                    ++shift[i];
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

            // Update interpolant basis

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*coeffs_intbas (note: these rows have degree at most
            // deg_intbas)
            // TODO possible small improvement for uniform shift: these rows
            // have degree less than deg_intbas, in this case (and deg_intbas
            // is reached on the diagonal, among the pivot degrees)
            for (long d = 0; d <= deg_intbas; ++d)
            {
                mul(kerint, kerbas, coeffs_intbas[d]);
                for (long i = 0; i < ker_dim; ++i)
                    coeffs_intbas[d][pivind[perm_rows_ker[i]]].swap(kerint[i]);
            }

            // rows with !is_pivind are multiplied by X-pts[offset+ord] (note:
            // these rows currently have degree less than deg_intbas)
            const zz_p pt(-pts[offset+ord]);
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                    coeffs_intbas[deg_intbas][i] = coeffs_intbas[deg_intbas-1][i];
            for (long d = deg_intbas-1; d > 0; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                    {
                        mul(coeffs_intbas[d][i], coeffs_intbas[d][i], pt);
                        add(coeffs_intbas[d][i], coeffs_intbas[d][i], coeffs_intbas[d-1][i]);
                    }
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                    mul(coeffs_intbas[0][i], coeffs_intbas[0][i], pt);
#ifdef MBASIS_PROFILE
            t_intbas += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Update the evaluations offset+ord+1...offset+order-1
            // We do not consider the evals offset+ord, since it is known that
            // its update will be zero (and it will not be used in further
            // iterations)

            // Submatrix of rows of evals[k] corresponding to pivind are
            // replaced by kerbas*evals[k]
            for (long k = offset+ord+1; k < offset+order; ++k)
            {
                mul(kerres, kerbas, evals[k]);
                for (long i = 0; i < ker_dim; ++i)
                    evals[k][pivind[perm_rows_ker[i]]].swap(kerres[i]);
            }

            // rows with !is_pivind are multiplied by pts[k]-pts[offset+ord]
            for (long k = offset+ord+1; k < offset+order; ++k)
            {
                const zz_p mulpt(pts[k]+pt);
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                        mul(evals[k][i], evals[k][i], mulpt);
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
        // else if (ker_dim==m)
        // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
        // --> interpolant basis is already correct for this order, no need to
        // change it or to change shift
        // --> we just advance to the next residual matrix
    }

#ifdef MBASIS_PROFILE
    t_now = GetWallTime();
#endif
    // Convert interpolant basis to polynomial matrix representation
    conv(intbas, coeffs_intbas);
#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif
#ifdef MBASIS_PROFILE
    double t_total = t_residual + t_intbas + t_kernel + t_others;
    std::cout << "~~mbasis_resupdate~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_residual/t_total << "," << t_intbas/t_total << "," <<
    t_kernel/t_total << "," << t_others/t_total << std::endl;
#endif
}

void popov_mbasis(
                  Mat<zz_pX> & intbas,
                  Vec<Mat<zz_p>> & evals,
                  const Vec<zz_p> & pts,
                  VecLong & shift
                 )
{
    // compute first basis, after saving the evaluations and the shift
    Vec<Mat<zz_p>> copy_evals(evals);
    VecLong rdeg(shift);
    pmbasis(intbas,copy_evals,pts,rdeg,0,pts.length());
    copy_evals.kill();

    // shift for second call: negated pivot degree
    VecLong popov_shift(intbas.NumCols());
    std::transform(shift.begin(), shift.end(), rdeg.begin(),
                   popov_shift.begin(), std::minus<long>());

    // output shifted row degree
    shift=rdeg;

    // save `popov_shift` using `rdeg` as a buffer
    rdeg=popov_shift;

    // second call, basis shifted Popov up to constant transformation
    pmbasis(intbas,evals,pts,popov_shift,0,pts.length());

    // perform the constant transformation
    Mat<zz_p> lmat;
    row_leading_matrix(lmat, intbas, rdeg);
    inv(lmat, lmat);
    mul(intbas,lmat,intbas);
}



// version with residual constant matrix computed at each iteration,
// intbas is stored as evaluations at the points in pts
// TODO currently experimental, not properly tested
// TODO deal with case where intbas reaches degree = nbpoints
//VecLong mbasis_rescomp_eval(
//                            Vec<Mat<zz_p>> & intbas,
//                            const Vec<Mat<zz_p>> & evals,
//                            const Vec<zz_p> & pts,
//                            const VecLong & shift,
//                            long offset,
//                            long order
//                           )
//{
//#ifdef MBASIS_PROFILE
//    double t_others=0.0,t_residual=0.0,t_intbas=0.0,t_kernel=0.0,t_now;
//    t_now = GetWallTime();
//#endif
//    // A. General
//
//    // A.1 dimensions of input matrix
//    const long m = evals[0].NumRows();
//    const long n = evals[0].NumCols();
//
//    // A.2 store iota since it will be used at each iteration
//    VecLong iota(m);
//    std::iota(iota.begin(), iota.end(), 0);
//
//    // B. Input representation; initialize output
//
//    // B.1 initially, intbas is the identity matrix --> all evaluations are zero
//    intbas.SetLength(order);
//    for (long k = 0; k < order; ++k)
//        ident(intbas[k], m);
//    // shifted row degree of intbas, initially equal to shift
//    VecLong rdeg(shift);
//
//    // C. Residual matrix (m x n constant matrix, next coefficient
//    // of intbas * pmat which we want to annihilate)
//
//    // C.1 stores the residual, initially evals[0]
//    Mat<zz_p> residuals(evals[offset]);
//
//    // C.2 permuted residual, used as input to the kernel at the "base case"
//    Mat<zz_p> p_residual(INIT_SIZE, m, n);
//
//    // D. Base case (essentially amounts to finding the left kernel of the
//    // permuted residual p_residual)
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
//    // E. Updating intbas
//    // stores the product "constant-kernel * coeffs_intbas[d]"
//    Mat<zz_p> kerint; 
//
//#ifdef MBASIS_PROFILE
//    t_others += GetWallTime()-t_now;
//#endif
//
//    for (long ord = offset; ord < offset+order; ++ord)
//    {
//#ifdef MBASIS_PROFILE
//        t_now = GetWallTime();
//#endif
//        // compute permutation which realizes stable sort of rdeg
//        // --> we need to permute things, to take into account the "priority"
//        // (i.e. "weights") indicated by the shift; at this stage, the input
//        // shift is the shift-row degree 'rdeg' of intbas
//        p_rdeg = iota;
//        stable_sort(p_rdeg.begin(), p_rdeg.end(),
//                    [&](const long& a, const long& b)->bool
//                    {
//                    return (rdeg[a] < rdeg[b]);
//                    } );
//
//        // permute rows of the residual accordingly
//        for (long i = 0; i < m; ++i)
//            p_residual[i].swap(residuals[p_rdeg[i]]);
//#ifdef MBASIS_PROFILE
//        t_others += GetWallTime()-t_now;
//        t_now = GetWallTime();
//#endif
//        // find the (permuted) left kernel basis, hopefully in row echelon form;
//        kernel(p_kerbas,p_residual);
//#ifdef MBASIS_PROFILE
//        t_kernel += GetWallTime()-t_now;
//#endif
//        const long ker_dim = p_kerbas.NumRows();
//
//        if (ker_dim==0)
//        {
//            // Exceptional case: the residual matrix has empty left kernel
//            // --> left-multiply intbas by (x-pts[ord])
//            // --> compute next residual
//
//            // update intbas
//            for (long k=0; k<ord-offset; ++k)
//                mul(intbas[k], intbas[k], pts[offset+k]-pts[ord]);
//            for (long k=ord-offset+1; k<order; ++k)
//                mul(intbas[k], intbas[k], pts[offset+k]-pts[ord]);
//
//            // update residual, if not the last iteration
//            if (ord<offset+order-1)
//                mul(residuals, intbas[ord-offset+1], evals[ord+1]);
//
//            // update rdeg
//            std::for_each(rdeg.begin(), rdeg.end(), [](long& a) { ++a; });
//        }
//
//        else if (ker_dim==m)
//        {
//            // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
//            // --> interpolant basis is already correct for this point, no need to
//            // change it or to change rdeg
//            // --> we just need to compute the next residual
//            // (unless we are at the last iteration, in which case the algorithm returns)
//            if (ord<offset+order-1)
//            {
//#ifdef MBASIS_PROFILE
//                t_now = GetWallTime();
//#endif
//                mul(residuals, intbas[ord-offset+1], evals[ord+1]);
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
//#ifdef MBASIS_PROFILE
//            t_now = GetWallTime();
//#endif
//            // Compute pivots indices (pivot = rightmost nonzero entry)
//            // Experiments show that:
//            //   * kernel is expected to be of the form [ K | Id ]
//            //   * in general it is a column-permutation of such a matrix
//            // However note that a column-permutation is not sufficient for our needs
//            // Another property: if pivots are in the expected location (diagonal of
//            // rightmost square submatrix), then the corresponding column is the identity column.
//            bool expected_pivots = true;
//            for (long i = 0; i<ker_dim; ++i)
//            {
//                p_pivind[i] = m-1;
//                while (IsZero(p_kerbas[i][p_pivind[i]]))
//                    --p_pivind[i];
//                if (p_pivind[i] != m-ker_dim+i)
//                    expected_pivots = false;
//            }
//
//            if (not expected_pivots)
//            {
//                // find whether p_pivind has pairwise distinct entries
//                // (use pivind as temp space)
//                pivind = p_pivind;
//                std::sort(pivind.begin(), pivind.end());
//                // if pairwise distinct, then fine, the basis will not
//                // be Popov but will be ordered weak Popov (the goal of
//                // expected_pivots above was just to avoid this call to
//                // sort in the most usual case)
//
//                if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
//                {
//                    // the kernel is not in a shape we can deduce the intbas from (some pivots collide)
//                    // --> let's compute its lower triangular row echelon form
//                    // (use kerbas as temporary space)
//                    kerbas.SetDims(ker_dim,m);
//                    for (long i = 0; i < ker_dim; ++i)
//                        for (long j = 0; j < m; ++j)
//                            kerbas[i][j] = p_kerbas[i][m-1-j];
//                    image(kerbas, kerbas);
//                    // now column_permuted_ker is in upper triangular row echelon form
//                    for (long i = 0; i < ker_dim; ++i)
//                        for (long j = 0; j < m; ++j)
//                            p_kerbas[i][j] = kerbas[ker_dim-i-1][m-1-j];
//                    // and now p_kerbas is the sought lower triangular row echelon kernel
//
//                    // compute the actual pivot indices
//                    for (long i = 0; i<ker_dim; ++i)
//                    {
//                        p_pivind[i] = m-1;
//                        while (IsZero(p_kerbas[i][p_pivind[i]]))
//                            --p_pivind[i];
//                    }
//                }
//            }
//
//            // up to row permutation, the kernel is in "lower triangular" row
//            // echelon form (almost there: we want the non-permuted one)
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
//            // Now, update shifted row degree:
//            // entries corresponding to kernel pivot indices are kept, others are +1
//            for (long i = 0; i < m; ++i)
//                if (not is_pivind[i])
//                    ++rdeg[i];
//#ifdef MBASIS_PROFILE
//            t_others += GetWallTime()-t_now;
//            t_now = GetWallTime();
//#endif
//
//            // Update interpolant basis
//
//            // Submatrix of rows corresponding to pivind are replaced by
//            // kerbas*intbas
//            for (long d = 0; d < order; ++d)
//            {
//                mul(kerint, kerbas, intbas[d]);
//                for (long i = 0; i < ker_dim; ++i)
//                    intbas[d][pivind[perm_rows_ker[i]]].swap(kerint[i]);
//            }
//
//            // rows with !is_pivind are multiplied by pts[k]-pts[ord]
//            for (long k = 0; k < ord-offset; ++k)
//            {
//                const zz_p pt(pts[k+offset]-pts[ord]);
//                for (long i = 0; i < m; ++i)
//                    if (not is_pivind[i])
//                        mul(intbas[k][i], intbas[k][i], pt);
//            }
//            for (long i = 0; i < m; ++i)
//                if (not is_pivind[i])
//                    clear(intbas[ord-offset][i]);
//            for (long k = ord-offset+1; k < order; ++k)
//            {
//                const zz_p pt(pts[offset+k]-pts[ord]);
//                for (long i = 0; i < m; ++i)
//                    if (not is_pivind[i])
//                        mul(intbas[k][i], intbas[k][i], pt);
//            }
//#ifdef MBASIS_PROFILE
//            t_intbas += GetWallTime()-t_now;
//            t_now = GetWallTime();
//#endif
//            // Find next residual: evaluation at pts[ord+1] of intbas*pmat
//            // (this is not necessary if we are at the last iteration, i.e. ord==order-1)
//            // we have finished: intbas*pmat is zero mod X^order)
//            if (ord<offset+order-1)
//            {
//                mul(residuals, intbas[ord-offset+1], evals[ord+1]);
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
//#ifdef MBASIS_PROFILE
//    double t_total = t_residual + t_intbas + t_kernel + t_others;
//    std::cout << "~~mbasis_rescomp~~\t (residuals,basis,kernel,others): \t ";
//    std::cout << t_residual/t_total << "," << t_intbas/t_total << "," <<
//    t_kernel/t_total << "," << t_others/t_total << std::endl;
//#endif
//
//    // deduce pivot degree
//    for (long i = 0; i < m; ++i)
//        rdeg[i] -= shift[i];
//    return rdeg;
//}


void pmbasis_geometric(
                       Mat<zz_pX> & intbas,
                       const Mat<zz_pX> & pmat,
                       const zz_p & r,
                       const long order,
                       VecLong & shift,
                       Vec<zz_p> & pts
                      )
{
    zz_pX_Multipoint_Geometric eval(r, order);
    eval.get_points(pts);
    Vec<Mat<zz_p>> evals;
    eval.evaluate_matrix(evals, pmat);

    pmbasis_geometric(intbas, evals, pts, r, shift, 0, order);
}

void pmbasis_geometric(
                       Mat<zz_pX> & intbas,
                       Vec<Mat<zz_p>> & evals,
                       const Vec<zz_p> & pts,
                       const zz_p & r,
                       VecLong & shift,
                       long offset,
                       long order
                      )
{
    if (order <= 32)
    {
        mbasis(intbas, evals, pts, shift, offset, order);
        return;
    }

    // orders for recursive calls
    const long order1 = order/2;
    const long offset2 = offset+order1;
    const long order2 = order-order1;

    // first recursive call
    pmbasis_geometric(intbas, evals, pts, r, shift, offset, order1);

    // get the product of evaluations intbas(x_i) * pmat(x_i)
    // for the second half of the points

    // geometric progression of r, starting at pts[offset2]
    // note that we have by construction order2 >= order1 >= deg(intbas),
    // hence we can construct the multipoint object with parameter order2+1
    zz_pX_Multipoint_Geometric ev(r, pts[offset2], order2+1);
    Vec<Mat<zz_p>> intbas_eval;
    ev.evaluate_matrix(intbas_eval, intbas);
    for (long k = offset2; k < offset2+order2; ++k)
        mul(evals[k], intbas_eval[k-offset2], evals[k]);

    // second recursive call
    Mat<zz_pX> intbas2;
    pmbasis_geometric(intbas2, evals, pts, r, shift, offset2, order2);

    multiply(intbas,intbas2,intbas);
}


void pmbasis(
             Mat<zz_pX> & intbas,
             const Mat<zz_pX> & pmat,
             const Vec<zz_p> & pts,
             VecLong & shift
            )
{
    zz_pX_Multipoint_General eval(pts);
    Vec<Mat<zz_p>> evals;
    eval.evaluate_matrix(evals, pmat);

    pmbasis(intbas, evals, pts, shift, 0, pts.length());
}

void pmbasis(
             Mat<zz_pX> & intbas,
             Vec<Mat<zz_p>> & evals,
             const Vec<zz_p> & pts,
             VecLong & shift,
             long offset,
             long order
            )
{
    if (order <= 32)
    {
        mbasis(intbas, evals, pts, shift, offset, order);
        return;
    }

    // orders for the recursive calls
    const long order1 = order/2;
    const long offset2 = offset+order1;
    const long order2 = order-order1;

    // first recursive call
    pmbasis(intbas, evals, pts, shift, offset, order1);

    // Residual:  update the evaluations, for the second part of the points
    // --> product of evaluations intbas(x_i) * pmat(x_i)
    zz_pX_Multipoint_General ev(pts, offset2, order2);
    Vec<Mat<zz_p>> intbas_eval;
    ev.evaluate_matrix(intbas_eval, intbas);
    for (long k = offset2; k < offset2+order2; ++k)
        mul(evals[k], intbas_eval[k-offset2], evals[k]);

    // second recursive call
    Mat<zz_pX> intbas2;
    pmbasis(intbas2, evals, pts, shift, offset2, order2);

    multiply(intbas,intbas2,intbas);
}

void popov_pmbasis(
                   Mat<zz_pX> & intbas,
                   const Mat<zz_pX> & pmat,
                   const Vec<zz_p> & pts,
                   VecLong & shift
                  )
{
    // calling pmbasis twice on pmat would involve evaluating
    // it twice at the point --> avoid this, rather evaluate once for all
    zz_pX_Multipoint_General eval(pts);
    Vec<Mat<zz_p>> evals;
    eval.evaluate_matrix(evals, pmat);

    // compute first basis, after saving the evaluations and the shift
    Vec<Mat<zz_p>> copy_evals(evals);
    VecLong rdeg(shift);
    pmbasis(intbas,copy_evals,pts,rdeg,0,pts.length());
    copy_evals.kill();

    // shift for second call: negated pivot degree
    VecLong popov_shift(pmat.NumRows());
    std::transform(shift.begin(), shift.end(), rdeg.begin(),
                   popov_shift.begin(), std::minus<long>());

    // output shifted row degree
    shift=rdeg;

    // save `popov_shift` using `rdeg` as a buffer
    rdeg=popov_shift;

    // second call, basis shifted Popov up to constant transformation
    pmbasis(intbas,evals,pts,popov_shift,0,pts.length());

    // perform the constant transformation
    Mat<zz_p> lmat;
    row_leading_matrix(lmat, intbas, rdeg);
    inv(lmat, lmat);
    mul(intbas,lmat,intbas);
}

//VecLong pmbasis_geometric_multithread(
//                          Mat<zz_pX> & intbas,
//                          const Vec<Mat<zz_p>> & evals,
//                          const Vec<zz_p> & pts,
//                          const zz_p & r,
//                          const VecLong & shift
//                         )
//{
//    const long order = pts.length();
//    zz_pContext context;
//
//    if (order <= 32)
//        return mbasis(intbas, evals, pts, shift);
//
//    VecLong pivdeg; // pivot degree, first call
//    VecLong pivdeg2; // pivot degree, second call
//    VecLong rdeg(evals[0].NumRows()); // shifted row degree
//    long order1 = order/2; // order of first call
//    long order2 = order-order1; // order of second call
//    Mat<zz_pX> intbas2; // basis for second call
//
//    // first recursive call
//    Vec<zz_p>  pts1;
//    pts1.SetLength(order1);
//    for (long i = 0; i < order1; i++)
//        pts1[i] = pts[i];
//    pivdeg = pmbasis_geometric(intbas, evals, pts1, r, shift);
//
//    long max_pivdeg = pivdeg[0];
//    for (auto i : pivdeg)
//        if (max_pivdeg < i) max_pivdeg = i;
//
//    // shifted row degree = shift for second call = pivdeg+shift
//    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());
//
//    // get the product of evaluations intbas(x_i) * pmat(x_i)
//    // for the second half of the points
//    Vec<zz_p> pts2;
//    pts2.SetLength(order2);
//    for (long i=0; i<order2; ++i)
//        pts2[i] = pts[order1+i];
//
//    // geometric progression of r, starting at pts2[0]
//    // first ensure we have enough points for the degree
//    max_pivdeg = max(order2, max_pivdeg); 
//    zz_pX_Multipoint_Geometric ev(r,pts2[0], max_pivdeg+1);
//
//    Vec<Mat<zz_p>> evals2;
//    evals2.SetLength(order2);
//    for (long i = 0; i < order2; i++)
//        evals2[i].SetDims(intbas.NumRows(),intbas.NumCols());
//
//    context.save();  
//
//    // evaluate and store
//    NTL_EXEC_RANGE(intbas.NumRows(),first,last)
//
//    context.restore();    
//    for (long r = first; r < last; r++)
//        for (long c = 0; c < intbas.NumCols(); c++)
//        {
//            Vec<zz_p> val;
//            ev.evaluate(val, intbas[r][c]);
//
//            for (long i = 0; i < order2; i++)
//                evals2[i][r][c] = val[i];
//
//        }
//    NTL_EXEC_RANGE_END
//
//    context.save();  
//
//    // multiply and store    
//    NTL_EXEC_RANGE(order2,first,last)
//    context.restore();    
//
//    for (long i = first; i < last; i++)
//    {
//        evals2[i] = evals2[i] * evals[order1+i];
//    }
//    NTL_EXEC_RANGE_END
//
//    // second recursive call
//    pivdeg2 = pmbasis_geometric(intbas2, evals2, pts2, r, rdeg);
//
//    multiply(intbas,intbas2,intbas);
//
//    // final pivot degree = pivdeg1+pivdeg2
//    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());
//
//    return pivdeg;    
//}

//
//VecLong pmbasis_multithread(
//                Mat<zz_pX> & intbas,
//                const Vec<Mat<zz_p>> & evals,
//                const Vec<zz_p> & pts,
//                const VecLong & shift
//               )
//{
//    const long order = pts.length();
//    zz_pContext context;
//
//    if (order <= 32)
//        return mbasis(intbas, evals, pts, shift);
//
//    VecLong rdeg(evals[0].NumRows()); // shifted row degree
//    long order1 = order/2; // order of first call
//    long order2 = order-order1; // order of second call
//    Mat<zz_pX> intbas2; // basis for second call
//
//    // first recursive call
//    Vec<zz_p>  pts1;
//    pts1.SetLength(order1);
//    for (long i = 0; i < order1; i++)
//        pts1[i] = pts[i];
//    VecLong pivdeg = pmbasis_multithread(intbas, evals, pts1, shift);
//
//    // shifted row degree = shift for second call = pivdeg+shift
//    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());
//
//    // get the product of evaluations intbas(x_i) * pmat(x_i)
//    // for the second half of the points
//    Vec<zz_p> pts2;
//    pts2.SetLength(order2);
//    for (long i=0; i<order2; ++i)
//        pts2[i] = pts[order1+i];
//
//    zz_pX_Multipoint_General ev(pts2);
//
//    Vec<Mat<zz_p>> evals2;
//    evals2.SetLength(order2);
//    for (long i = 0; i < order2; i++)
//        evals2[i].SetDims(intbas.NumRows(),intbas.NumCols());
//
//    context.save();  
//
//    // evaluate and store
//    NTL_EXEC_RANGE(intbas.NumRows(),first,last)
//
//    context.restore();    
//    for (long r = first; r < last; r++)
//        for (long c = 0; c < intbas.NumCols(); c++)
//        {
//            Vec<zz_p> val;
//            ev.evaluate(val, intbas[r][c]);
//            for (long i = 0; i < order2; i++)
//                evals2[i][r][c] = val[i];
//
//        }
//    NTL_EXEC_RANGE_END
//
//    context.save();  
//
//    // multiply and store    
//    NTL_EXEC_RANGE(order2,first,last)
//    context.restore();    
//
//    for (long i = first; i < last; i++)
//    {
//        evals2[i] = evals2[i] * evals[order1+i];
//    }
//    NTL_EXEC_RANGE_END
//
//    // second recursive call
//    VecLong pivdeg2 = pmbasis_multithread(intbas2, evals2, pts2, rdeg);
//
//    multiply(intbas,intbas2,intbas);
//
//    // final pivot degree = pivdeg1+pivdeg2
//    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());
//
//    return pivdeg;
//}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
