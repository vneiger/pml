#include <flint/nmod_mat.h>
#include <flint/perm.h>
#include <flint/profiler.h>
#include "nmod_mat_poly.h"
#include "nmod_mat_extra.h"

//#define MBASIS1_PROFILE
//#define MBASIS_PROFILE
//#define PMBASIS_PROFILE


/* type for stable sort while retaining the permutation */
typedef struct
{
    slong value;
    slong index;
} slong_pair;

/* comparator for quicksort, lexicographic (total) order to ensure stable sort */
static inline int _slong_pair_compare(const void * a, const void * b)
{
    slong_pair aa = * (const slong_pair *) a;
    slong_pair bb = * (const slong_pair *) b;
    if (aa.value == bb.value)
    {
        if (aa.index < bb.index)
            return -1;
        else if (aa.index > bb.index)
            return 1;
        else // aa.index == bb.index
            return 0;
    }
    else if (aa.value < bb.value)
        return -1;
    else // aa.value > bb.value
        return 1;
}

/** Creates a permutation from the sorting of a shift
 * After running this, perm is the unique list of integers which sorts
 * the pairs (shift,index) increasingly, i.e.
 * shift[perm[0]] <= shift[perm[1]] < ... < shift[perm[n-1]]
 *
 * \param perm permutation (list of integers), length n
 * \param shift list of integer to be sorted nondecreasingly, length n
 * \param pair_tmp temporary storage, length n
 * \param n length
 *
 */
static inline void _find_shift_permutation(slong * perm,
                                           const slong * shift,
                                           slong n,
                                           slong_pair * pair_tmp)
{
    for (slong i = 0; i < n; i++)
    {
        pair_tmp[i].value = shift[i];
        pair_tmp[i].index = i;
    }

    qsort(pair_tmp, n, sizeof(slong_pair), _slong_pair_compare);

    for (slong i = 0; i < n; i++)
        perm[i] = pair_tmp[i].index;
}





void nmod_mat_poly_mbasis(nmod_mat_poly_t appbas,
                          slong * shift,
                          const nmod_mat_poly_t matp,
                          ulong order)
{
#ifdef MBASIS_PROFILE
    timeit_t t_others,t_residual,t_appbas,t_kernel,t_now;
    timeit_start(t_now);
#endif
    // dimensions of input matrix
    const slong m = matp->r;
    const slong n = matp->c;
    const slong len = matp->length;

    // initialize output approximant basis with identity
    nmod_mat_poly_one(appbas);

    // residual matrix:
    // m x n constant matrix, next coefficient of appbas * matp to be
    // annihilated, initially coeffs_matp[0]
    nmod_mat_t res;
    nmod_mat_init_set(res, nmod_mat_poly_coeff(matp,0));

    // temporary matrix used during the computation of residuals
    nmod_mat_t res_tmp;
    nmod_mat_init(res_tmp, m, n, res->mod.n);

    // will hold the stable permutation making the shift nondecreasing
    slong * perm = _perm_init(m);
    slong_pair * pair_tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));

    // to store the left constant nullspace bases, their ranks,
    // and corresponding pivot information
    nmod_mat_t nsbas;
    long nullity;
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));

    // OTHERS?

    //// D.1 pivot indices in kernel basis (which is in row echelon form)
    //// Note: length is probably overestimated (usually kernel has m-n rows),
    //// but this avoids reallocating the right length at each iteration
    //VecLong pivind(m-1);
    //// Vector indicating if a given column index appears in this pivot index
    //// i.e. is_pivind[pivind[i]] = true and others are false
    //std::vector<bool> is_pivind(m, false);

    //// D.2 permutation for the rows of the constant kernel
    //VecLong perm_rows_ker;
    //// pivot indices of row echelon form before permutation
    //VecLong p_pivind(m-1);

    // and its permuted version
    //Mat<zz_p> p_nsbas;

    //// E. Updating appbas
    //// stores the product "constant-kernel * coeffs_appbas[d]"
    //Mat<zz_p> kerapp; 

#ifdef MBASIS_PROFILE
    timeit_stop(t_now);
    t_others->cpu += t_now->cpu;
    t_others->wall += t_now->wall;
#endif

    // Note: iterations will guarantee that `shift` (initially, this is the
    // input shift) holds the `shift`-shifted row degree of appbas; in
    // particular this is the case when the algorithm returns
    for (slong ord = 1; ord <= (slong)order; ++ord)
    {
#ifdef MBASIS_PROFILE
        timeit_start(t_now);
#endif
        // compute stable permutation which makes the shift nondecreasing
        // --> we will need to permute things, to take into account the
        // "priority" indicated by the shift
        _find_shift_permutation(perm, shift, m, pair_tmp);

        // permute rows of the residual accordingly
        nmod_mat_permute_rows(res, perm, NULL);

#ifdef MBASIS_PROFILE
        timeit_stop(t_now);
        t_others->cpu += t_now->cpu;
        t_others->wall += t_now->wall;
        timeit_start(t_now);
#endif

        // find the left nullspace basis of the (permuted) residual, in reduced
        // row echelon form and compact storage (see the documentation)
        nullity = nmod_mat_left_nullspace_compact(nsbas,pivots,res);
#ifdef MBASIS_PROFILE
        timeit_stop(t_now);
        t_kernel->cpu += t_now->cpu;
        t_kernel->wall += t_now->wall;
#endif

        if (nullity==0)
        {
            // Exceptional case: the residual matrix has empty left kernel
            // --> no need to compute more: the final basis is X^(order-ord+1)*appbas
            nmod_mat_poly_shift_left(appbas, appbas, order-ord+1);
            for (long i = 0; i < m; ++i)
                shift[i] += order-ord+1;
            return;
        }
    }

    // TODO CLEAR EVERYTHING
}


//        else if (ker_dim==m)
//        {
//            // Exceptional case: residual coeff was zero, and kernel 'nsbas' is identity
//            // --> approximant basis is already correct for this order, no need to
//            // change it or to change shift
//            // --> we just need to compute the next residual
//            // (unless ord == order, in which case the algorithm returns)
//            // this "residual" is the coefficient of degree ord in appbas * pmat:
//            // Note: at this point, res==0
//            if (ord<order)
//            {
//#ifdef MBASIS_PROFILE
//                t_now = GetWallTime();
//#endif
//                for (long d = std::max<long>(0,ord-coeffs_pmat.length()+1); d <= deg_appbas; ++d)
//                {
//                    mul(res_tmp, coeffs_appbas[d], coeffs_pmat[ord-d]);
//                    add(res, res, res_tmp);
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
//                while (IsZero(p_nsbas[i][p_pivind[i]]))
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
//                    // the kernel is not in a shape we can deduce the appbas from (some pivots collide)
//                    // --> let's compute its lower triangular row echelon form
//                    // (use nsbas as temporary space)
//                    nsbas.SetDims(ker_dim,m);
//                    for (long i = 0; i < ker_dim; ++i)
//                        for (long j = 0; j < m; ++j)
//                            nsbas[i][j] = p_nsbas[i][m-1-j];
//                    image(nsbas, nsbas);
//                    // now column_permuted_ker is in upper triangular row echelon form
//                    for (long i = 0; i < ker_dim; ++i)
//                        for (long j = 0; j < m; ++j)
//                            p_nsbas[i][j] = nsbas[ker_dim-i-1][m-1-j];
//                    // and now p_nsbas is the sought lower triangular row echelon kernel
//
//                    // compute the actual pivot indices
//                    for (long i = 0; i<ker_dim; ++i)
//                    {
//                        p_pivind[i] = m-1;
//                        while (IsZero(p_nsbas[i][p_pivind[i]]))
//                            --p_pivind[i];
//                    }
//                }
//            }
//
//
//            // up to row permutation, the kernel is in "lower triangular" row
//            // echelon form (almost there: we want the non-permuted one)
//            // prepare kernel permutation by permuting kernel pivot indices;
//            // also record which rows are pivot index in this kernel
//            // (note that before this loop, is_pivind is filled with 'false')
//            for (long i = 0; i < ker_dim; ++i)
//            {
//                pivind[i] = p_shift[p_pivind[i]];
//                is_pivind[pivind[i]] = true;
//            }
//
//            // perm_rows_ker = [0 1 2 ... ker_dim-1]
//            perm_rows_ker.resize(ker_dim);
//            std::copy_n(iota.begin(), ker_dim, perm_rows_ker.begin());
//            // permutation putting the pivot indices pivind in increasing order
//            sort(perm_rows_ker.begin(), perm_rows_ker.end(),
//                    [&](const long& a, const long& b)->bool
//                    {
//                    return (pivind[a] < pivind[b]);
//                    } );
//
//            // permute rows and columns of kernel back to original order
//            nsbas.SetDims(ker_dim,m);
//            for (long i = 0; i < ker_dim; ++i)
//                for (long j = 0; j < m; ++j)
//                    nsbas[i][p_shift[j]] = p_nsbas[perm_rows_ker[i]][j];
//#ifdef MBASIS_PROFILE
//            t_others += GetWallTime()-t_now;
//            t_now = GetWallTime();
//#endif
//
//            // Now, update shifted row degree:
//            // entries corresponding to kernel pivot indices are kept, others are +1
//            // Also, deduce the degree of appbas
//            bool deg_updated=false;
//            for (long i = 0; i < m; ++i)
//                if (not is_pivind[i])
//                {
//                    ++shift[i];
//                    if (not deg_updated && not IsZero(coeffs_appbas[deg_appbas][i]))
//                    { ++deg_appbas; deg_updated=true; }
//                }
//
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
//            // nsbas*coeffs_appbas (note: these rows currently have degree
//            // at most deg_appbas)
//            // TODO possible small improvement for uniform shift: these rows
//            // have degree less than deg_appbas, in this case (and deg_appbas
//            // is reached on the diagonal, among the pivot degrees)
//            for (long d = 0; d <= deg_appbas; ++d)
//            {
//                mul(kerapp, nsbas, coeffs_appbas[d]);
//                for (long i = 0; i < ker_dim; ++i)
//                    coeffs_appbas[d][pivind[perm_rows_ker[i]]].swap(kerapp[i]);
//            }
//
//            // rows with !is_pivind are multiplied by X (note: these rows
//            // currently have degree less than deg_appbas)
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
//            // Note: at this point, res==0
//            if (ord<order)
//            {
//                long dmin=std::max<long>(0,ord-coeffs_pmat.length()+1);
//                for (long d = dmin; d < deg_appbas+1; ++d) // we have deg_appbas <= ord
//                {
//                    mul(res_tmp, coeffs_appbas[d], coeffs_pmat[ord-d]);
//                    add(res, res, res_tmp);
//                }
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
//    std::cout << "~~mbasis_rescomp~~\t (residual,basis,kernel,others): \t ";
//    std::cout << t_residual/t_total << "," << t_appbas/t_total << "," <<
//    t_kernel/t_total << "," << t_others/t_total << std::endl;
//#endif
//}
//


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
