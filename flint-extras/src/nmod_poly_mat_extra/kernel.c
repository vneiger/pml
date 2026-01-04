/*
    Copyright (C) 2025 Gilles Villard
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h>
#include <math.h>
#include <stdlib.h>

#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_multiply.h"

/* general interface:
 * A/ compute rdeg and cdeg
 * ->suppress zero columns
 * ->trivialize zero rows
 * B/ handle case of constant matrices
 * C/ go for zls
 * */

slong nmod_poly_mat_kernel_via_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      const nmod_poly_mat_t pmat)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const slong d = nmod_poly_mat_degree(pmat);

    /* pmat == 0 => kernel is identity*/
    if (d == -1)
    {
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        return m;
    }

    /* FIXME may be reasonable to also have a case for d == 0 */

    /* pmat full row rank => empty kernel */
    /* full row rank implied by m == 1 or (m <= n && constant coefficient has full rank) */
    if (m == 1)
        return 0;

    if (m <= n)
    {
        nmod_mat_t coeff0;
        nmod_mat_init(coeff0, m, n, pmat->modulus);
        nmod_poly_mat_get_coeff_mat(coeff0, pmat, 0);
        /* could use: const slong r = nmod_mat_rank(coeff0);                */
        /* but let's call more directly its core computation                */
        /* (this uses pivind as temporary of length m for row permutation)  */
        const slong r = nmod_mat_lu(pivind, coeff0, 1);
        nmod_mat_clear(coeff0);
        if (r == m)
            return 0;
    }

    /* compute amplitude of `shift`, make `pivind` a temporary copy of `shift` */
    slong min = shift[0];
    slong amp = shift[0];
    for (slong i = 0; i < m; i++)
    {
        slong s = shift[i];
        pivind[i] = s;
        if (s < min)
            min = s;
        else if (s > amp)
            amp = s;
    }
    amp = amp - min;

    /* order for approximation: bound on rdeg_s(ker) + deg(pmat) + 1 */
    /* (see notes at end of file for rdeg_s(ker)) */
    const slong order = FLINT_MIN(m, n) * d + amp + d + 1;

    /* compute approximant basis and shifted row degree */
    nmod_poly_mat_pmbasis(ker, pivind, pmat, order);

    /* gather information for rows which belong to the kernel */
    /* note: at this stage, pivind == rdeg_s(appbas) */
    slong nz = 0;
    for (slong i = 0; i < m; i++)
    {
        if (pivind[i] - shift[i] + amp < order - d)
        {
            shift[nz] = pivind[i];
            pivind[nz] = i;
            for (slong j = 0; j < m; j++)
                FLINT_SWAP(nmod_poly_struct,
                           *nmod_poly_mat_entry(ker, nz, j),
                           *nmod_poly_mat_entry(ker, i, j));
            nz += 1;
        }
    }

    return nz;
}

/* Follows the description of Zhou-Labahn-Storjohann algorithm in [LNVZ22, Algo.1]  */
/* [LNVZ22] Labahn-Neiger-Vu-Zhou, Proceedings ISSAC 2022, arxiv.org/pdf/2202.09329 */
slong nmod_poly_mat_kernel_zls_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      nmod_poly_mat_t pmat)
{
    /* flint_printf("call %ld, %ld\n", pmat->r, pmat->c); */
    /* nmod_poly_mat_degree_matrix_print_pretty(pmat); */
    const slong m = pmat->r;
    const slong n = pmat->c;

    /* empty matrices, or single row/column vectors that are zero */
    if (m == 0 || n == 0
        || ((m == 1 || n == 1) && (nmod_poly_mat_max_length(pmat) == 0)))
    {
        /* flint_printf("base case empty||zero!\n"); */
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        return m;
    }

    /* one (nonzero) row vector */
    if (m == 1)
        return 0;

    /* one (nonzero) column vector */
    /* TODO */
    /* if (n == 1) */
    /*     return 0; */

    /* if n > m/2 (with here m >= 2), straightforward recursion by splitting */
    /* the matrix into two submatrices with n/2 columns                      */
    if (2*n > m)
    {
        /* we split pmat == [ n - n/2 cols | n/2 cols]                       */
        /* (we want more columns in left part, for the residual computation) */
        slong nz;

        nmod_poly_mat_t submat;    /* window only */
        nmod_poly_mat_t residual;  /* window only */
        nmod_poly_mat_t ker1;      /* window only */
        nmod_poly_mat_t ker2;      /* actual init */
        nmod_poly_mat_t ker3;      /* actual init */

        /* first recursive call on left submatrix */
        nmod_poly_mat_window_init(submat, pmat, 0, 0, m, n - n/2);
        nz = nmod_poly_mat_kernel_zls_approx(ker, pivind, shift, submat);
        nmod_poly_mat_window_clear(submat);

        /* if first kernel is empty, early exit */
        if (nz == 0)
            return 0;

        /* actual kernel: nz first rows */
        nmod_poly_mat_window_init(ker1, ker, 0, 0, nz, m);

        /* residual: (first ker) * (right submatrix) */
        /* overwrite left columns of pmat */
        nmod_poly_mat_window_init(submat, pmat, 0, n - n/2, m, n);
        nmod_poly_mat_window_init(residual, pmat, 0, 0, nz, n/2);
        nmod_poly_mat_mul_classical(residual, ker1, submat);
        nmod_poly_mat_window_clear(submat);

        /* recursive call 2, on residual */
        nmod_poly_mat_init(ker2, nz, nz, pmat->modulus);
        slong * pivind2 = FLINT_ARRAY_ALLOC(nz, slong);
        nz = nmod_poly_mat_kernel_zls_approx(ker2, pivind2, shift, residual);
        nmod_poly_mat_window_clear(residual);

        /* multiply bases and update pivind */
        nmod_poly_mat_init(ker3, nz, m, pmat->modulus);
        nmod_poly_mat_window_init(submat, ker2, 0, 0, nz, ker2->c);
        nmod_poly_mat_mul(ker3, submat, ker1);
        nmod_poly_mat_window_clear(submat);
        nmod_poly_mat_window_clear(ker1);

        for (slong i = 0; i < nz; i++)
            for (slong j = 0; j < m; j++)
                FLINT_SWAP(nmod_poly_struct,
                           *nmod_poly_mat_entry(ker, i, j),
                           *nmod_poly_mat_entry(ker3, i, j));

        for (slong i = 0; i < nz; i++)
            pivind[i] = pivind[pivind2[i]];

        nmod_poly_mat_clear(ker2);
        nmod_poly_mat_clear(ker3);
        flint_free(pivind2);

        return nz;
    }

    /* here, we are in the case 1 <= n <= m/2 */

    /* find rdeg, maxdeg, and diff_shift = min(shift - max(0, rdeg(pmat))) */
    slong * buf = FLINT_ARRAY_ALLOC(m, slong);
    nmod_poly_mat_row_degree(buf, pmat, NULL);
    slong maxdeg = buf[0];
    slong diff_shift = shift[0] - FLINT_MAX(0, buf[0]);
    for (slong i = 1; i < m; i++)
    {
        slong d = buf[i];
        slong diff = shift[i] - FLINT_MAX(0, d);
        if (d > maxdeg)
            maxdeg = d;
        if (diff_shift > diff)
            diff_shift = diff;
    }

    /* early exit: pmat == 0 => kernel is identity */
    if (maxdeg == -1)
    {
        /* flint_printf("second base case zero!\n"); */
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        flint_free(buf);
        return m;
    }

    /* build order for approximation (this choice diverges from [LNVZ22, Algo.1]): */
    /*     order = 1 + max(maxdeg, 1 + floor((sum(shift) - 1) / (m - n)))          */
    /* -> approx basis "generically" captures the whole kernel (see end of file)   */
    /* FIXME could be improved in case of unbalanced cdeg(pmat),                   */
    /* based on [LNVZ22, Theorem 3.3, Proof of Item (iii)]                         */
    slong order = - m * diff_shift - 1;
    for (slong i = 0; i < m; i++)
        order += shift[i];
    order = 1 + order / (m - n);
    order = 1 + FLINT_MAX(maxdeg, order);

    nmod_poly_mat_pmbasis(ker, shift, pmat, order);

    /* run an easy degree-based detection of rows in the kernel */
    /* permute ker into [kernel rows \\ other rows], preserving increasing pivots */
    slong nz = 0;
    for (slong i = 0; i < m; i++)
    {
        if (shift[i] < order + diff_shift)
        {
            pivind[nz] = i;
            nz += 1;
        }
        else
            buf[i - nz] = i;
    }
    for (slong i = nz; i < m; i++)
        pivind[i] = buf[i - nz];

    nmod_poly_mat_permute_rows(ker, pivind, shift);
    /* flint_printf("pmbasis ok, found %ld rows\n", nz); */
    /* flint_printf("pivind: %{slong*}\n", pivind, m); */
    /* nmod_poly_mat_degree_matrix_print_pretty(ker); */

    /* early exit */
    /* TODO */
//    // if the sum of pivot degrees in the part of kernel basis computed is
//    // equal to one of the upper bounds (here, sum of row degrees of pmat for
//    // non-pivot rows), then we know that we have captured the whole kernel
//    // --> this is often the case (e.g. in the generic situation mentioned
//    // above), and testing this avoids going through a few multiplications and
//    // the two recursive calls (which only lead to the conclusion that there is
//    // no new kernel row to be found)
//    long sum_matdeg = 0;
//    for (long i = 0; i < m2; ++i)
//        sum_matdeg += rdeg[rows_app_notker[i]];
//    long sum_kerdeg = 0;
//    for (long i = 0; i < m1; ++i)
//        sum_kerdeg += deg(appbas[rows_app_ker[i]][rows_app_ker[i]]);
//
//    const bool early_exit = (sum_kerdeg == sum_matdeg);
//
//    // if the whole kernel was captured according to the above test, or if the kernel is full (matrix was zero),
//    // then we have the whole kernel: just copy and return
//    if (early_exit)
//    {
//        kerbas.SetDims(m1, m);
//        for (long i = 0; i < m1; ++i)
//            kerbas[i].swap(appbas[rows_app_ker[i]]);
//        shift.swap(rdeg1);
//        return;
//    }

    /* proceed with remaining rows */
    nmod_poly_mat_t residual;  /* actual init */
    nmod_poly_mat_t ker2;      /* actual init */
    nmod_poly_mat_t prod;      /* actual init */
    nmod_poly_mat_t approx;    /* window only */
    nmod_poly_mat_t ker2nz;    /* window only */

    nmod_poly_mat_init(residual, m - nz, n, pmat->modulus);
    nmod_poly_mat_window_init(approx, ker, nz, 0, m, m);
    nmod_poly_mat_middle_product_naive(residual, approx, pmat, order, order + maxdeg + 1);
    /* FIXME provide (and use) general middle_product interface */

    /* TODO handle zero rows in G? */

    /* kernel of residual */
    nmod_poly_mat_init(ker2, m - nz, m - nz, pmat->modulus);
    const slong nz2 = nmod_poly_mat_kernel_zls_approx(ker2, buf, shift + nz, residual);
    /* flint_printf("second kernel ok, found %ld rows\n", nz2); */
    /* flint_printf("pivind: %{slong*}\n", buf, m - nz); */

    /* if nz2 == 0, whole kernel already computed, just return */
    if (nz2 == 0)
    {
        nmod_poly_mat_clear(residual);
        nmod_poly_mat_window_clear(approx);
        nmod_poly_mat_clear(ker2);
        flint_free(buf);
        return nz;
    }

    /* multiply to get missing part of kernel */
    /* flint_printf("nz2 == %ld\n", nz2); */
    /* nmod_poly_mat_degree_matrix_print_pretty(ker2); */
    /* nmod_poly_mat_degree_matrix_print_pretty(approx); */
    nmod_poly_mat_window_init(ker2nz, ker2, 0, 0, nz2, ker2->c);
    nmod_poly_mat_init(prod, nz2, m, pmat->modulus);
    nmod_poly_mat_mul(prod, ker2nz, approx);
    /* flint_printf("before swap::\n"); */
    /* nmod_poly_mat_degree_matrix_print_pretty(prod); */
    /* nmod_poly_mat_degree_matrix_print_pretty(ker); */
    for (slong i = 0; i < nz2; i++)
        for (slong j = 0; j < m; j++)
            FLINT_SWAP(nmod_poly_struct,
                       *nmod_poly_mat_entry(prod, i, j),
                       *nmod_poly_mat_entry(approx, i, j));
    /* flint_printf("after swap::\n"); */
    /* nmod_poly_mat_degree_matrix_print_pretty(ker); */

    /* update pivind */
    for (slong i = 0; i < nz2; i++)
        pivind[nz + i] = pivind[nz + buf[i]];

    /* permute rows to get increasing pivot indices */
    nmod_poly_mat_window_clear(approx);
    nmod_poly_mat_window_init(approx, ker, 0, 0, nz + nz2, m);
    slong i1 = 0;
    slong i2 = 0;
    for (slong i = 0; i < nz + nz2; i++)
    {
        if (i2 == nz2 || (i1 < nz && pivind[i1] < pivind[nz + i2]))
        {
            buf[i] = i1;
            i1++;
        }
        else
        {
            buf[i] = nz + i2;
            i2++;
        }
    }
    /* flint_printf("permutation buf: %{slong*}\n", buf, nz + nz2); */
    nmod_poly_mat_permute_rows(approx, buf, NULL);
    slong * tmp = FLINT_ARRAY_ALLOC(nz + nz2, slong);  /* use TMP_ALLOC */
    for (slong i = 0; i < nz + nz2; i++)
        tmp[i] = pivind[i];
    for (slong i = 0; i < nz + nz2; i++)
        pivind[i] = tmp[buf[i]];
    for (slong i = 0; i < nz + nz2; i++)
        tmp[i] = shift[i];
    for (slong i = 0; i < nz + nz2; i++)
        shift[i] = tmp[buf[i]];
    flint_free(tmp);

    nmod_poly_mat_window_clear(ker2nz);
    nmod_poly_mat_window_clear(approx);
    nmod_poly_mat_clear(residual);
    nmod_poly_mat_clear(ker2);
    nmod_poly_mat_clear(prod);
    flint_free(buf);

    /* nmod_poly_mat_degree_matrix_print_pretty(ker); */
    return nz + nz2;
}





/**
 *
 * The content of the file has not yet been tested very extensively
 *
 */

/**
 *
 * In place, M is m x n, shifted column degrees
 *   permute the columns in nondecreasing order
 *
 * Todo/to see: zero is set to degree 0 (not -1), for specific use
 *     in ZLS algorithm for the kernel
 *
 * Input:
 *  M  m x n
 *  ishift[m]
 *
 * Output:
 *  M is modified in place
 *  perm[n],initialized outside, the carried out permutation
 *  sdeg[n] initialized outside
 *
 */

// Dummy function for qsort
 int cmp(const void *a, const void *b)
{
    const slong *x = a;
    const slong *y = b;

    if (*x < *y) return -1;
    if (*x > *y) return 1;
    return 0;
}

void _nmod_poly_mat_sort_permute_columns_zls(nmod_poly_mat_t M, slong *sdeg, \
                                            slong *perm, const slong *ishift)
{

    slong n = M->c;
    slong j;

    nmod_poly_mat_column_degree(sdeg, M, ishift);

    for (j=0; j<n; j++)
        if (sdeg[j] < 0) sdeg[j]=0;

    _nmod_poly_mat_permute_columns_by_sorting_vec(M, n, sdeg, perm);

}

/**
 *
 *  Internal function
 *
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been
 *     sorted by shifted degree
 *
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 *
 *  TODO/TO SEE:
 *
 *  Input:
 *    A in m x n
 *    ishift[n], NULL (the degrees are computed) or initialized outside,
 *      the shift for the kernel
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the
 *       corresponding column degrees of A."
 *    kappa, a double >= 2, for the order of the order bases
 *              kappa * s instead of 3 *s in ZLS
 *
 *  Output:
 *    returns the dimension w of the kernel, which may be zero
 *    N, n x w, is initialized here only when w > 0,
 *           should be freed outside in that case only,
 *            gives a minimal basis of the right kernel
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shift degrees of the kernel basis
 *
 */

int nmod_poly_mat_zls_sorted(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t A, \
                                 const slong *ishift, const double kappa)
{

    slong i,j,k;

    slong m = A->r;
    slong n = A->c;


    slong min_mn;
    slong rho=0;

    // Tuning the order of the subsequent approximant
    // ----------------------------------------------

    long shift[n]; // temporary variable


    // In order to sort for computing the order of the approximant, for the test m=1 to be ok
    for (i=0; i<n; i++)
    {
        shift[i]=ishift[i];
    }


    qsort(shift, n, sizeof(slong), cmp);

    if (m <= n)
    {
        min_mn=m;
        for (i=n-m; i<n; i++)
            rho+=shift[i];
    }
    else
    {   // No sort needed, we consider all degrees
        min_mn=n;
        for (i=0; i<n; i++)
            rho+=ishift[i];
        // printf("\n m > n  %ld  %ld",m,n);   // To see, rectangular case
    }


    slong s;
    s = ceil((double) rho/min_mn);

    // Transposed A for the approximant PT on the left
    // PT will then be use without transposing

    nmod_poly_mat_t AT;
    nmod_poly_mat_init(AT, n, m, A->modulus);
    nmod_poly_mat_transpose(AT,A);

    nmod_poly_mat_t PT;
    nmod_poly_mat_init(PT, n, n, A->modulus);

    // Approximant PT
    // shift is modified in place
    // --------------------------

    for (i=0; i<n; i++)
    {
        shift[i]=ishift[i];
    }

    slong ks;
    ks=ceil((double) kappa*s);
    nmod_poly_mat_pmbasis(PT, shift, AT, ks+1);

    // Looking for zero residues and non zero residues
    //  the global residue is reused later
    // -----------------------------------------------

    nmod_poly_mat_t RT;
    nmod_poly_mat_init(RT, n, m, A->modulus);
    nmod_poly_mat_mul(RT,PT,AT);

    slong cdeg[n];
    nmod_poly_mat_row_degree(cdeg, RT, NULL);

    nmod_poly_mat_clear(AT);

    // Zero residue matrix P1  n x n1
    // nonzero residues matrix P2  n x n2
    //-----------------------------------

    slong n1=0;
    slong n2;

    for (j=0; j<n; j++) {
        if (cdeg[j]<0)
            n1+=1;
    }

    nmod_poly_mat_t P1;

    if (n1>0) {
        nmod_poly_mat_init(P1, n, n1, A->modulus);

        k=0;
        for (j = 0; j < n; j++)
        {
            if (cdeg[j]<0) {

                for (i = 0; i < n; i++)
                    nmod_poly_set(nmod_poly_mat_entry(P1, i, k), nmod_poly_mat_entry(PT, j, i));
                k+=1;
            }
        }
    }

    n2=n-n1;

    // the kernel is found
    if (n2==0) {
        nmod_poly_mat_init_set(N,P1);
        nmod_poly_mat_column_degree(degN, P1, ishift);

        nmod_poly_mat_clear(P1);
        return n1;
    }

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, n, n2, A->modulus);

    k=0;
    for (j = 0; j < n; j++)
    {
        if (cdeg[j]>=0) {

            for (i = 0; i < n; i++)
                nmod_poly_set(nmod_poly_mat_entry(P2, i, k), nmod_poly_mat_entry(PT, j, i));
            k+=1;
        }

    }

    //nmod_poly_mat_clear(PT);

    //  Special case m=1
    //  before the divide and conquer on m
    //  the kernel is found
    // -----------------------------------

    if (m==1){

        if (n1==0)
            return 0;
        else {

            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(degN, P1, ishift);

            nmod_poly_mat_clear(P1);
            return n1;
        }
    }

    //  Now n2 <> 0 and m > 1
    //    one can proceed to the divide and conquer from the rows
    //    of the nonzero residue P2, and of A.P2
    // --------------------------------------------------------

    slong * perm = flint_malloc(n2 * sizeof(slong));


    slong degP2[n2];

    // perm will be used also later for G from RT
    _nmod_poly_mat_sort_permute_columns_zls(P2,degP2,perm,ishift);

    for (i = 0; i < n2; i++) {
        degP2[i]=degP2[i]-ks;  // used below for the recursive calls
    }


    //+++++++++ OLD

    // nmod_poly_mat_t G;
    // nmod_poly_mat_init(G, m, n2, A->modulus);

    // nmod_poly_mat_mul(G, A, P2);

    // nmod_poly_mat_t TT;
    // nmod_poly_mat_init(TT, m, n2, A->modulus);


    // nmod_poly_mat_shift_right(TT,G,ks);


    // slong new_m=floor((double) m/2);

    // nmod_poly_mat_t G1;
    // nmod_poly_mat_init(G1, new_m, n2, A->modulus);

    // for (i = 0; i < new_m; i++){
    //     for (j = 0; j < n2; j++) {
    //         nmod_poly_set(nmod_poly_mat_entry(G1, i, j), nmod_poly_mat_entry(TT, i, j));
    //     }
    // }

    // nmod_poly_mat_t G2;
    // nmod_poly_mat_init(G2, m-new_m, n2, A->modulus);

    // for (i = 0; i < m-new_m; i++) {
    //     for (j = 0; j < n2; j++) {
    //         nmod_poly_set(nmod_poly_mat_entry(G2, i, j), nmod_poly_mat_entry(TT, i+new_m, j));
    //     }
    // }


    //++++++++  END OLD


    //++++++++++ NEW ++++++++++++++
    nmod_poly_mat_t TT;
    nmod_poly_mat_init(TT, m, n2, A->modulus);

    // We extract G (temporary TT) from RT as we have been extracting P2 from PT
    //   and apply above ordering of P2
    k=0;
    for (j = 0; j < n; j++)
    {
        if (cdeg[j]>=0) {

            for (i = 0; i < m; i++)
                nmod_poly_set(nmod_poly_mat_entry(TT, i, k), nmod_poly_mat_entry(RT, j, i));
            k+=1;
        }

    }
    nmod_poly_mat_clear(RT);

    nmod_poly_mat_permute_columns(TT, perm, NULL);

    nmod_poly_mat_t G;
    nmod_poly_mat_init(G, m, n2, A->modulus);

    nmod_poly_mat_shift_right(G,TT,ks);


    // We split G for the recursive call
    // ---------------------------------

    slong new_m=floor((double) m/2);

    nmod_poly_mat_t G1;
    nmod_poly_mat_init(G1, new_m, n2, A->modulus);

    for (i = 0; i < new_m; i++){
        for (j = 0; j < n2; j++) {
            nmod_poly_set(nmod_poly_mat_entry(G1, i, j), nmod_poly_mat_entry(G, i, j));
        }
    }

    nmod_poly_mat_t G2;
    nmod_poly_mat_init(G2, m-new_m, n2, A->modulus);

    for (i = 0; i < m-new_m; i++) {
        for (j = 0; j < n2; j++) {
            nmod_poly_set(nmod_poly_mat_entry(G2, i, j), nmod_poly_mat_entry(G, i+new_m, j));
        }
    }

    nmod_poly_mat_clear(TT);
    nmod_poly_mat_clear(G);
    //++++++++++++++  END NEW ++++++++++++++++++

    // Recursive calls
    // ---------------

    nmod_poly_mat_t N1;
    nmod_poly_mat_t N2;

    slong c1=0;
    slong c2=0;


    c1=nmod_poly_mat_zls_sorted(N1, degN, G1, degP2, kappa);


    if (c1 != 0) {

        nmod_poly_mat_t G3;
        nmod_poly_mat_init(G3, m-new_m, c1, A->modulus);

        nmod_poly_mat_mul(G3, G2, N1);

        for (i=0; i<c1; i++) {
            shift[i]=degN[i];
        }

        c2=nmod_poly_mat_zls_sorted(N2, degN, G3, shift, kappa);
        nmod_poly_mat_clear(G3);

    }


    nmod_poly_mat_clear(G1);
    nmod_poly_mat_clear(G2);

    // the recursive calls did not provide anything more

    if ((c1==0) || (c2==0)){

        if (n1==0) {
            return 0;
        }
        else {
            nmod_poly_mat_init_set(N,P1);
            nmod_poly_mat_column_degree(degN, P1, ishift);

            nmod_poly_mat_clear(P1);
            nmod_poly_mat_clear(P2);
            return n1;
        }
    }

    // A new part of the kernel has been found by the recursive calls, Q
    // -----------------------------------------------------------------


    nmod_poly_mat_t Q1;
    nmod_poly_mat_init(Q1, n, c1, A->modulus);

    nmod_poly_mat_mul(Q1, P2, N1);

    nmod_poly_mat_t Q;
    nmod_poly_mat_init(Q, n, c2, A->modulus);

    nmod_poly_mat_mul(Q, Q1, N2);

    nmod_poly_mat_clear(N1);
    nmod_poly_mat_clear(N2);
    nmod_poly_mat_clear(Q1);
    nmod_poly_mat_clear(P2);

    if (n1 ==0) {

        nmod_poly_mat_init_set(N,Q); // We should not need to copy
        nmod_poly_mat_column_degree(degN, Q, ishift);

        nmod_poly_mat_clear(Q);
        return c2;

    }
    else {
        nmod_poly_mat_init(N, n, n1+c2, A->modulus);

        for (i = 0; i < n; i++) {
            for (j = 0; j < n1; j++) {
                nmod_poly_set(nmod_poly_mat_entry(N, i, j), nmod_poly_mat_entry(P1,i,j));
            }
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < c2; j++) {
                nmod_poly_set(nmod_poly_mat_entry(N, i, j+n1), nmod_poly_mat_entry(Q, i, j));
            }
        }

        slong odeg[n1+c2];

        nmod_poly_mat_column_degree(odeg, P1, ishift);

        for (i=0; i<n1; i++) {
            degN[i]=odeg[i];
        }

        nmod_poly_mat_column_degree(odeg, Q, ishift);

        for (i=0; i<c2; i++) {
            degN[i+n1]=odeg[i];
        }

        nmod_poly_mat_clear(P1);
        nmod_poly_mat_clear(Q);

        return n1+c2;

    }

    return 0;
}


/**
 *
 *  Right shifted kernel of a polynomial matrix, assuming that the columns has been
 *     sorted by shifted degree
 *
 *  Algorithm of Wei Zhou, George Labahn, and Arne Storjohann
 *   "Computing Minimal Nullspace Bases"
 *    ISSAC 2012, https://dl.acm.org/doi/abs/10.1145/2442829.2442881
 *
 *  Calls nmod_poly_mat_zls_sorted after an initial sorting
 *
 *  TODO/TO SEE:
 *
 * Input:
 *    iA in m x n
 *     ishift[n], NULL (the degrees are computed) or initialized outside,
 *      the shift for the kernel
 *      values should be at least 0 (even for zero columns in A
 *      "with entries arranged in non-decreasing order and bounding the
 *       corresponding column degrees of A."
 *    kappa, a double >= 2, for the order of the order bases
 *              kappa * s instead of 3 *s in ZLS
 *
 *  Output:
 *    returns the dimension w of the kernel, which may be zero
 *    N, is initialized  n x n outside
 *       its first w columns give a minimal basis of the kernel
 *    degN[n], initialized outside, its first w entries are concerned,
 *        they are the ishift shifted degrees of the kernel basis
 *
 */

/**
 *  Calls nmod_poly_mat_zls_sorted after an initial sorting
 */

int nmod_poly_mat_kernel_zls(nmod_poly_mat_t N, slong *degN, const nmod_poly_mat_t iA, \
                         const slong *ishift, const double kappa)
{

    slong j,k;

    slong n = iA->c;


    nmod_poly_mat_t A;
    nmod_poly_mat_init_set(A, iA);

    slong * perm = flint_malloc(n * sizeof(slong));
    slong sdeg[n];

    // No input shift, simply the column degrees, then ordered
    //   a permutation is carried out in place
    // -------------------------------------------------------

    if (ishift == NULL) {
        nmod_poly_mat_column_degree(sdeg, A, NULL);

        for (j=0; j<n; j++)
            if (sdeg[j] < 0) sdeg[j]=0;

        _nmod_poly_mat_permute_columns_by_sorting_vec(A, n, sdeg, perm);
    }
    // Input shift, we sort
    // --------------------
    else {

        for (j=0; j<n; j++)
            sdeg[j]=ishift[j];

        _nmod_poly_mat_permute_columns_by_sorting_vec(A, n, sdeg, perm);
    }

    // Call to ZLS
    // -----------

    slong nz;
    slong tdeg[n];

    nmod_poly_mat_t NT;
    nz=nmod_poly_mat_zls_sorted(NT, tdeg, A, sdeg, kappa);


    // Undo the permutation for the kernel of the input matrix
    // -------------------------------------------------------

    if (nz !=0) {

        for (k = 0; k < n; k++) {
            for (j = 0; j < nz; j++){
                nmod_poly_set(nmod_poly_mat_entry(N, perm[k], j), nmod_poly_mat_entry(NT,k,j));
                degN[perm[k]]=tdeg[k];
            }
        }

        nmod_poly_mat_clear(NT); // Since nz > 0

        return nz;
    }

    flint_free(perm);
    nmod_poly_mat_clear(A);

    return 0;

}


/**
 * Experimental, should not be really considered
 *
 */

int nmod_poly_mat_approximant_kernel(nmod_poly_mat_t N, slong * FLINT_UNUSED(degN), const nmod_poly_mat_t A, \
                                 const slong *ishift)
{

    slong i,j,k;

    slong m = A->r;
    slong n = A->c;

    slong deg[n];

    nmod_poly_mat_column_degree(deg, A, NULL);

    slong degmax=deg[0];
    for (i=1; i<n; i++)
        if (degmax < deg[i])
            degmax=deg[i];

    slong min_nm=n;
    if (n >m)
        min_nm=m;

    nmod_poly_mat_t AT;
    nmod_poly_mat_init(AT, n, m, A->modulus);
    nmod_poly_mat_transpose(AT,A);

    nmod_poly_mat_t PT;
    nmod_poly_mat_init(PT, n, n, A->modulus);


    // Approximant PT
    // shift is modified in place
    // --------------------------

    slong shift[n];
    for (i=0; i<n; i++)
    {
        shift[i]=ishift[i];
    }

    nmod_poly_mat_pmbasis(PT, shift, AT, degmax*(min_nm +1)+1);


    // Looking for zero residues and non zero residues
    //  the global residue is reused later
    // -----------------------------------------------

    nmod_poly_mat_t RT;
    nmod_poly_mat_init(RT, n, m, A->modulus);
    nmod_poly_mat_mul(RT,PT,AT);

    slong cdeg[n];
    nmod_poly_mat_row_degree(cdeg, RT, NULL);

    nmod_poly_mat_clear(AT);
    nmod_poly_mat_clear(RT);


    // Zero residue matrix P1
    //-----------------------

    slong n1=0;

    for (j=0; j<n; j++) {
        if (cdeg[j]<0) {
            n1+=1;
        }
    }


    if (n1 > 0) {

        k=0;

        for (j=0; j<n; j++) {

            if (cdeg[j]<0) {
                for (i = 0; i < n; i++)
                    nmod_poly_set(nmod_poly_mat_entry(N, i, k), nmod_poly_mat_entry(PT, j, i));
                k+=1;
            }
        }

        nmod_poly_mat_clear(PT);
        return n1;
    }

    nmod_poly_mat_clear(PT);
    return 0;
}




/*------------------------------------------------------------*/
/*                          NOTES                             */
/*------------------------------------------------------------*/

/** Degree bounds for kernel bases.
 *
 * Assume pmat is m x n and has rank r. Let P be some shifted ordered weak
 * Popov kernel basis P of `pmat` (for a given, arbitrary shift).
 *
 * From Theorem 3.3 in
 *    Labahn, Neiger, Vu, Zhou
 *    Proceedings ISSAC 2022 (https://arxiv.org/pdf/2202.09329)
 * we get non-strict upper bounds on the sum of pivot degree of P:
 *   - it is <= sum of degrees of the r largest-degree rows of pmat
 *   - in particular, <= r * deg(pmat) <= min(m, n) * deg(pmat)
 *   - it is <= sum of degrees of the r largest-degree columns of pmat
 *   - in particular, <= sum of cdeg(pmat)
 * The third item above is not found directly in the mentioned Theorem 3.3,
 * but is deduced from its second item (part with deg(det(S)))
 *
 * The maximum of pivot degrees can be equal to their sum (although this is not
 * expected generically). This gives bounds on the non-pivot part of P: pick
 * one of the bounds above, and add the amplitude `max(shift) - min(shift)`.
 * Note that this might be pessimistic for shifts of large amplitude, since we
 * also have bounds coming from formula with adjugate of some submatrix of `pmat`.
 * For example, one can prove that the non-pivot part of the kernel basis in
 * shifted Popov form has degree at most n * deg(pmat)  (cf Lemma 12 of a research
 * internship report by Vu Thi Xuan).
 *
 * The maximum pivot degree also gives a bound on max(rdeg_s(ker)): pick one of
 * the bounds above for sum of pivot degrees, and add max(s) - min(s).
 */
