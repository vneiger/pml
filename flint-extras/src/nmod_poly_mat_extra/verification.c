/*
    Copyright (C) 2026 Vincent Neiger, Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_extra.h"  /* for prototypes */
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_interpolant.h"

/* TODO currently specialized to ROW_LOWER (or at least ROW_stuff) */
int nmod_poly_mat_is_approximant_basis(const nmod_poly_mat_t appbas,
                                       const nmod_poly_mat_t pmat,
                                       slong order,
                                       const slong * shift,
                                       orientation_t orient)
{
    const slong rdim = pmat->r;
    const slong cdim = pmat->c;
    const ulong prime = pmat->modulus;

    nmod_poly_mat_t residual;
    nmod_poly_t pol;
    nmod_mat_t CP0;

    nmod_poly_init(pol, prime);
    nmod_poly_mat_init(residual, rdim, cdim, prime);
    nmod_mat_init(CP0, rdim, rdim+cdim, prime);

    int success = 1;

    /* check appbas is square with the right dimension */
    if (appbas->r != rdim || appbas->c != rdim)
    {
        printf("basis has wrong row dimension or column dimension\n");
        success = 0;
    }

    /* check appbas is shifted reduced */
    if (!nmod_poly_mat_is_ordered_weak_popov(appbas, shift, orient))
    {
        printf("basis is not shifted-weak Popov\n");
        success = 0;
    }

    /* compute residual, check rows of appbas are approximants */
    nmod_poly_mat_multiply(residual, appbas, pmat);

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            nmod_poly_set_trunc(pol, nmod_poly_mat_entry(residual, i, j), order);
            if (!nmod_poly_is_zero(pol))
            {
                printf("entry %ld, %ld of residual has valuation less than the order\n",i,j);
                success = 0;
            }
        }
    }

    /* check generation: follows ideas from Algorithm 1 in Giorgi-Neiger, ISSAC 2018 */

    /* generation, test 1: check determinant of appbas is lambda * x**D      *
     * since ordered weak Popov, deg(det(appbas)) is sum of diagonal degrees *
     */
    slong D = 0;
    for (slong i = 0; i < rdim; i++)
        D += nmod_poly_degree(nmod_poly_mat_entry(appbas, i, i));
    nmod_poly_mat_det(pol, appbas);
    if (nmod_poly_degree(pol) != D)
    {
        printf("determinant is not lambda * x**(sum(diag-deg))");
        success = 0;
    }

    /* generation, test 2: check that [P(0)  C] has full rank *
     * where C = (appbas * pmat * X^{-order})  mod X          *
     * (coefficient "C" of degree order of the residual)      *
     **/
    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < rdim; j++)
        {
            ulong c = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(appbas, i, j), 0);
            nmod_mat_set_entry(CP0, i, j, c);
        }
        for (slong j = 0; j < cdim; j++)
        {
            ulong c = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(residual, i, j), order);
            nmod_mat_set_entry(CP0, i, rdim+j, c);
        }
    }

    slong rank = nmod_mat_rank(CP0);
    if (rank != rdim)
    {
        printf("generation test (step: [C P(0)] has full rank) failed");
        success = 0;
    }

    nmod_poly_clear(pol);
    nmod_poly_mat_clear(residual);
    nmod_mat_clear(CP0);

    return success;
}

/** First case: pairwise distinct points.
 *   
 * Mirrors is_approximant basis structure closely (same three checks: shape,
 * shift-reduced form, membership, generation), swapping the two pieces
 * that are genuinely order-truncation-specific for their points-based
 * counterparts. 
 *
 * - Membership: nmod_poly_mat_is_approximant_basis truncates the product
 *   appbas*pmat mod x^order and checks it is zero. Here, intbas is
 *   evaluated at each of the d points and multiplied by the corresponding
 *   evaluation E_k -- not a truncation, a genuine pointwise product), and
 *   each of the d products must be zero.
 *
 * - Generation (does intbas generate the whole interpolant module, not
 *   just a valid sub-basis?): nmod_poly_mat_is_approximant_basis runs two
 *   checks here:
 *     (1) deg(det(appbas)) == sum of diagonal degrees of appbas. 
           Ported here verbatim below, purely as the same kind of
 *         implementation-bug safety net (an independent nmod_poly_mat_det
 *         call cross-checked against the diagonal-degree sum). 
 *     (2) a Monte Carlo full-rank test on [P(0) | C] (C being the
 *         order-truncated residual coefficient at degree order) -- this is
 *         the test that actually certifies generation. 
 *   For interpolants, check (2)'s role is instead filled directly, with no
 *   Monte Carlo needed and no genericity caveat: since the points are
 *   pairwise distinct, K[x]/(prod_k(x-pts_k)) splits via CRT into
 *   direct_sum_k K[x]/(x-pts_k) = direct_sum_k K (one field per point).
 *   The interpolant module I contains M(x)*K[x]^{1xm} unconditionally
 *   (M(x) = prod_k(x-pts_k) vanishes at every point, so M(x)*e_i is always
 *   an interpolant), and under the CRT splitting, I / (M(x)*K[x]^{1xm})
 *   corresponds exactly to direct_sum_k ker_left(E_k), of K-dimension
 *   sum_k (m - rank(E_k)). Combined with the standard facts
 *   dim_K(K[x]^{1xm}/I) = deg(det(any basis of I)) and
 *   dim_K(K[x]^{1xm}/(M(x)*K[x]^{1xm})) = m*d, a short exact sequence gives,
 *   deg(det(any basis of I)) = sum_k rank(E_k).
 *
 *   So "deg(det(intbas)) == sum_k rank(E_k)" is an exact formula for
 *   the same quantity check (2).
 */


/* TODO currently specialized to ROW_LOWER (or at least ROW_stuff), same
 * limitation as nmod_poly_mat_is_approximant_basis in this same file. */
static int nmod_poly_mat_is_interpolant_basis_distinct(const nmod_poly_mat_t intbas,
                                       const ulong * pts,
                                       const nmod_mat_struct * E,
                                       slong d,
                                       const slong * shift,
                                       orientation_t orient)
{
    const slong rdim = intbas->r;
    const slong cdim = (d > 0) ? E[0].c : 0;
    const ulong prime = intbas->modulus;

    int success = 1;

    /* check intbas is square */
    if (intbas->r != intbas->c)
    {
        printf("basis has wrong row dimension or column dimension\n");
        success = 0;
    }

    /* check intbas is shifted reduced */
    if (!nmod_poly_mat_is_ordered_weak_popov(intbas, shift, orient))
    {
        printf("basis is not shifted-weak Popov\n");
        success = 0;
    }

    /* check rows of intbas are interpolants: intbas(pts[k]) * E_k == 0
     * for every one of the d points -- E_k is used directly (E + k) */
    nmod_mat_t Ik, prod;
    nmod_mat_init(Ik, rdim, rdim, prime);
    nmod_mat_init(prod, rdim, cdim, prime);

    for (slong k = 0; k < d; k++)
    {
        nmod_poly_mat_evaluate_nmod(Ik, intbas, pts[k]);
        nmod_mat_mul(prod, Ik, E + k);
        if (!nmod_mat_is_zero(prod))
        {
            printf("intbas(pts[%ld]) * E_%ld is not zero\n", (long) k, (long) k);
            success = 0;
        }
    }

    /** check generation, test 1 (self-consistency, mirrors PML's own
     * nmod_poly_mat_is_approximant_basis verbatim): deg(det(intbas)),
     * must equal D,
     * the sum of diagonal degrees of intbas (valid since an ordered weak
     * Popov basis has its row-i pivot exactly on the diagonal). Not a
     * distinct minimality criterion on its own -- a cheap cross-check 
     * between two different code paths, catching bugs in either 
     * nmod_poly_mat_is_ordered_weak_popov or this
     * diagonal-degree summation, kept purely as the same defense-in-depth
     * PML's own verifier already pays for. */
    slong D = 0;
    for (slong i = 0; i < rdim; i++)
        D += nmod_poly_degree(nmod_poly_mat_entry(intbas, i, i));

    nmod_poly_t det;
    nmod_poly_init(det, prime);
    nmod_poly_mat_det(det, intbas);
    if (nmod_poly_degree(det) != D)
    {
        printf("determinant degree (%ld) != sum of diagonal degrees (%ld)\n", (long) nmod_poly_degree(det), (long) D);
        success = 0;
    }
    nmod_poly_clear(det);

    /* check generation, test 2 (the actual generation/minimality test, see
     * this file's header comment for the full derivation): D must equal
     * sum_k rank(E_k), the points-based analogue of PML's own [P(0) | C]
     * Monte Carlo full-rank test -- here an exact, unconditional formula,
     * not an approximation. E_k used directly (E + k). */
    slong target_D = 0;
    for (slong k = 0; k < d; k++)
        target_D += nmod_mat_rank(E + k);

    if (D != target_D)
    {
        printf("sum of diagonal degrees (%ld) != sum of ranks of evaluation matrices (%ld)\n", (long) D, (long) target_D);
        success = 0;
    }

    nmod_mat_clear(Ik);
    nmod_mat_clear(prod);

    return success;
}


/* sum over distinct point values of the rank of the horizontally-stacked
   E's sharing that value -- see below for the
   full derivation. Collapses to sum_k rank(E_k) when every point is
   distinct (every group has multiplicity 1). */
static slong sum_rank_stacked_per_distinct_point(const nmod_mat_struct * E,
                                                 const ulong * pts, slong d,
                                                 slong m, slong n, ulong prime)
{
    slong target = 0;
    for (slong k = 0; k < d; k++)
    {
        int first = 1;
        for (slong j = 0; j < k; j++)
            if (pts[j] == pts[k]) { first = 0; break; }
        if (!first)
            continue;

        slong mult = 0;
        for (slong j = k; j < d; j++)
            if (pts[j] == pts[k])
                mult++;

        nmod_mat_t stacked;
        nmod_mat_init(stacked, m, n * mult, prime);
        slong col = 0;
        for (slong j = k; j < d; j++)
            if (pts[j] == pts[k])
            {
                for (slong i = 0; i < m; i++)
                    for (slong c = 0; c < n; c++)
                        nmod_mat_entry(stacked, i, col + c) = nmod_mat_entry(E + j, i, c);
                col += n;
            }
        target += nmod_mat_rank(stacked);
        nmod_mat_clear(stacked);
    }
    return target;
}

/** Formalizes the exploration of the non-distinct points case.
 *
 * In nmod_poly_mat_is_interpolant_basis_distinct, three checks,
 * only one depends on distinctness:
 *   1. shape (square) -- property of `intbas` alone, untouched.
 *   2. shift-ordered-weak-Popov form -- also a property of `intbas` alone
 *      (nmod_poly_mat_is_ordered_weak_popov, reused verbatim).
 *   3. membership: intbas(pts[k])*E_k == 0 for every k -- evaluate-and-
 *      check per k, independent of whether some k's share a point value.
 *      Untouched.
 *   4. generation: deg(det(intbas)) == sum_k rank(E_k). This is the one
 *      that changes. That formula is derived from K[x]/(prod_k(x-pts_k))
 *      splitting via CRT into direct_sum_k K[x]/(x-pts_k) -- valid only
 *      when the pts_k are pairwise distinct (coprime linear factors). A
 *      repeated point v makes (x-v) appear more than once in the product;
 *      repeated linear factors are not coprime, so that direct-sum
 *      decomposition does not hold, and summing rank(E_k) over every k
 *      independently overcounts whenever two E_k's sharing a point v have
 *      overlapping column space.
 *
 * Generalized formula TO CHECK.
 * Group the k's by their common point value v; for each distinct v,
 * horizontally concatenate (stack as columns) the E_k's sharing that v
 * into one m x (n*mult(v)) matrix, and take that matrix's rank -- not the
 * sum of the individual rank(E_k)'s, which can overcount when those E_k's
 * are not independent. Then
 *     deg(det(intbas)) == sum_{distinct v} rank(stacked E's at v).
 * This collapses to the ordinary sum_k rank(E_k) exactly when every point
 * is distinct (each group has multiplicity 1), so it is a strict
 * generalization of the real verifier's own generation check.
 */
static int nmod_poly_mat_is_interpolant_basis_nondistinct(const nmod_poly_mat_t intbas,
                                            const ulong * pts,
                                            const nmod_mat_struct * E,
                                            slong d,
                                            const slong * shift,
                                            orientation_t orient)
{
    slong m = intbas->r;
    slong n = (d > 0) ? E[0].c : 0;
    ulong prime = intbas->modulus;
    int success = 1;

    if (intbas->r != intbas->c)
    {
        printf("basis has wrong row dimension or column dimension\n");
        success = 0;
    }

    if (!nmod_poly_mat_is_ordered_weak_popov(intbas, shift, orient))
    {
        printf("basis is not shifted-weak Popov\n");
        success = 0;
    }

    nmod_mat_t Ik, prod;
    nmod_mat_init(Ik, m, m, prime);
    nmod_mat_init(prod, m, n, prime);
    for (slong k = 0; k < d; k++)
    {
        nmod_poly_mat_evaluate_nmod(Ik, intbas, pts[k]);
        nmod_mat_mul(prod, Ik, E + k);
        if (!nmod_mat_is_zero(prod))
        {
            printf("intbas(pts[%ld]) * E_%ld is not zero\n", (long) k, (long) k);
            success = 0;
        }
    }
    nmod_mat_clear(Ik);
    nmod_mat_clear(prod);

    slong D = 0;
    for (slong i = 0; i < m; i++)
        D += nmod_poly_degree(nmod_poly_mat_entry(intbas, i, i));

    nmod_poly_t det;
    nmod_poly_init(det, prime);
    nmod_poly_mat_det(det, intbas);
    if (nmod_poly_degree(det) != D)
    {
        printf("determinant degree (%ld) != sum of diagonal degrees (%ld)\n",
               (long) nmod_poly_degree(det), (long) D);
        success = 0;
    }
    nmod_poly_clear(det);

    slong target = sum_rank_stacked_per_distinct_point(E, pts, d, m, n, prime);
    if (D != target)
    {
        printf("sum of diagonal degrees (%ld) != sum of ranks of stacked evaluation matrices per distinct point (%ld)\n",
               (long) D, (long) target);
        success = 0;
    }

    return success;
}


/** Dispatcher between the distinct and non-distinct cases. 
 *  The two branches are mathematically equivalent on 
 * distinct points, so this is a performance dispatch, not a correctness fork.
 */
int nmod_poly_mat_is_interpolant_basis(const nmod_poly_mat_t intbas,
                                         const ulong * pts,
                                         const nmod_mat_struct * E,
                                         slong d,
                                         const slong * shift,
                                         orientation_t orient)
{
    int distinct = 1;
    for (slong k = 0; k < d && distinct; k++)
        for (slong j = 0; j < k; j++)
            if (pts[j] == pts[k])
            {
                distinct = 0;
                break;
            }

    if (distinct)
        return nmod_poly_mat_is_interpolant_basis_distinct(intbas, pts, E, d, shift, orient);
    else
        return nmod_poly_mat_is_interpolant_basis_nondistinct(intbas, pts, E, d, shift, orient);
}


/* TODO currently does not check generation */
int nmod_poly_mat_is_kernel(const nmod_poly_mat_t ker,
                            const slong * shift,
                            const nmod_poly_mat_t mat,
                            poly_mat_form_t form,
                            orientation_t orient)
{
    if (orient == ROW_UPPER || orient == ROW_LOWER)
    {
        const slong rdim = mat->r;
        const slong cdim = mat->c;
        const ulong prime = mat->modulus;

        nmod_poly_mat_t residual;
        nmod_poly_mat_init(residual, rdim, cdim, prime);

        int success = 1;

        /* check kernel has the right column dimension */
        if (ker->c != rdim)
        {
            printf("basis has wrong column dimension\n");
            success = 0;
        }

        /* check rank */
        slong rk = nmod_poly_mat_rank(mat);
        if (ker->r != rdim - rk)
        {
            printf("number of rows does not equal nullity\n");
            success = 0;
        }

        /* check kernel has the right form */
        if (!nmod_poly_mat_is_form(ker, form, shift, orient))
        {
            printf("basis does not have the required form\n");
            success = 0;
        }

        /* compute residual, check rows of ker are in the kernel */
        nmod_poly_mat_multiply(residual, ker, mat);
        if (!nmod_poly_mat_is_zero(residual))
        {
            printf("not all rows are in the kernel\n");
            success = 0;
        }

        nmod_poly_mat_clear(residual);

        return success;
    }

    if (orient == COL_UPPER || orient == COL_LOWER)
    {
        const slong rdim = mat->r;
        const slong cdim = mat->c;
        const ulong prime = mat->modulus;

        nmod_poly_mat_t residual;
        nmod_poly_mat_init(residual, rdim, cdim, prime);

        int success = 1;

        /* check kernel has the right column dimension */
        if (ker->r != cdim)
        {
            printf("basis has wrong row dimension\n");
            success = 0;
        }

        /* check rank */
        slong rk = nmod_poly_mat_rank(mat);
        if (ker->c != cdim - rk)
        {
            printf("number of rows does not equal nullity\n");
            success = 0;
        }

        /* check kernel has the right form */
        if (!nmod_poly_mat_is_form(ker, form, shift, orient))
        {
            printf("basis does not have the required form\n");
            success = 0;
        }

        /* compute residual, check rows of ker are in the kernel */
        nmod_poly_mat_multiply(residual, mat, ker);
        if (!nmod_poly_mat_is_zero(residual))
        {
            printf("not all rows are in the kernel\n");
            success = 0;
        }

        nmod_poly_mat_clear(residual);

        return success;
    }

    flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_is_kernel). Requested orientation not implemented.");
}
