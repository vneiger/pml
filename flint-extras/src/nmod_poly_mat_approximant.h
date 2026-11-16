/*
    Copyright (C) 2025 Vincent Neiger, Kevin Tran
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/


#include "pml.h"

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TODO                                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// ** shift: modify them in place or keep them as it is done now?

// ** mbasis: better compare different version, and clean/refactor

// ** pmbasis: handle different orders via max(order)-order shifting of pmat

// ** partial linearization for the case of few columns, at least in somehow
// generic case (way to verify result?)

// ** specific case of Pade approx [[F],[-1]]: there is some gain which can be
// obtained by tuning the higher recursion level, at least

// ** for uniform shift, implement generic case (not handling shifts/not doing
// permutations at all) and see if it is substantially faster than at least for
// some ranges of parameters

// ** shifted generic: from shift, deduce the pivot degree expected generically
// and use this as a shift instead of 'shift', then obtain directly Popov
// approx basis; skipping the first call to find the pivot degree (in addition,
// will this be more efficient because shift is nicer?)

// ** thresholds:
//   -- mbasis_resupdate / mbasis_rescomp depending on m/n
//   -- in pmbasis: base case vs recursion

// iterative version / mbasis / pmbasis: compare different ways of obtaining Popov
// (recompute with new shift, or maintain normal form, or work on output basis)


/** \file nmod_poly_mat_approximant.h
 * Definition (approximant basis).
 * -------------------------------
 * Consider:
 *   - an m x n matrix of univariate polynomials F,
 *   - an approximation order d (list of n positive integers).
 *
 * Then an approximant basis for (F,d) is a matrix over the univariate
 * polynomials whose rows form a basis for the following module:
 * { p in K[X]^{1 x m}  |  the column j of p F is 0 modulo X^{d[j]} }.
 * Note that such a matrix is square, m x m, and nonsingular. Its determinant
 * has the form X^k for some 0 <= k <= d[0] + d[1] + .. + d[n-1].
 */

/** \file nmod_poly_mat_approximant.h
 * Definition (shifted minimal approximant basis).
 * -----------------------------------------------
 * Starting from the definition of an approximant basis, consider further:
 *   - a degree shift s (a list of m integers).
 *
 * Then an approximant basis for (F,d) is said to be <em>a shift-minimal</em>
 * (resp. <em>a shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>)
 * approximant basis if it is in shift-reduced form (resp. in shift-ordered
 * weak Popov form, resp. in shift-Popov form). See nmod_poly_mat_forms.h
 * for definitions of these forms.
 */

/** \file nmod_poly_mat_approximant.h
 * Conventions.
 * ------------
 * Apart from the general interfaces (TODO) which offer the choice between left or
 * right approximants, all other functions compute left approximant bases
 * (approximants operate on the left of the matrix F; the basis elements are
 * the rows of the matrix).
 *
 * Most functions below use the following parameters.
 *
 * \param[out] appbas the output approximant basis (cannot alias `pmat`)
 * \param[in] pmat the input polynomial matrix (no restriction)
 * \param[in] order the input order (list of strictly positive integers, length must be the number of columns of `pmat`)
 * \param[in,out] shift in: the input shift; and out: the output shifted row degree of `appbas` (list of integers, length must be the number of rows of `pmat`)
 *
 * Note that the latter two restrictions on the lengths of the lists are
 * assuming left approximants; for right approximants, they are swapped.
 */


#ifndef NMOD_POLY_MAT_APPROXIMANT_H
#define NMOD_POLY_MAT_APPROXIMANT_H

#define PMBASIS_THRES 32

#include <flint/nmod_types.h>

#include "nmod_poly_mat_forms.h" // for orientation_t
#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name General interfaces for approximant basis computation and verification */
//@{

/** General interface for approximant basis computation.
 *
 * Computes a `shift`-ordered weak Popov approximant basis `P` of `F` at
 * order `orders[b]` for column `b`, or the `shift`-Popov one if `form` is
 * `POPOV`; with `form` among `NONE`, `REDUCED`, `WEAK_POPOV`,
 * `ORD_WEAK_POPOV`, the basis is `shift`-ordered weak Popov, which implies
 * these forms. Left approximants (rows), ROW_LOWER pivots.
 *
 * With `D` the sum of the orders and `n` the number of columns of positive
 * order, the algorithm is a divide and conquer on `D` (see
 * nmod_poly_mat_extra/approximant_basis_general.c):
 * - when all these columns have the same order, ::nmod_poly_mat_pmbasis
 *   (with its Krylov leaves, see @ref pmbasis);
 * - otherwise, when `n' < m` and `D <= c m / n'^(1/3)` for
 *   `n' = ceil(D / max(orders))` (the average number of columns per order,
 *   at most `n`), one call to ::nmod_poly_mat_approximant_basis_krylov (the
 *   `shift`-Popov basis);
 * - otherwise, when the orders are close to each other (`n max(orders) <=
 *   (1 + eps) D`), ::nmod_poly_mat_pmbasis on the input whose column `b` is
 *   multiplied by `x^(max(orders) - orders[b])`, at order `max(orders)`;
 * - otherwise, a first basis at the orders `min(orders[b], tau)`, with tau
 *   such that they sum to about `D/2`, a second one for the residual at the
 *   remaining orders (middle products, columns with windows of similar
 *   lengths grouped together), and their product.
 * The values of c and eps are NMOD_POLY_MAT_APPROXIMANT_BASIS_KRYLOV_CUTOFF
 * and NMOD_POLY_MAT_APPROXIMANT_BASIS_RAISE_SLACK (nmod_poly_mat_extra/impl.h).
 * For `POPOV`, unless the basis found is already in `shift`-Popov form, a
 * second basis is computed with the shift `-delta`, where `delta` is the
 * `shift`-pivot degree, and normalized by the inverse of its leading matrix,
 * which takes about the time of the first one.
 *
 * \param[out] P approximant basis, initialized with dimensions `m x m`
 * \param[out] rdeg the `shift`-row degree of `P` (length `m`), or `NULL`; it may alias `shift`
 * \param[in] F input polynomial matrix, `m x n'`
 * \param[in] orders list of `n'` nonnegative integers (0: no condition on that column)
 * \param[in] shift input shift (length `m`), or `NULL` for the zero shift
 * \param[in] form required form for the output basis (see #poly_mat_form_t), at most `POPOV`
 *
 * Size of the integers (see @ref IntegerSizes): it is assumed that
 * `max(shift) - min(shift) + D <= WORD_MAX`.
 */
void nmod_poly_mat_approximant_basis(nmod_poly_mat_t P,
                                     slong * rdeg,
                                     const nmod_poly_mat_t F,
                                     const slong * orders,
                                     const slong * shift,
                                     poly_mat_form_t form);

/** Verifying if a matrix is a minimal approximant basis.
 *
 * This checks whether the matrix `appbas` is a `shift`-minimal approximant
 * basis for (`pmat`,`order`) for the required form `form`. One may consider
 * left approximants (default, with `row_wise` set to `true`) or right
 * approximants (with `row_wise` set to `false`).
 *
 * \param[in] appbas approximant basis
 * \param[in] pmat polynomial matrix
 * \param[in] order order
 * \param[in] shift shift
 * \param[in] form required form for `appbas` (see #poly_mat_form_t)
 * \param[in] orient indicates the orientation (left/right approximants) and the definition of pivots
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo update documentation
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 * \todo WARNING! for the moment, hardcoded to check for ordered weak Popov
 */
int nmod_poly_mat_is_approximant_basis(const nmod_poly_mat_t appbas,
                                       const nmod_poly_mat_t pmat,
                                       slong order,
                                       const slong * shift,
                                       orientation_t orient);

//@} // doxygen group: General interfaces for approximant basis computation and verification


/** @name Iterative algorithms
 *
 * These functions implement iterative approximant basis algorithms inspired
 * from:
 * - B. Beckermann. J. Comput. Appl. Math. 40 (1992) 19-42
 * - B. Beckermann, G. Labahn. SIAM J. Matrix Anal. Appl. 15 (1994) 804-823
 * - M. Van Barel, A. Bultheel. Numer. Algorithms 3 (1992) 451-462
 */
//@{

/** Computes a `shift`-ordered weak Popov approximant basis for
 * `(pmat,order)`. At the end of the computation, the vector `shift` contains
 * the shifted row degree of `appbas`, for the input shift. The parameter
 * `order_wise` allows one to choose between two strategies for choosing the
 * next coefficient to deal with during the iteration, which may have some
 * impact on the timings depending on the input:
 * - process `pmat` order-wise (choose column with largest order)
 * - process `pmat` column-wise (choose leftmost column not yet completed). */
// TODO
//void appbas_iterative( ... );

/** Computes a `shift`-Popov approximant basis for `(pmat,order)`. At the end
 * of the computation, the vector `shift` contains the shifted row degree of
 * `appbas`, for the input shift. The parameter `order_wise` allows one to
 * choose between two strategies for choosing the next coefficient to deal with
 * during the iteration, which may have some impact on the timings depending on
 * the input:
 * - process `pmat` order-wise (choose column with largest order)
 * - process `pmat` column-wise (choose leftmost column not yet completed).
 *
 * \todo currently, this is based on the approach "compute a first basis which
 * gives the pivot degrees and re-compute with these degrees as the shift to
 * obtain the normal form". It would be interesting to implement the strategy
 * of Beckermann-Labahn 2000 which ensures Popov normal form, and compare it
 * with the current technique.
 */
// TODO
//void popov_appbas_iterative( ...);
//@} // doxygen group: Iterative algorithms


/** @name Approximant basis via linear algebra
 * \anchor approx_linalg
 *
 *  These functions compute a shifted minimal or shifted Popov approximant
 *  basis using fast linear algebra on constant matrices.
 *
 *  \todo implement the general Krylov-based approach from JNSV17, which
 *  efficiently solves the problem when the sum of the orders is small
 */
//@{

//@} // doxygen group: Approximant basis via linear algebra

/** @name M-Basis algorithm (uniform approximant order)
 *
 * The core of these functions is implemented with `nmod_mat_poly_t` type,
 * for efficiency reasons. The bulk of the documentation can be found
 * in .nmod_mat_poly.h .
 *
 */
//@{

/** Computes the `shift`-Popov approximant basis for `(pmat,order)`,
 * relying on the `mbasis` approach.
 *
 * \todo currently, this is based on the approach "compute a first basis which
 * gives the pivot degrees and re-compute with these degrees as the shift to
 * obtain the normal form". It would be interesting to implement a strategy
 * which maintains a shifted Popov form along the iterations, and compare it
 * with the current technique. Another reasonable strategy for practical
 * improvement is to make sure almost no computation is necessary in "generic
 * cases" (only for uniform shift?).
 * \todo integrate
 */
//void popov_mbasis( ... );

// TODO DOC
// appbas must be initialized with right dimensions
/** Size of the integers (see @ref IntegerSizes): as for
 * ::nmod_mat_poly_mbasis, it is assumed that `order <= WORD_MAX` and
 * `max(shift) + order <= WORD_MAX` for the input `shift`. */
void nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                          slong * shift,
                          const nmod_poly_mat_t pmat,
                          ulong order);

//@} // doxygen group: M-Basis algorithm (uniform approximant order)


/** @name PM-Basis algorithm (uniform approximant order)
 * \anchor pmbasis
 *
 * These functions compute a `shift`-minimal ordered weak Popov approximant
 * basis for `(pmat,orders)` in the case where `orders` is given by a single
 * integer `orders = (order,...order)`. They use a divide and conquer approach,
 * computing a first basis at order `order/2`, finding the so-called _residual
 * matrix_, computing a second basis at order `order/2`, and deducing the
 * sought basis by multiplying the two obtained bases.
 *
 * The first recursive call returns an approximant basis `appbas1` such that
 * `appbas1*pmat = 0 mod x^{order/2}`, and the residual matrix has the same
 * dimensions as `pmat` and is defined by the matrix middle product
 * `(x^{-order/2} appbas1*pmat) mod x^{order/2}`.
 *
 * At the end of the computation, the vector `shift` contains the shifted row
 * degree of `appbas`, for the input shift. 
 *
 * Krylov leaves: at each node of the recursion, with `m x n` the dimensions
 * of `pmat` and `order` the order at this node, if `n < m` and `n order <=
 * c m / n^(1/3)` (for a constant c, see below), the basis of this node
 * is computed directly by ::nmod_poly_mat_approximant_basis_krylov (linear
 * algebra on a Krylov matrix, about `D^2 (D + m)` operations for `D = n
 * order`), which is faster than the recursion for inputs with few columns and
 * an order that is not too large compared to `m / n`. This gives the
 * `shift`-Popov basis of the node; when such a leaf is used at the top of the
 * recursion, `appbas` is the `shift`-Popov basis, and in any case it is
 * `shift`-ordered weak Popov. The value of `c` is
 * NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF (nmod_poly_mat_extra/impl.h).
 *
 * Size of the integers (see @ref IntegerSizes): as in `mbasis`, each entry of
 * `shift` increases by at most `order` during the computation; it is assumed
 * that `max(shift) + order <= WORD_MAX` for the input `shift`.
 *
 * This is inspired from the algorithm _pmbasis_ described in
 *  - P. Giorgi, C.-P. Jeannerod, G. Villard. Proceeding ISSAC 2003,
 *  - P. Giorgi, R. Lebreton. Proceedings ISSAC 2014.
 */
//@{

/** Computes a `shift`-ordered weak Popov approximant basis for `(pmat,order)`
 * using the algorithm PM-Basis (see @ref pmbasis) */

// TODO doc
// TODO middle_product currently naive


void nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           slong order);


/** Computes a `shift`-Popov approximant basis for `(pmat,order)` using the
 * algorithm PM-Basis (see @ref pmbasis) twice: the first call yields an
 * ordered weak Popov form which indicates the `shift`-pivot degree, which is
 * then used as a shift in the second call to obtain the sought canonical
 * basis.
 *
 * \todo often (generic case), the second call is not necessary and one just
 * has to multiply by the inverse of the leading matrix. Implement this(?). See
 * also related comment for popov_mbasis.
 **/
// TODO integrate
//void popov_pmbasis( ... );

//@} // doxygen group: PM-Basis algorithm (uniform approximant order)


/** @name Approximant and M-Pade bases via Krylov matrices
 *
 * These functions compute the basis in `shift`-Popov form (row-wise, with
 * ROW_LOWER pivots, see @ref MatrixForms), which is unique, of the module of
 * vectors `p` of `m` polynomials such that, for each column `b` of the
 * `m x n` matrix `F`, with `d_b = orders[b]`:
 * - for ::nmod_poly_mat_approximant_basis_krylov, the product `p F[:,b]` is
 *   zero modulo `x^{d_b}`;
 * - for ::nmod_poly_mat_mpade_basis_krylov, the product `p F[:,b]` is zero
 *   modulo `(x - pts[b])^{d_b}` (points may be repeated).
 *
 * With `D = d_0 + ... + d_{n-1}`, they form the `m x D` matrix `V` whose
 * block `b` of `d_b` columns holds the first `d_b` coefficients of the
 * expansion of `F[:,b]` at `pts[b]` (at 0 for approximants), and compute the
 * kernel of its Krylov matrix with respect to the `D x D` Jordan matrix with
 * `n` blocks, the block `b` having eigenvalue `pts[b]` and size `d_b` (see
 * @ref KrylovMatrices, where `D` and the block sizes `D_b = d_b` have the
 * same meaning, and ::nmod_poly_mat_krylov_kernel_basis_jordan, which
 * computes general relation bases, in nmod_poly_mat_relation_basis.h).
 *
 * M-Pade approximation covers the two other problems: approximants are the
 * case of all points zero, and interpolants
 * (::nmod_poly_mat_interpolant_basis_krylov) the case of a constant matrix
 * `F = V` with all orders 1. The order `d_b` at `pts[b]` involves the
 * derivatives of `p F[:,b]` at `pts[b]`, whereas in the interpolant function
 * a point repeated in several columns gives several conditions of order 1
 * at that point, with different vectors of values.
 *
 * The output `P` must be initialized with dimensions `m x m`; `shift == NULL`
 * means the zero shift; the orders must be nonnegative, and the modulus
 * prime. Size of the integers (see @ref IntegerSizes): there is no condition
 * on the shift, whose entries may be arbitrary `slong` values; `D` is the
 * number of columns of `V`, so it is bounded by the available memory.
 *
 * The cost is about `D^2 (D + m)` operations, whatever the shift: this is
 * suited to `D` not much larger than `m`, with few columns. For example, for
 * one column with `D` up to `8 m`, this is 2 to 28 times faster than
 * ::nmod_poly_mat_mbasis (which gives an ordered weak Popov basis), whereas
 * for many columns of small order (`n = m/2` columns at order 8 or more)
 * ::nmod_poly_mat_mbasis and ::nmod_poly_mat_pmbasis are faster (Intel Xeon
 * with AVX-512, 30- and 60-bit primes).
 */
//@{

/** Approximant basis in `shift`-Popov form of `F` at order `orders`, by
 * linear algebra on Krylov matrices. */
void nmod_poly_mat_approximant_basis_krylov(nmod_poly_mat_t P,
                                            const nmod_poly_mat_t F,
                                            const slong * orders,
                                            const slong * shift);

/** M-Pade basis in `shift`-Popov form of `F` at the points `pts` and orders
 * `orders`, by linear algebra on Krylov matrices. */
void nmod_poly_mat_mpade_basis_krylov(nmod_poly_mat_t P,
                                      const nmod_poly_mat_t F,
                                      nn_srcptr pts,
                                      const slong * orders,
                                      const slong * shift);

//@} // doxygen group: Approximant and M-Pade bases via Krylov matrices

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_APPROXIMANT_H
