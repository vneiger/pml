/**
 * \file nmod_poly_mat_approximant.h
 * \brief Minimal approximant bases
 * \version 0.0
 * \date 2022-06-25
 *
 * \todo CLEAN DOC
 */


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

#include "nmod_poly_mat_forms.h" // for testing form of approx basis, for orientation_t

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h> 

#ifdef __cplusplus
extern "C" {
#endif

/** @name General interfaces for approximant basis computation and verification */
//@{

/** General interface for approximant basis computation.
 *
 * Computes a `shift`-minimal approximant basis `appbas` for (`pmat`,`order`).
 * If the user requires a specific shifted form (minimal, ordered weak Popov,
 * Popov), the algorithm returns an approximant basis which has the required
 * form (possibly a stronger one). Approximants can be considered either as
 * left approximants (`row_wise` set to `true`) or as right approximants
 * (`row_wise` set to `false`).
 * TODO everywhere: add support for orientation
 *
 * If the user knows that the input matrix is sufficiently generic (for
 * example, it was randomly generated, working over a large prime field), then
 * it may be beneficial to set `generic` to `true`, which indicates that
 * specific code for this case may be used. For the moment, specific code has
 * been written only for matrices of dimensions n = 2m, with uniform shift and
 * uniform order; this does speed up computations by a factor ranging from
 * about 1x to about 2x (larger speed-up when the order is not very large).
 * Setting the generic flag in other contexts will have no effect (falling back
 * to the general algorithms).
 *
 * \todo With `generic` set to `true`, correctness is currently not guaranteed
 * (Monte Carlo algorithm). The user may want to use is_approximant_basis to
 * verify correctness.
 *
 * \param[out] appbas output approximant basis
 * \param[out] rdeg the `shift`-row degree of `appbas`
 * \param[in] pmat input polynomial matrix
 * \param[in] order input order
 * \param[in] shift input shift
 * \param[in] form indicates the required form for the output basis (see #poly_mat_form_t)
 * \param[in] row_wise indicates whether to compute left approximants (working row-wise) or right approximants (working column-wise)
 * \param[in] generic if `true`, specific code for generic instances may be used
 *
 * \todo The `row_wise` and `generic` flags are currently not supported.
 * \todo
 */
//void approximant_basis(
//                       appbas,
//                       rdeg,
//                       pmat,
//                       order,
//                       shift = vec of long's,
//                       form = ORD_WEAK_POPOV,
//                       row_wise = true,
//                       generic = false
//                      );

/** General user-friendly interface for approximant basis computation.
 *
 * Same function as above, but taking a single positive integer for the
 * parameter `order` instead of a list. This simply interprets this as
 * specifying a list of integers all equal to `order` (with the right length),
 * and then calls the function above.
 *
 * \todo `row_wise` false not handled
 * \todo
 */



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
 * \param[in] row_wise indicates whether to compute left approximants (working row-wise) or right approximants (working column-wise)
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo add parameter row_wise
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 * \todo WARNING! for the moment, does not really check generation!
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
 * \anchor mbasis1
 *
 *  These functions compute a shifted minimal or shifted Popov approximant
 *  basis using fast linear algebra on constant matrices.
 *
 *  Currently, this is only implemented for the uniform order `(1,...,1)`, following
 *  the algorithm _mbasis_ described in
 *  - P. Giorgi, C.-P. Jeannerod, G. Villard. Proceeding ISSAC 2003,
 *  - P. Giorgi, R. Lebreton. Proceedings ISSAC 2014,
 *  - C.-P. Jeannerod, V. Neiger, G. Villard. Preprint 2018.
 *
 *  The latter reference explicitly shows how to ensure that we obtain the
 *  canonical s-Popov approximant basis.
 *
 *  \todo implement the general Krylov-based approach from JNSV17, which
 *  efficiently solves the problem when the sum of the orders is small
 */
//@{

/** Computes a shifted Popov approximant basis at order `(1,...,1)` using fast
 * linear algebra on constant matrices. Thus, in this context, the input matrix
 * `pmat` is a constant matrix. This approximant basis consists of two sets of
 * rows:
 * - rows forming a left kernel basis in reduced row echelon form of some
 *   row-permutation of `pmat` (the permutation depends on the shift)
 * - the other rows are coordinate vectors multiplied by the variable `x`;
 *   there is one such row at each index `i` which is not among the pivot
 *   indices of the left kernel basis above.
 *
 * Here the whole approximant basis is stored in the OUT parameter `appbas`.
 *
 * \param[out] appbas polynomial matrix
 * \param[in] pmat constant matrix
 * \param[in] shift shift (vector of integers)
 *
 * \todo modify shift in place
 * \todo does it return Popov? if yes name it popov_mbasis1 and write mbasis1
 * if that makes sense; otherwise write other version popov_mbasis1
 * \todo safety correctness checks
 * \todo is output Popov by this method? check
 */
void mbasis1(nmod_poly_mat_t appbas,
             slong * res_shift,
             const nmod_mat_t mat,
             const slong * shift);

/** Computes a shifted Popov approximant basis at order `(1,...,1)` using fast
 * linear algebra on constant matrices. Thus, in this context, the input matrix
 * `pmat` is a constant matrix. This approximant basis consists of two sets of
 * rows:
 * - rows forming a left kernel basis in reduced row echelon form of some
 *   row-permutation of `pmat` (the permutation depends on the shift)
 * - the other rows are coordinate vectors multiplied by the variable `x`;
 *   there is one such row at each index `i` which is not among the pivot
 *   indices of the left kernel basis above.
 *
 * Here only the first set of rows is stored in the OUT parameter `kerbas`,
 * along with permutation...
 * \todo improve doc, to say exactly what this computes (and what is the return
 * value?)
 *
 * \param[out] kerbas constant matrix
 * \param[in] pmat constant matrix
 * \param[in] shift shift (vector of integers)
 *
 * \todo modify shift in place
 * \todo safety correctness checks
 * \todo is output Popov by this method? check
 */
slong mbasis1_for_mbasis(nmod_mat_t kerbas,
                         slong * res_shift,
                         slong * res_perm,
                         const nmod_mat_t mat,
                         const slong * shift);

//@} // doxygen group: Approximant basis via linear algebra

/** @name M-Basis algorithm (uniform approximant order)
 *
 * The core of these functions are implemented with `nmod_mat_poly_t` type,
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

/**********************************************************
*  TODO MBASIS TO BE FURTHER BENCH'D AND THEN INTEGRATED  *
**********************************************************/
/** Compute the minimal approximant basis of F for the order order and the shift shift
 *
 * It used the structured multiplication blocks to compute x^{-k} P_{k-1} F mod x and to compute
 * P_{k} = M P_{k-1}
 * Use polynomial matrix multiplication
 *
 */
void mbasis(nmod_poly_mat_t appbas,
            slong * res_shift,
            const nmod_poly_mat_t pmat,
            ulong order,
            const slong * shift);

/**********************************************************
*  TODO MBASIS TO BE FURTHER BENCH'D AND THEN INTEGRATED  *
**********************************************************/
/** Compute the minimal approximant basis of F for the order order and the shift shift
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] will be FIX and
 * compute iteratively P_{k} \in K[x]^{mxm}, then will transform to P_prime_{k} \in K^{mxm}[x] 
 * and compute x^{-k} P_{k-1} F mod x with a naive polynomial multiplication
 *
 * Use naive polynomial multiplication
 */
void mbasisII(nmod_poly_mat_t appbas,
               slong * res_shift,
               const nmod_poly_mat_t pmat,
               ulong order,
               const slong * shift);

/**********************************************************
*  TODO MBASIS TO BE FURTHER BENCH'D AND THEN INTEGRATED  *
**********************************************************/
/** Compute the minimal approximant basis of F for the order order and the shift shift
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_mat_poly.h)
 * Compute P_{k-1} \in K[x]^{mxm} with the structured_multiplication_blocks
 *
 * Use structured_multiplication_blocks and list_structured_multiplication_blocks 
 */
void mbasisIII(nmod_poly_mat_t appbas,
                slong * res_shift,
                const nmod_poly_mat_t pmat,
                ulong order,
                const slong * shift);

/**********************************************************
*  TODO MBASIS TO BE FURTHER BENCH'D AND THEN INTEGRATED  *
**********************************************************/
/** Compute the minimal approximant basis of F for the order order and the shift shift
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_mat_poly.h)
 * Compute P_{k-1} \in K^{mxm}[x] with the list_structured_multiplication_blocks
 * Transform P_{order - 1} as a poly_mat_t 
 *
 * Use list_structured_multiplication_blocks 
 */
void mbasisIV(nmod_poly_mat_t appbas,
               slong * res_shift,
               const nmod_poly_mat_t pmat,
               ulong order,
               const slong * shift);

/**********************************************************
*  TODO MBASIS TO BE FURTHER BENCH'D AND THEN INTEGRATED  *
**********************************************************/
/** Compute the minimal approximant basis of F for the order order and the shift shift
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] FIX and
 * Compute P_{k-1} \in K[x]^{mxm} with the list structured_multiplication_blocks
 * Compute x^{-k} P_{k-1} F mod x iteratively with a naive polynomial multiplication  
 *
 * Use naive polynomial multiplication and list_structured_multiplication_blocks 
 */
void mbasisV(nmod_poly_mat_t appbas,
              slong * res_shift,
              const nmod_poly_mat_t pmat,
              ulong order,
              const slong * shift);

// TODO DOC
// appbas must be initialized with right dimensions
FLINT_DLL void
nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
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
 * This is inspired from the algorithm _pmbasis_ described in
 *  - P. Giorgi, C.-P. Jeannerod, G. Villard. Proceeding ISSAC 2003,
 *  - P. Giorgi, R. Lebreton. Proceedings ISSAC 2014.
 */
//@{

/** 
 * TODO old doc, to be cleaned (new doc below, maybe not complete enough)
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] FIX and
 * Compute P_{k-1} \in K[x]^{mxm} with the list structured_multiplication_blocks
 * Compute x^{-k} P_{k-1} F mod x iteratively with a naive polynomial multiplication  
 *
 * Use naive polynomial multiplication and list_structured_multiplication_blocks 
 */
/** Computes a `shift`-ordered weak Popov approximant basis for `(pmat,order)`
 * using the algorithm PM-Basis (see @ref pmbasis) */
// TODO modify shift in place
void pmbasis(nmod_poly_mat_t appbas,
             slong * res_shift,
             const nmod_poly_mat_t pmat,
             ulong order,
             const slong * shift);

// TODO doc
// TODO middle_product currently naive
FLINT_DLL void
nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
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


#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_APPROXIMANT_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
