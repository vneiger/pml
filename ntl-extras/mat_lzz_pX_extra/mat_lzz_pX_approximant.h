#ifndef MAT_LZZ_PX_APPROXIMANT__H
#define MAT_LZZ_PX_APPROXIMANT__H

/** \brief Minimal approximant bases.
 *
 * \file mat_lzz_pX_approximant.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-07
 *
 */

/** \file mat_lzz_pX_approximant.h
 * Definition (approximant basis).
 * -------------------------------
 * Consider
 *   - an m x n matrix of univariate polynomials F,
 *   - an approximation order d (list of n positive integers).
 *
 * Then an approximant basis for (F,d) is a matrix over the univariate
 * polynomials whose rows form a basis for the following module:
 * { p in K[X]^{1 x m}  |  the column j of p F is 0 modulo X^{d[j]} }.
 * Note that such a matrix is square, m x m, and nonsingular.
 */

/** \file mat_lzz_pX_approximant.h
 * Definition (shifted minimal approximant basis).
 * -----------------------------------------------
 * Considering furthermore:
 *   - a degree shift s (a list of m integers).
 *
 * Then an approximant basis for (F,d) is said to be <em>a shift-minimal</em>
 * (resp. <em>a shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>)
 * approximant basis if it is in shift-reduced form (resp. in shift-ordered
 * weak Popov form, resp. in shift-Popov form). See mat_lzz_pX_forms.h
 * for definitions of these forms.
 */

/** \file mat_lzz_pX_approximant.h
 * Conventions.
 * ------------
 * Apart from the general interfaces which offer the choice between left or
 * right approximants, all other functions compute left approximant bases
 * (approximants operate on the left of the matrix F; the basis elements are
 * the rows of the matrix).
 *
 * Most functions below use the following parameters.
 *
 * \param[out] appbas the output approximant basis (cannot alias `pmat`)
 * \param[in] pmat the input polynomial matrix (no restriction)
 * \param[in] order the input order (list of strictly positive integers, length must be `pmat.NumCols()`)
 * \param[in,out] shift the input shift and the output shifted row degree of `appbas` (list of integers, length must be `pmat.NumRows()`)
 *
 * Note that the latter two restrictions on the lengths of the lists are
 * assuming left approximants; for right approximants, they are swapped.
 */

#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_forms.h"
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

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
 * \param[out] rdeg the `shifts`-row degree of `appbas`
 * \param[in] pmat input polynomial matrix
 * \param[in] order input order
 * \param[in] shift input shift
 * \param[in] form indicates the required form for the output basis (see #PolMatForm)
 * \param[in] row_wise indicates whether to compute left approximants (working row-wise) or right approximants (working column-wise)
 * \param[in] generic if `true`, specific code for generic instances may be used
 *
 * \todo The `row_wise` and `generic` flags are currently not supported.
 */
void approximant_basis(
                       Mat<zz_pX> & appbas,
                       VecLong & rdeg,
                       const Mat<zz_pX> & pmat,
                       const VecLong & order,
                       const VecLong & shift = VecLong(),
                       const PolMatForm form = ORD_WEAK_POPOV,
                       const bool row_wise = true,
                       const bool generic = false
                      );

/** General user-friendly interface for approximant basis computation.
 *
 * Same function as above, but taking a single positive integer for the
 * parameter `order` instead of a list. This simply interprets this as
 * specifying a list of integers all equal to `order` (with the right length),
 * and then calls the function above.
 *
 * \todo `row_wise` false not handled
 */
inline void approximant_basis(
                              Mat<zz_pX> & appbas,
                              VecLong & rdeg,
                              const Mat<zz_pX> & pmat,
                              const long order,
                              const VecLong & shift = VecLong(),
                              const PolMatForm form = ORD_WEAK_POPOV,
                              const bool row_wise = true,
                              const bool generic = false
                             )
{
    VecLong orders(pmat.NumCols(),order);
    return approximant_basis(appbas,rdeg,pmat,orders,shift,form,row_wise,generic);
}

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
 * \param[in] form required form for `appbas` (see #PolMatForm)
 * \param[in] row_wise indicates whether to compute left approximants (working row-wise) or right approximants (working column-wise)
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo add parameter row_wise
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 */
bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const VecLong & order,
                          const VecLong & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool randomized = false
                         );

/** Verifying if a matrix is a minimal approximant basis.
 *
 * Same function as above, but taking a single positive integer for the
 * parameter `order` instead of a list. This simply interprets this as
 * specifying a list of integers all equal to `order` (with the right length),
 * and then calls the function above.
 */
inline bool is_approximant_basis(
                                 const Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order,
                                 const VecLong & shift,
                                 const PolMatForm & form = ORD_WEAK_POPOV,
                                 const bool randomized = false
                                )
{
    VecLong orders(pmat.NumCols(),order);
    return is_approximant_basis(appbas,pmat,orders,shift,form,randomized);
}

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
void appbas_iterative(
                      Mat<zz_pX> & appbas,
                      const Mat<zz_pX> & pmat,
                      const VecLong & order,
                      VecLong & shift,
                      bool order_wise=true
                     );

/** Computes a `(s1,s2)`-ordered weak Popov approximant basis for
 * `([f,g]^t,order)`. This is similar to `appbas_iterative` but specific for
 * the case of an input matrix of dimensions 2 x 1. At the end of the
 * computation, the vector `(s1,s2)` contains the shifted row degree of
 * the matrix formed by `p11,p12,p21,p22`, with respect to the input shift.
 **/
// TODO try several options:
// -- rescomp vs resupdate
// -- currently it is resupdate but done naively, storing the residual
// could probably be improved (for trying to reduce shifts/copies)
void appbas_iterative_2x1(
                           zz_pX & p00,
                           zz_pX & p01,
                           zz_pX & p10,
                           zz_pX & p11,
                           const zz_pX & f0,
                           const zz_pX & f1,
                           long order,
                           long & s0,
                           long & s1
                          );

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
void popov_appbas_iterative(
                            Mat<zz_pX> & appbas,
                            const Mat<zz_pX> & pmat,
                            const VecLong & order,
                            VecLong & shift,
                            bool order_wise=true
                           );
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
 * The first set of rows is what is stored in the OUT parameter `kerbas`.
 *
 * \param[out] kerbas matrix over the base field `zz_p`
 * \param[in] pmat matrix over the base field `zz_p`
 * \param[in] shift shift (vector of integers)
 *
 * \todo This is currently not optimized when matrices with non-generic rank
 * profiles are given in input (pmat). For this, it would be much easier to
 * rely on a library which provides precisely the needed decomposition (such as
 * PLE or PLUQ), e.g. LinBox but this becomes less efficient if working over a
 * field with a prime of 30 bits or more.
 */
VecLong popov_mbasis1(
                      Mat<zz_p> & kerbas,
                      const Mat<zz_p> & pmat,
                      const VecLong & shift
                     );

/** Computes a shifted-ordered weak Popov approximant basis at order
 * `(1,...,1)` using fast linear algebra on constant matrices. Thus, in this
 * context, the input matrix `pmat` is a constant matrix. This approximant
 * basis consists of two sets of rows:
 * - rows forming a left kernel basis in (non-reduced) row echelon form of some
 *   row-permutation of `pmat` (the permutation depends on the shift)
 * - the other rows are coordinate vectors multiplied by the variable `x`;
 *   there is one such row at each index `i` which is not among the pivot
 *   indices of the left kernel basis above.
 *
 * The first set of rows is what is stored in the OUT parameter `kerbas`.
 *
 * \param[out] kerbas matrix over the base field `zz_p`
 * \param[in] pmat matrix over the base field `zz_p`
 * \param[in] shift shift (vector of integers)
 *
 * \todo This is currently not optimized when matrices with non-generic rank
 * profiles are given in input (pmat). For this, it would be much easier to
 * rely on a library which provides precisely the needed decomposition (such as
 * PLE or PLUQ), e.g. LinBox but this becomes less efficient if working over a
 * field with a prime of 30 bits or more.
 */
VecLong mbasis1(
                Mat<zz_p> & kerbas,
                const Mat<zz_p> & pmat,
                const VecLong & shift
               );

//@} // doxygen group: Approximant basis 

/** @name M-Basis algorithm (uniform approximant order)
 * \anchor mbasis
 *
 * These functions compute a `shift`-minimal ordered weak Popov approximant
 * basis for `(pmat,orders)` in the case where `orders` is given by a single
 * integer `orders = (order,...order)`. They iterate from `1` to `order`,
 * computing at each step a basis at order `1` (see @ref mbasis1) and using it
 * to update the output `appbas`, the so-called _residual matrix_, and the
 * considered shift. At step `d`, we have `appbas*pmat = 0 mod x^{d-1}`, and we
 * want to update `appbas` so that this becomes zero modulo `x^d`.
 *
 * In this context, the residual matrix is a constant matrix with the same
 * dimensions as `pmat` which, at the iteration `d`, is equal to the
 * coefficient of degree `d` of `appbas*pmat` (the coefficients of lower degree
 * being already zero).
 *
 * At the end of the computation, the vector `shift` contains the shifted row
 * degree of `appbas`, for the input shift. 
 *
 * This is inspired from the algorithm _mbasis_ described in
 *  - P. Giorgi, C.-P. Jeannerod, G. Villard. Proceeding ISSAC 2003,
 *  - P. Giorgi, R. Lebreton. Proceedings ISSAC 2014.
 */
//@{

/** Plain version of `mbasis` (see @ref mbasis), where the input `pmat` is
 * represented as `Vec<Vec<zz_pX>>`. This is almost always less efficient than
 * other provided variants, but is kept here for legacy, being a direct
 * implementation of the algorithm from the references in @ref mbasis. */
void mbasis_plain(
                  Mat<zz_pX> & appbas,
                  const Mat<zz_pX> & pmat,
                  const long order,
                  VecLong & shift
                 );

/** Variant of `mbasis` (see @ref mbasis) which first converts `pmat` to its
 * representation by a vector of constant matrices `Vec<Mat<zz_p>>`, then
 * performs the computations with this storage, and eventually converts back to
 * the polynomial matrix `Mat<zz_pX>` representation. In this variant, the
 * residual matrix is computed from `appbas` and `pmat` at each iteration.
 *
 * \todo improve when `deg(pmat) << order` 
 */
// Complexity: pmat is m x n
//   - 'order' calls to mbasis1 with dimension m x n, each one gives a
//   constant matrix K which is generically m-n x m  (may have more rows in
//   exceptional cases)
//   - order products (X Id + K ) * appbas to update the approximant basis
//   - order computations of "coeff k of appbas*pmat" to find residuals
// Assuming the degree of appbas at iteration 'ord' is m 'ord' / n (it is at
// least this almost always; and for the uniform shift it is equal to this for
// generic pmat), then the third item costs O(m n^2 order^2 / 2) operations,
// assuming cubic matrix multiplication over the field.
void mbasis_rescomp(
                    Mat<zz_pX> & appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    VecLong & shift
                   );

/** Variant of `mbasis` (see @ref mbasis) which first converts `pmat` to its
 * representation by a vector of constant matrices `Vec<Mat<zz_p>>`, then
 * performs the computations with this storage, and eventually converts back to
 * the polynomial matrix `Mat<zz_pX>` representation. In this variant, we store
 * a vector of residual matrices, initially the coefficients of `pmat`, and we
 * update all of them at each iteration; at the iteration `d` we use the `d`-th
 * matrix in this vector.
 *
 * \todo improve when `deg(pmat) << order` 
 */
// Variant which first converts to vector of constant matrices,
// performs the computations with this storage, and eventually
// converts back to polynomial matrices
// Residual (X^-d appbas*pmat mod X^(order-d)) is continuously updated along
// the iterations
// Complexity: pmat is m x n
//   - 'order' calls to mbasis1 with dimension m x n, each one gives a
//   constant matrix K which is generically m-n x m  (may have more rows in
//   exceptional cases)
//   - order products (X Id + K ) * appbas to update the approximant basis
//   - order-1 products (X Id + K ) * (matrix of degree order-ord) to update
//   the residual, for ord=1...order-1
// Assuming cubic matrix multiplication over the field, the third item costs
// O(m n (m-n) order^2/2) operations
void mbasis_resupdate(
                      Mat<zz_pX> & appbas,
                      const Mat<zz_pX> & pmat,
                      const long order,
                      VecLong & shift
                     );


/** Main `mbasis` function which chooses the most efficient variant depending
 * on the parameters (dimensions and order).
 *
 * \todo current choice is not perfectly tuned
 *    - might be wrong for a very small number of columns (but this is not the
 *    best algorithm in that case; one should rely on partial linearization,
 *    not implemented yet)
 *    - might be wrong for large orders (but one should call `pmbasis` in this case)
 *    - might not be right when the shift is far from uniform (todo: perform
 *    some tests; again, this is not the best known algorithm in this case)
 */
inline void mbasis(
                   Mat<zz_pX> & appbas,
                   const Mat<zz_pX> & pmat,
                   const long order,
                   VecLong & shift
                  )
{
    if (pmat.NumCols() > pmat.NumRows()/2 + 1)
        mbasis_resupdate(appbas, pmat, order, shift);
    else
        mbasis_rescomp(appbas, pmat, order, shift);
    // To understand the threshold (cdim > rdim/2 + 1), see the complexities
    // mentioned above for these two variants of mbasis
}

/** Computes the `shift`-Popov approximant basis for `(pmat,order)`,
 * relying on the `mbasis` approach.
 *
 * \todo currently, this is based on the approach "compute a first basis which
 * gives the pivot degrees and re-compute with these degrees as the shift to
 * obtain the normal form". It would be interesting to implement a strategy
 * which maintains a shifted Popov form along the iterations, and compare it
 * with the current technique.
 */
void popov_mbasis(
                  Mat<zz_pX> &appbas,
                  const Mat<zz_pX> & pmat,
                  const long order,
                  VecLong & shift
                 );
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

/** Computes a `shift`-ordered weak Popov approximant basis for `(pmat,order)`
 * using the algorithm PM-Basis (see @ref pmbasis) */
void pmbasis(
             Mat<zz_pX> & appbas,
             const Mat<zz_pX> & pmat,
             const long order,
             VecLong & shift
            );

/** Computes a `shift`-Popov approximant basis for `(pmat,order)` using the
 * algorithm PM-Basis (see @ref pmbasis) twice: the first call yields an
 * ordered weak Popov form which indicates the `shift`-pivot degree, which is
 * then used as a shift in the second call to obtain the sought canonical
 * basis.
 *
 * \todo often (generic case), the second call is not necessary and one just
 * has to multiply by the inverse of the leading matrix. Implement this(?).  */
void popov_pmbasis(
                   Mat<zz_pX> &appbas,
                   const Mat<zz_pX> & pmat,
                   const long order,
                   VecLong & shift
                  );

/** Computes a `(s1,s2)`-ordered weak Popov approximant basis for
 * `([f,g]^t,order)`. This is the same as `pmbasis` but specific for the case
 * of a 2x1 input matrix (see @ref pmbasis) */
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
                 long threshold=32
                );

//@} // doxygen group: PM-Basis algorithm (uniform approximant order)


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MBASIS -- GENERIC INPUT -- UNIFORM SHIFT                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name M-Basis/PM-Basis for generic input and uniform shift
 * \anchor pmbasis_generic
 *
 * The functions below take only the parameters `(appbas,pmat,order)` and they
 * compute an ordered weak Popov approximant basis `appbas` for `pmat` and the
 * uniform order `order` (here, we use the uniform shift `(0,...,0)`). It is
 * required that `pmat` is generic, in a sense detailed below; this holds with
 * high probability for a random matrix `pmat` if the field has sufficiently
 * large cardinality.
 *
 * They share the following requirements and properties:
 * - Requirement: the input matrix `pmat` is `m x n` with `m = 2*n`
 * - Requirement: the order is at least 2
 * - Requirement: the input matrix `pmat` has some genericity property: the
 *   algorithms compute kernels of constant `2n x n` matrices, which should all
 *   have the form `[ * | I ]` ; in other words the bottom `n x n` submatrix of
 *   these matrices should be invertible. (This can be rewritten as an
 *   invertibility condition on the principal minors of a block-Hankel matrix
 *   of dimension `n*order x n*order` built from `pmat`.)
 * - Property: the output basis `appbas` is in ordered weak Popov form with
 *   degree `ceil(order/2)`. More precisely:
 *     - if `order = 2d`, then `appbas = [[X^d I + P00,  P01], [X P10, X^d I +
 *     X P11]]` where `P00`, `P01`, `P10` have degree `d-1` and `P11` has
 *     degree `d-2`
 *     - if `order = 2d+1`, then `appbas = [[X^{d+1} I + X P00,  X P01], [P10,
 *     X^d I + P11]]` where `P00`, `P01`, `P11` have degree `d-1` and `P10` has
 *     degree `d`
 *   In particular, in both cases, the leading and trailing principal `n x n`
 *   submatrices of appbas are in Popov form.
 *
 * \todo Currently, no check is performed concerning these requirements and the
 * algorithm may throw an error or have undefined behavior if `pmat` does not
 * have the required genericity property.
 */
//@{

/** Computes an ordered weak Popov approximant basis `appbas` for
 * `(pmat,order)` iteratively, similarly to #mbasis_rescomp, under the
 * requirements and with the properties detailed in @ref pmbasis_generic.
 *
 * \todo throw away if never faster than the _resupdate_ version.
 */
void mbasis_generic_2n_n_rescomp(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order
                                );

/** Computes an ordered weak Popov approximant basis `appbas` for
 * `(pmat,order)` iteratively, similarly to #mbasis_resupdate, under the
 * requirements and with the properties detailed in @ref pmbasis_generic. */
void mbasis_generic_2n_n_resupdate(
                                   Mat<zz_pX> & appbas,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );

/** Computes an ordered weak Popov approximant basis `appbas` for
 * `(pmat,order)` by a divide and conquer algorithm, similar to #pmbasis, under
 * the requirements and with the properties detailed in @ref pmbasis_generic.
 *
 * Considering the two forms of bases described in @ref pmbasis_generic, which
 * we call form 1 and form 2 respectively, we note that:
 *  - the product of two bases having form 1 (with respective degrees `d1` and
 *  `d2`) also has form 1 (with degree `d1 + d2`)
 *  - the product of a basis having form 2 (with degree `d1`) by one having
 *  form 1 (with degree `d2`) is a basis with form 2 (and degree `d1+d2`)
 *
 * We use this remark to choose specific orders `order1` and `order2` close to
 * `order/2` for the recursive calls, allowing us to never have to deal with
 * shifts (which usually appear during the computations even when one starts
 * with the uniform shift `(0,...,0)`):
 *   - if `order` is even, all orders of recursive calls are even (as a result,
 *   the final basis is a product of forms 1, and thus has form 1 itself)
 *   - if `order` is odd, then only the very last leaf of the recursive tree
 *   will be with odd order, the others will be with even order; since the last
 *   leaf gives the leftmost basis in the product which yields the final basis,
 *   this means all bases will have form 1 except the leftmost one which has
 *   form 2, hence the final one has form 2.
 * */
void pmbasis_generic_2n_n(
                          Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order
                         );

/** Computes the top half `appbas` of the rows of an ordered weak Popov
 * approximant basis for `(pmat,order)`, under the requirements and with the
 * properties detailed in @ref pmbasis_generic. Fundamentally relies on
 * #pmbasis_generic_2n_n, with a slightly modified recursion tree on top (the
 * branch from the root to the last leaf computes only the top rows, the others
 * are unchanged). */
void pmbasis_generic_2n_n_top_rows(
                                   Mat<zz_pX> & appbas,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );

//@} // doxygen group: M-Basis/PM-Basis for generic input and uniform shift

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MATRIX-PADE APPROXIMATION -- GENERIC INPUT -- no shift     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Matrix-Padé approximation
 * \anchor matrix_pade
 *
 * _Example use_: matrix fraction reconstruction in the context of
 * block-Wiedemann-like algorithms.
 *
 * The problem is the following:
 * - Input: a square `n x n` matrix `pmat` of degree less than `order`,
 * - Output: a square `n x n` matrix `den` of degree at most `order/2` and in
 *   row reduced form such that `num = den * pmat` has degree less than
 *   `order/2`.
 *
 * Here, `den` stands for denominator, `num` stands for numerator, and if
 * the order is large enough, `pmat = den^{-1} num` is a proper irreducible
 * left fraction description.
 *
 * Note that this can be solved by taking the leading principal `n x n`
 * submatrix of an ordered weak Popov approximant basis for `[[pmat], [-Id]]`
 * at order `order`. Here we aim at deriving algorithms close to `mbasis` and
 * `pmbasis` but which exploit the fact that there is this identity matrix in
 * the input and that we only want a submatrix of the final approximant basis.
 *
 * \todo implement general Matrix-Pade (for the moment, only the generic case
 * is implemented).
 *
 * The functions below share the following requirements:
 * - Requirement: the input matrix `pmat` is square, `n x n`
 * - Requirement: the input matrix `pmat` has some genericity property: the
 *   algorithms compute kernels of constant `2n x n` matrices, which should all
 *   have the form `[ * | I ]` ; in other words the bottom `n x n` submatrix of
 *   these matrices should be invertible. (This can be rewritten as an
 *   invertibility condition on the principal minors of a block-Hankel matrix
 *   of dimension `n*order x n*order` built from `pmat`.)
 *
 * The latter holds with high probability for a random matrix `pmat`, if the
 * field is sufficiently large. Note that here no check is performed and the
 * algorithm may throw an error if `pmat` does not have the required genericity
 * property.
 * 
 * Some of these functions return only the denominator `den`, in Popov form.
 * Precisely it has the form `X^dd I + D`, where `deg(D) < dd` and `dd =
 * ceil(order/2)`, and it is exactly the leading principal `n x n` submatrix of
 * the approximant basis computed by the functions #pmbasis_generic_2n_n above,
 * on input `order` and `[[pmat], [-I]]`.
 *
 * The other functions return a couple `(den1,den2)`, where:
 *    - `den1` is as in the previous paragraph,
 *    - `den2` is a matrix of degree `floor(order/2)` such that `[[den1],
 *    [den2]]` forms the left `2n x n` submatrix of the approximant basis
 *    computed by the functions `mbasis_generic_2n_n` above, on input `order` and
 *    `[[pmat], [-I]]` (in particular, if `order` is even then `den2(0)` is zero).
 *  
 */
//@{

/** Computes the denominators `(den1,den2)` for `(pmat,order)` iteratively,
 * with the properties and under the requirements detailed in @ref
 * matrix_pade. This algorithm may be seen as Berlekamp-Massey performed on
 * square matrices, and reversed (it goes from low degree to high degree).
 * */
void matrix_pade_generic_iterative(
                                   Mat<zz_pX> & den1,
                                   Mat<zz_pX> & den2,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );

/** Computes the denominator `den` for `(pmat,order)` iteratively, with the
 * properties and under the requirements detailed in @ref matrix_pade. This
 * simply calls the function of the same name which computes both `den1` and
 * `den2`, and then discards `den2`. */
inline void matrix_pade_generic_iterative(
                                          Mat<zz_pX> & den,
                                          const Mat<zz_pX> & pmat,
                                          const long order
                                         )
{
    Mat<zz_pX> den2;
    matrix_pade_generic_iterative(den, den2, pmat, order);
}

/** Computes the denominators `(den1,den2)` for `(pmat,order)` by a divide and
 * conquer approach (the matrix `den` is `2n x n` and stores `[[den1],
 * [den2]]`), with the properties and under the requirements detailed in @ref
 * matrix_pade. This is essentially a specialization of `pmbasis` to the
 * present situation, where we need just the left half of the columns of the
 * output approximant basis, and where the input contains an identity matrix.
 */
void matrix_pade_generic_recursion(
                                   Mat<zz_pX> & den,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );

/** Computes the denominator `den` for `(pmat,order)` by a divide and conquer
 * algorithm, with the properties and under the requirements detailed in @ref
 * matrix_pade. This fundamentally relies on #matrix_pade_generic_recursion
 * with the recursion tree slightly modified with the branch from the root to
 * the last leaf computing just the top denominator instead of both `den1` and
 * `den2`, thanks to a call to #pmbasis_generic_2n_n_top_rows. */
void matrix_pade_generic(
                         Mat<zz_pX> & den,
                         const Mat<zz_pX> & pmat,
                         const long order
                        );

//@} // doxygen group: Matrix-Padé approximation



/*------------------------------------------------------------*/
/* TODO in progress: (P)MBASIS, generic case, few columns     */
/*------------------------------------------------------------*/

/** Using partial linearization to handle the case where the number of columns
 * is quite smaller than the number of rows. Iterative.
 *
 * \todo Currently a prototype. Not tuned and barely tested. */
VecLong mbasis_generic_onecolumn(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order,
                                 const VecLong & shift
                                );

/** Using partial linearization to handle the case where the number of columns
 * is quite smaller than the number of rows. Divide and conquer.
 *
 * \todo Currently a prototype. Not tuned and barely tested. */
VecLong pmbasis_generic_onecolumn(
                                  Mat<zz_pX> & appbas,
                                  const Mat<zz_pX> & pmat,
                                  const long order,
                                  const VecLong & shift
                                 );


/** \todo doc */
void popov_mbasis_rescomp_generic(
                                  Mat<zz_pX> & appbas,
                                  const Mat<zz_pX> & pmat,
                                  const long order
                                 );

void popov_pmbasis_generic(
                           Mat<zz_pX> &appbas,
                           const Mat<zz_pX> & pmat,
                           const long order
                          );



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TODO                                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// ** pmbasis: handle different orders via max(order)-order shifting of pmat

// ** generic (uniform shift): just do everything without shifts as would
// happen generically with the uniform shift; might need assumptions on n
// divide m, etc  (--> does not look very promising, from the 2nxn case)

// ** shifted generic: from shift, deduce the pivot degree expected generically
// and use this as a shift instead of 'shift', then obtain directly Popov
// approx basis; skipping the first call to find the pivot degree (in addition,
// will this be more efficient because shift is nicer?)

// ** thresholds:
//   -- mbasis_resupdate / mbasis_rescomp depending on m/n and nthreads
//   -- in pmbasis: base case

// ** threads:
//   -- in mbasis_resupdate
//   -- in pmbasis?

// iterative version / mbasis: compare different ways of obtaining Popov
// (recompute with new shift, or maintain normal form)


#endif // MAT_LZZ_PX_APPROXIMANT__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
