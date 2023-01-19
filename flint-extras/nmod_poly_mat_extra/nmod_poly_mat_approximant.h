/**
 * \file nmod_poly_mat_approximant.h
 * \brief Minimal approximant bases
 * \version 0.0
 * \date 2022-06-25
 *
 */

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
 * TODO rather use orientation
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
 * \param[in] form indicates the required form for the output basis (see #poly_mat_form_t)
 * \param[in] row_wise indicates whether to compute left approximants (working row-wise) or right approximants (working column-wise)
 * \param[in] generic if `true`, specific code for generic instances may be used
 *
 * \todo The `row_wise` and `generic` flags are currently not supported.
 * \todo
 */
//void approximant_basis(
//                       Mat<zz_pX> & appbas,
//                       VecLong & rdeg,
//                       const Mat<zz_pX> & pmat,
//                       const VecLong & order,
//                       const VecLong & shift = VecLong(),
//                       const poly_mat_form_t form = ORD_WEAK_POPOV,
//                       const bool row_wise = true,
//                       const bool generic = false
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
//inline void approximant_basis(
//                              Mat<zz_pX> & appbas,
//                              VecLong & rdeg,
//                              const Mat<zz_pX> & pmat,
//                              const long order,
//                              const VecLong & shift = VecLong(),
//                              const poly_mat_form_t form = ORD_WEAK_POPOV,
//                              const bool row_wise = true,
//                              const bool generic = false
//                             )
//{
//    VecLong orders(pmat.NumCols(),order);
//    return approximant_basis(appbas,rdeg,pmat,orders,shift,form,row_wise,generic);
//}

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
 * \todo
 */
//bool is_approximant_basis(
//                          const Mat<zz_pX> & appbas,
//                          const Mat<zz_pX> & pmat,
//                          const VecLong & order,
//                          const VecLong & shift,
//                          const poly_mat_form_t & form = ORD_WEAK_POPOV,
//                          const bool randomized = false
//                         );

/** Verifying if a matrix is a minimal approximant basis.
 *
 * Same function as above, but taking a single positive integer for the
 * parameter `order` instead of a list. This simply interprets this as
 * specifying a list of integers all equal to `order` (with the right length),
 * and then calls the function above.
 * \todo
 */
//inline bool is_approximant_basis(
//                                 const Mat<zz_pX> & appbas,
//                                 const Mat<zz_pX> & pmat,
//                                 const long order,
//                                 const VecLong & shift,
//                                 const poly_mat_form_t & form = ORD_WEAK_POPOV,
//                                 const bool randomized = false
//                                )
//{
//    VecLong orders(pmat.NumCols(),order);
//    return is_approximant_basis(appbas,pmat,orders,shift,form,randomized);
//}

/** TODO temporary fix; and not finished yet **/
int is_minimal_approximant_basis(const nmod_poly_mat_t base,
                                 const nmod_mat_t mat,
                                 slong order,
                                 const slong * shifts);

//@} // doxygen group: General interfaces for approximant basis computation and verification






















/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MBASIS1                                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/**
 * \fn void Basis(nmod_poly_mat_t res, slong * res_shift, nmod_mat_t mat, slong * shift)
 *
 * \brief store in res the shifts-minimal approximant basis of mat at order 1,
 * store in res_shift the shifts-row degree of res
 * TODO if output Popov, say it
 * 
 * All parameters are supposed initialized.
 *
 * \param res a polynomial matrix of dimensions (rdim, rdim)
 * \param res_shift a vector of length rdim
 * \param mat a polynomial matrix
 * \param shift a vector of length rdim
 *
 */
void Basis(nmod_poly_mat_t res,
           slong * res_shifts,
           const nmod_mat_t mat,
           const slong * shifts);

slong Basis_for_M_basis(nmod_mat_t res,
                        slong * res_shifts,
                        slong * res_perm,
                        const nmod_mat_t mat,
                        const slong * shifts);

/** static void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
 *                                               const slong *perm, slong rank)
 * 
 * This function computes the multiplication of specific polynomial matrix
 * It takes an integer r, A a nmod_mat_t, res a nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * M = perm^(-1) * [[x, 0], [A, 1]] * perm and res = [[R1],[R2]] 
 * Stores the result in res
 */
void structured_multiplication_blocks(nmod_poly_mat_t res,
                                      const nmod_mat_t A,
                                      const slong * perm,
                                      slong rank);

/** void M_basis(nmod_poly_mat_t res, slong * res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const slong * shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * It used the structured multiplication blocks to compute x^{-k} P_{k-1} F mod x and to compute
 * P_{k} = M P_{k-1}
 * Use polynomial matrix multiplication
 *
 */
void M_basis(nmod_poly_mat_t res,
             slong * res_shifts,
             const nmod_poly_mat_t F,
             ulong sigma,
             const slong * shifts);


/** void M_basisII(nmod_poly_mat_t res, slong * res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const slong * shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] will be FIX and
 * compute iteratively P_{k} \in K[x]^{mxm}, then will transform to P_prime_{k} \in K^{mxm}[x] 
 * and compute x^{-k} P_{k-1} F mod x with a naive polynomial multiplication
 *
 * Use naive polynomial multiplication
 */
void M_basisII(nmod_poly_mat_t res,
               slong * res_shifts,
               const nmod_poly_mat_t F,
               ulong sigma,
               const slong * shifts);

/** void M_basisIII(nmod_poly_mat_t res, slong * res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const slong * shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_mat_poly.h)
 * Compute P_{k-1} \in K[x]^{mxm} with the structured_multiplication_blocks
 *
 * Use structured_multiplication_blocks and list_structured_multiplication_blocks 
 */
void M_basisIII(nmod_poly_mat_t res,
                slong * res_shifts,
                const nmod_poly_mat_t F,
                ulong sigma,
                const slong * shifts);

/** void M_basisIV(nmod_poly_mat_t res, slong * res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const slong * shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_mat_poly.h)
 * Compute P_{k-1} \in K^{mxm}[x] with the list_structured_multiplication_blocks
 * Transform P_{sigma - 1} as a poly_mat_t 
 *
 * Use list_structured_multiplication_blocks 
 */
void M_basisIV(nmod_poly_mat_t res,
               slong * res_shifts,
               const nmod_poly_mat_t F,
               ulong sigma,
               const slong * shifts);


/** void M_basisIII(nmod_poly_mat_t res, slong * res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const slong * shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] FIX and
 * Compute P_{k-1} \in K[x]^{mxm} with the list structured_multiplication_blocks
 * Compute x^{-k} P_{k-1} F mod x iteratively with a naive polynomial multiplication  
 *
 * Use naive polynomial multiplication and list_structured_multiplication_blocks 
 */
void M_basisV(nmod_poly_mat_t res,
              slong * res_shifts,
              const nmod_poly_mat_t F,
              ulong sigma,
              const slong * shifts);


/** void PM_basis(nmod_poly_mat_t res, slong * res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const slong * shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] FIX and
 * Compute P_{k-1} \in K[x]^{mxm} with the list structured_multiplication_blocks
 * Compute x^{-k} P_{k-1} F mod x iteratively with a naive polynomial multiplication  
 *
 * Use naive polynomial multiplication and list_structured_multiplication_blocks 
 */
void PM_basis(nmod_poly_mat_t res,
              slong * res_shifts,
              const nmod_poly_mat_t F,
              ulong sigma,
              const slong * shifts);


#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_APPROXIMANT_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
