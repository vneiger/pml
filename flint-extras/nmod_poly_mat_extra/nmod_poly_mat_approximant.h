/**
 * \file nmod_poly_mat_approximant.h
 * \brief Minimal approximant bases
 * \author Vincent Neiger, Kevin Tran
 * \version 0.0
 * \date 2022-06-25
 *
 *
 */

#ifndef NMOD_POLY_MAT_APPROXIMANT_H
#define NMOD_POLY_MAT_APPROXIMANT_H

// TODO [urgent] change name
#define INITIAL_CASE_BOUND 32

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h> 

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \fn void Basis(nmod_poly_mat_t res, int64_t *res_shift,
	   nmod_mat_t mat, int64_t *shift);
 * \brief set on res the minimal approximant basis of mat for order 1
 and set res_shift the shift row degree of res
 * 
 *
 * \param All parameters are supposed init
 * \param res a polynomial matrix of dimensions (rdim, rdim)
 * \param res_shift a vector of length rdim
 * \param mat a polynomial matrix
 * \param shift a vector of length rdim
 *
 */
void Basis(nmod_poly_mat_t res, int64_t *res_shifts,
	   const nmod_mat_t mat, const int64_t *shifts);

slong Basis_for_M_basis(nmod_mat_t res, int64_t *res_shifts, slong *res_perm,
			const nmod_mat_t mat, const int64_t *shifts);


/**
 * \brief Do the multiplication permutation matrix (P*A) for the function Basis
 *
 * \param A, the input and result matrix
 * \param P, the permutation 
 * \param n, length of P and row dimension of A
 */
void apply_perm_rows_to_matrix(nmod_mat_t mat, const slong * perm, slong rdim);

void apply_perm_cols_to_poly_matrix(nmod_poly_mat_t mat, const slong * perm, slong cdim);

void apply_perm_rows_to_poly_matrix(nmod_poly_mat_t mat, const slong *perm, slong rdim);

void apply_perm_to_vector(int64_t *res, const int64_t *initial_vect,
			  const slong *perm, slong length);

void sort_and_create_perm(slong *perm, const int64_t *vec, slong n);















/** static void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
 *                                               const slong *perm, slong rank)
 * 
 * This function compute the multiplication of specific polynomials matrix
 * It takes an intrger r, A a nmod_mat_t, res a nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * M = perm^(-1) * [[x, 0], [A, 1]] * perm and res = [[R1],[R2]] 
 * Stocks in res the result
 *
 */
void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
				      const slong *perm, slong rank);

/** void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * It used the structured multiplication blocks to compute x^{-k} P_{k-1} F mod x and to compute
 * P_{k} = M P_{k-1}
 * Use polynomial matrix multiplication
 *
 */
void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	     const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);


/** void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] will be FIX and
 * compute iteratively P_{k} \in K[x]^{mxm}, then will transform to P_prime_{k} \in K^{mxm}[x] 
 * and compute x^{-k} P_{k-1} F mod x with a naive polynomial multiplication
 *
 * Use naive polynomial multiplication
 */
void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);

/** void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_list_poly_mat.h)
 * Compute P_{k-1} \in K[x]^{mxm} with the structured_multiplication_blocks
 *
 * Use structured_multiplication_blocks and list_structured_multiplication_blocks 
 */
void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
		const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);

/** void M_basisIV(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_list_poly_mat.h)
 * Compute P_{k-1} \in K^{mxm}[x] with the list_structured_multiplication_blocks
 * Transform P_{sigma - 1} as a poly_mat_t 
 *
 * Use list_structured_multiplication_blocks 
 */
void M_basisIV(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);


/** void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] FIX and
 * Compute P_{k-1} \in K[x]^{mxm} with the list structured_multiplication_blocks
 * Compute x^{-k} P_{k-1} F mod x iteratively with a naive polynomial multiplication  
 *
 * Use naive polynomial multiplication and list_structured_multiplication_blocks 
 */
void M_basisV(nmod_poly_mat_t res, int64_t *res_shifts,
	      const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);















void PM_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	      const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts);

















#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_APPROXIMANT_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
