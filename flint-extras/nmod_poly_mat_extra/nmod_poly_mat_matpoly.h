/**
 * \file nmod_mat_poly.h
 * \brief Matrices with univariate polynomial entries modulo word-size prime, represented as list of matrices
 * \author Vincent Neiger, Kevin Tran
 * \version 0.0
 * \date 2022-06-25
 *
 * \todo doc
 *
 */

#ifndef NMOD_MAT_POLY_H
#define NMOD_MAT_POLY_H

#include <flint/perm.h>
#include <flint/nmod_poly_mat.h>

typedef struct
{
  slong degree;
  slong length;
  nmod_mat_struct *mat;
  slong rdim;
  slong cdim;
  mp_limb_t modulus;
} nmod_list_poly_mat_struct;

typedef nmod_list_poly_mat_struct nmod_list_poly_mat_t[1];

slong nmod_list_poly_mat_nrows(const nmod_list_poly_mat_t A);

slong nmod_list_poly_mat_ncols(const nmod_list_poly_mat_t A);

slong nmod_list_poly_mat_degree(const nmod_list_poly_mat_t A);

mp_limb_t nmod_list_poly_mat_modulus(const nmod_list_poly_mat_t A);

void nmod_list_poly_mat_print(const nmod_list_poly_mat_t A);

void nmod_list_poly_mat_get_coef(nmod_mat_t res, const nmod_list_poly_mat_t F,
				 slong k);

void nmod_list_poly_mat_init(nmod_list_poly_mat_t res, slong degree,
			     slong length,
			     slong rdim, slong cdim, mp_limb_t modulus);

void nmod_list_poly_mat_clear(nmod_list_poly_mat_t A);


/** void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res,
 *				 const nmod_poly_mat_t F)
 *
 * F \in K[x]^{mxn} <-> res \in K^{mxn}[x]
 * With res->degree = F.degree(), and res will have only res->degree + 1 nmod_mat_t pointer
 * 
 */
void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res,
				 const nmod_poly_mat_t F);

/** void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res,
 *				 const nmod_poly_mat_t F)
 *
 * F \in K[x]^{mxn} <-> res \in K^{mxn}[x]
 * With res->degree = F.degree(), and res will have only length nmod_mat_t pointer
 * Improve for the memory space of the computation of x^{-k} P_{k-1} F mod x in M_basisIII and IV 
 *
 */
void nmod_list_poly_mat_init_setII(nmod_list_poly_mat_t res,
				   const nmod_poly_mat_t F, slong length);


/** void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res,
 *				 const nmod_poly_mat_t F)
 *
 * F \in K[x]^{mxn} <-> res \in K^{mxn}[x]
 * With res->degree = F.degree(), and res will have res->degree + length nmod_mat_t pointer
 * To stock the computation of P_{k-1} F in M_basisV 
 *
 */
void nmod_list_poly_mat_init_setIII(nmod_list_poly_mat_t res,
				    const nmod_poly_mat_t F, slong length);

void nmod_list_poly_mat_set(nmod_list_poly_mat_t res,
			    const nmod_poly_mat_t F);

void nmod_list_poly_mat_to_poly_mat(nmod_poly_mat_t res,
				    const nmod_list_poly_mat_t F);


/** void nmod_list_poly_mat_naive_mul_coef(nmod_mat_t res,
 *				       const nmod_list_poly_mat_t A,
 *				       const nmod_list_poly_mat_t B, slong k);
 *
 * A = sum^{deg_A}_{i=0} a_i x^i and B = sum^{deg_B}_{i=0} b_i x^i
 * Compute the coefficient k of the product AB 
 * C = AB = sum^{deg_A + deg_B}_{i=0} c_i x^i
 * res = c_k = sum^{k}_{i=0} a_i b_{k-i}
 *
 */
void nmod_list_poly_mat_naive_mul_coef(nmod_mat_t res,
				       const nmod_list_poly_mat_t A,
				       const nmod_list_poly_mat_t B, slong k);

/** static void list_structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
 *                                               const slong *perm, slong rank,
 *                                               slong k, slong sigma)
 * 
 * This function compute the multiplication of specific polynomials matrix
 * A a nmod_mat_t, res a list_nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * P = perm^(-1) * [[x, 0], [A, 1]] * perm \in K[x]^{mxm}  and 
 * res = sum_{i=0}^{sigma - 1} r_i x^i \in K^{mxn}[x]
 * To do that we will permutate r_i's rows, shift the top and mul by A and add the bottom on each i
 * It will not take in count the k - 2 first matrix of res, 
 * and the degree(res) - sigma last matrix 
 * because we will only need the k-th coefficients of res
 * and k will go to sigma - 1 => we need to update the sigma - 1 coefficients
 *
 */
void structured_list_multiplication_blocks(nmod_list_poly_mat_t res,
					   const nmod_mat_t A,
					   const slong *perm,
					   slong rank, slong k, slong sigma);

void nmod_list_poly_mat_to_poly_mat(nmod_poly_mat_t res,
				    const nmod_list_poly_mat_t F);


/** static void list_structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
 *                                               const slong *perm, slong rank,
 *                                               slong k, slong sigma)
 * 
 * This function compute the multiplication of specific polynomials matrix
 * A a nmod_mat_t, res a list_nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * P = perm^(-1) * [[x, 0], [A, 1]] * perm \in K[x]^{mxm}  and 
 * res = sum_{i=0}^{degree(res)} r_i x^i \in K^{mxn}[x]
 * To do that we will permutate r_i's rows, shift the top and mul by A and add the bottom on each i
 * and update the degree of res
 */
void structured_list_multiplication_blocks_full(nmod_list_poly_mat_t res,
						const nmod_mat_t A,
						const slong *perm, slong rank);

#endif /* NMOD_LIST_POLY_MAT_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
