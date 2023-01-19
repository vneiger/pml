/**
 * \file nmod_mat_poly.h
 * \brief Matrices with univariate polynomial entries modulo word-size prime, represented as list of matrices
 * \version 0.0
 * \date 2022-06-25
 *
 * \todo doc
 *
 */

#ifndef NMOD_POLY_MAT_MAT_POLY_H
#define NMOD_POLY_MAT_MAT_POLY_H

#include <flint/perm.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    slong degree;
    slong length;
    nmod_mat_struct *mat;
    slong r;
    slong c;
    mp_limb_t mod;
} nmod_mat_poly_struct;
// TODO use length and alloc as for Flint's poly

typedef nmod_mat_poly_struct nmod_mat_poly_t[1];

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* GETTERS                                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

NMOD_POLY_MAT_INLINE slong
nmod_mat_poly_nrows(const nmod_mat_poly_t matp)
{
    return matp->r;
}

NMOD_POLY_MAT_INLINE slong
nmod_mat_poly_ncols(const nmod_mat_poly_t matp)
{
	return matp->c;
}

NMOD_POLY_MAT_INLINE mp_limb_t
nmod_mat_poly_modulus(const nmod_mat_poly_t matp)
{
	return matp->mod;
}

void nmod_mat_poly_get_coef(nmod_mat_t res,
                            const nmod_mat_poly_t F,
                            slong k);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SETTERS                                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


// TODO getter for pointer to i-th coeff

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MEMORY MANAGEMENT                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO realloc?
// TODO fit_length?

void nmod_mat_poly_clear(nmod_mat_poly_t A);

void nmod_mat_poly_init(nmod_mat_poly_t res, slong degree,
                        slong length,
                        slong r, slong c, mp_limb_t mod);


/** void nmod_mat_poly_init_set(nmod_mat_poly_t res,
 *                              const nmod_poly_mat_t F)
 *
 * F \in K[x]^{mxn} <-> res \in K^{mxn}[x]
 * With res->degree = F.degree(), and res will have only res->degree + 1 nmod_mat_t pointer
 *
 */
void nmod_mat_poly_init_set(nmod_mat_poly_t res,
                            const nmod_poly_mat_t F);

/** void nmod_mat_poly_init_set(nmod_mat_poly_t res,
 *				                const nmod_poly_mat_t F)
 *
 * F \in K[x]^{mxn} <-> res \in K^{mxn}[x]
 * With res->degree = F.degree(), and res will have only length nmod_mat_t pointer
 * Improves memory space for x^{-k} P_{k-1} F mod x in M_basis{III,IV}
 *
 */
void nmod_mat_poly_init_setII(nmod_mat_poly_t res,
                              const nmod_poly_mat_t F,
                              slong length);


/** void nmod_mat_poly_init_set(nmod_mat_poly_t res,
 *				                const nmod_poly_mat_t F)
 *
 * F \in K[x]^{mxn} <-> res \in K^{mxn}[x]
 * With res->degree = F.degree(), and res will have res->degree + length nmod_mat_t pointer
 * To store the computation of P_{k-1} F in M_basisV
 *
 */
void nmod_mat_poly_init_setIII(nmod_mat_poly_t res,
                               const nmod_poly_mat_t F,
                               slong length);



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TODO to categorize                                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/** void nmod_mat_poly_naive_mul_coef(nmod_mat_t res,
 *				       const nmod_mat_poly_t A,
 *				       const nmod_mat_poly_t B, slong k);
 *
 * A = sum^{deg_A}_{i=0} a_i x^i and B = sum^{deg_B}_{i=0} b_i x^i
 * Compute the coefficient k of the product AB
 * C = AB = sum^{deg_A + deg_B}_{i=0} c_i x^i
 * res = c_k = sum^{k}_{i=0} a_i b_{k-i}
 *
 */
void nmod_mat_poly_naive_mul_coef(nmod_mat_t res,
                                  const nmod_mat_poly_t A,
                                  const nmod_mat_poly_t B,
                                  slong k);

void nmod_mat_poly_print(const nmod_mat_poly_t A);

#ifdef __cplusplus
}
#endif

#endif /* NMOD_POLY_MAT_MAT_POLY_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
