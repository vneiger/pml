#include "nmod_mat_extra.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_mat_poly.h" // TODO remove
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

/****************************************************************
*  TODO next functions for structured blocks mul to be cleaned  *
****************************************************************/
/** 
 * This function computes the multiplication of specific polynomial matrix
 * It takes rank, mat, pmat, perm. 
 * It will compute the mutiplication of
 * perm^(-1) * [[x, 0], [mat, 1]] * perm and pmat = [[pmat1],[pmat2]] 
 * Stores the result in pmat
 */
void structured_multiplication_blocks(nmod_poly_mat_t pmat,
                                      const nmod_mat_t mat,
                                      const slong * perm,
                                      slong rank);



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
void structured_list_multiplication_blocks(nmod_mat_poly0_t res,
					   const nmod_mat_t A,
					   const slong *perm,
					   slong rank, slong deb, slong sigma)
{
	slong r = res->r, c = res->c;
	mp_limb_t mod = res->mod;
	slong i;
	nmod_mat_t previous_R1, R1, R2, R1_cp, R2_cp;
	nmod_mat_struct *r_i;
	slong *inv_perm = _perm_init(r);

	/** init **/
	nmod_mat_init(R1_cp, rank, c, mod);
	nmod_mat_init(previous_R1, rank, c, mod);
	nmod_mat_init(R2_cp, r - rank, c, mod);

	for (i = deb; i < sigma; i++)
	{
		r_i = res->mat + i;
		nmod_mat_permute_rows(r_i, perm, NULL);
		nmod_mat_window_init(R1, r_i, 0, 0, rank, c);
		nmod_mat_set(R1_cp, R1);

		nmod_mat_set(R1, previous_R1); //it will set the top block of r_{i-1} on r_{i}

	nmod_mat_set(previous_R1, R1_cp);

	nmod_mat_window_init(R2, r_i, rank, 0, r, c);
	nmod_mat_set(R2_cp, R2);

	nmod_mat_mul(R2, A, R1_cp);
	nmod_mat_add(R2, R2, R2_cp);

	nmod_mat_window_clear(R1);
	nmod_mat_window_clear(R2);
	}

	/** apply perm^(-1) **/
	_perm_inv(inv_perm, perm, r);
	for (i = deb + 1; i < sigma; i++)
		nmod_mat_permute_rows(res->mat + i, inv_perm, NULL);

	/** clear **/
	nmod_mat_clear(R1_cp);
	nmod_mat_clear(R2_cp);
	nmod_mat_clear(previous_R1);

	_perm_clear(inv_perm);
}

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
void structured_list_multiplication_blocks_full(nmod_mat_poly0_t res,
						const nmod_mat_t A,
						const slong *perm, slong rank)
{
	slong r = res->r, c = res->c;
	mp_limb_t mod = res->mod;
	slong i;
	nmod_mat_t previous_R1, R1, R2, R1_cp, R2_cp;
	nmod_mat_struct *r_i;
	slong *inv_perm = _perm_init(r);

	/** init **/
	nmod_mat_init(R1_cp, rank, c, mod);
	nmod_mat_init(previous_R1, rank, c, mod);
	nmod_mat_init(R2_cp, r - rank, c, mod);

	res->degree += 1;

	for (i = 0; i <= res->degree; i++)
	{
		r_i = res->mat + i;
		nmod_mat_permute_rows(r_i, perm, NULL);
		nmod_mat_window_init(R1, r_i, 0, 0, rank, c);
		nmod_mat_set(R1_cp, R1);

		nmod_mat_set(R1, previous_R1); //it will set the top block of r_{i-1} on r_{i}

	nmod_mat_set(previous_R1, R1_cp);

	nmod_mat_window_init(R2, r_i, rank, 0, r, c);
	nmod_mat_set(R2_cp, R2);

	nmod_mat_mul(R2, A, R1_cp);
	nmod_mat_add(R2, R2, R2_cp);

	nmod_mat_window_clear(R1);
	nmod_mat_window_clear(R2);
	}

	/** apply perm^(-1) **/
	_perm_inv(inv_perm, perm, r);
	for (i = 0; i <= res->degree; i++)
		nmod_mat_permute_rows(res->mat + i, inv_perm, NULL);

	/** clear **/
	nmod_mat_clear(R1_cp);
	nmod_mat_clear(R2_cp);
	nmod_mat_clear(previous_R1);

	_perm_clear(inv_perm);
}

void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
                                      const slong *perm, slong rank)
{
    slong alloc, rdim = res->r, cdim = res->c, prime = res->modulus;
    nmod_poly_mat_t R1, R2, R1_cp, R2_cp, A_poly;
    nmod_poly_t P;

    /** init **/
    slong *inv_perm = _perm_init(rdim);
    nmod_poly_init(P, prime);
    nmod_poly_mat_init(R1_cp, rank, cdim, prime);
    nmod_poly_mat_init(R2_cp, rdim - rank, cdim, prime);

    /** apply perm **/
    nmod_poly_mat_permute_rows(res, perm, NULL);

    /** work on R1 (top) **/
    nmod_poly_mat_window_init(R1, res, 0, 0, rank, cdim);
    nmod_poly_mat_set(R1_cp, R1);

    nmod_poly_mat_shift_left(R1, R1, 1); //matpol function

    /** work on R2 (bottom) **/
    nmod_poly_mat_window_init(R2, res, rank, 0, rdim, cdim);
    nmod_poly_mat_set(R2_cp, R2);

    nmod_poly_mat_init(A_poly, rdim - rank, rank, prime);
    for (slong i = 0; i < rdim - rank; i++)
        for (slong j = 0; j < rank; j++)
        {
            alloc = (slong) nmod_mat_get_entry(A, i, j);
            if (alloc != 0)
            {
                nmod_poly_set_coeff_ui(P, 0, alloc);
                nmod_poly_set(nmod_poly_mat_entry(A_poly, i, j), P);
            }
        }

    nmod_poly_mat_mul(R2, A_poly, R1_cp);
    nmod_poly_mat_add(R2, R2, R2_cp);

    /** apply perm^(-1) **/
    _perm_inv(inv_perm, perm, rdim);
    nmod_poly_mat_permute_rows(res, inv_perm, NULL);

    /** clear **/
    nmod_poly_clear(P);

    nmod_poly_mat_clear(A_poly);
    nmod_poly_mat_clear(R1_cp);
    nmod_poly_mat_clear(R2_cp);

    nmod_poly_mat_window_clear(R1);
    nmod_poly_mat_window_clear(R2);

    _perm_clear(inv_perm);
}

void mbasis(nmod_poly_mat_t res, slong *res_shifts,
             const nmod_poly_mat_t F, ulong sigma, const slong *shifts)

{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;

    nmod_poly_mat_t F_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    /** init **/
    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_poly_mat_init_set(F_prime, F);
    nmod_poly_mat_get_coeff_mat(constant_mat, F_prime, 0); //Compute F mod x

    rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank); //Compute P0
    for (ulong k = 1; k < sigma; k++)
    {
        // doing the operation x^(-k) F_prime mod x
        structured_multiplication_blocks(F_prime, A_k, perm, rank);
        nmod_poly_mat_get_coeff_mat(constant_mat, F_prime, k);

        nmod_mat_clear(A_k); /** clear A_k before compute a new A_k because
                              * mbasis1_for_mbasis init A_k
                              */

        rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_multiplication_blocks(res, A_k, perm, rank);
    }

    /** clear **/
    nmod_poly_mat_clear(F_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void mbasisII(nmod_poly_mat_t res, slong *res_shifts,
               const nmod_poly_mat_t F, ulong sigma, const slong *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;

    nmod_mat_poly0_t F_prime, res_list_repr;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    /** init **/
    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_mat_poly0_init_set(F_prime, F);
    nmod_mat_poly0_get_coef(constant_mat, F_prime, 0); //Compute F mod x

    rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank); //Compute P0

    for (ulong k = 1; k < sigma; k++)
    {
        nmod_mat_poly0_init_set(res_list_repr, res);

        nmod_mat_poly0_naive_mul_coef(constant_mat, res_list_repr, F_prime, k);

        nmod_mat_clear(A_k);
        rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_multiplication_blocks(res, A_k, perm, rank);

        nmod_mat_poly0_clear(res_list_repr);
    }

    /** clear **/
    nmod_mat_poly0_clear(F_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void mbasisIII(nmod_poly_mat_t res, slong *res_shifts,
                const nmod_poly_mat_t F, ulong sigma, const slong *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_mat_poly0_t F_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    /** init **/
    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_mat_poly0_init_setII(F_prime, F, sigma);

    nmod_mat_poly0_get_coef(constant_mat, F_prime, 0);
    rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank);
    //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

    for (ulong k = 1; k < sigma; k++)
    {
        structured_list_multiplication_blocks(F_prime, A_k, perm, rank, k - 1, sigma);
        nmod_mat_poly0_get_coef(constant_mat, F_prime, k);

        nmod_mat_clear(A_k);
        rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_multiplication_blocks(res, A_k, perm, rank);
    }

    /** clear **/
    nmod_mat_poly0_clear(F_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void mbasisIV(nmod_poly_mat_t res, slong *res_shifts,
               const nmod_poly_mat_t F, ulong sigma, const slong *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_mat_poly0_t F_prime, res_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_mat_poly0_init_setII(F_prime, F, sigma);
    nmod_mat_poly0_get_coef(constant_mat, F_prime, 0);

    rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank);
    nmod_mat_poly0_init_setIII(res_prime, res, sigma);
    //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

    for (ulong k = 1; k < sigma; k++)
    {
        structured_list_multiplication_blocks(F_prime, A_k, perm, rank, k - 1, sigma);
        nmod_mat_poly0_get_coef(constant_mat, F_prime, k);

        nmod_mat_clear(A_k);
        rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_list_multiplication_blocks_full(res_prime, A_k, perm, rank);
    }

    nmod_poly_mat_set_from_mat_poly0(res, res_prime); // can be improved

    nmod_mat_poly0_clear(F_prime);
    nmod_mat_poly0_clear(res_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void mbasisV(nmod_poly_mat_t res, slong *res_shifts,
              const nmod_poly_mat_t F, ulong sigma, const slong *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_mat_poly0_t F_prime, res_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_mat_poly0_init_set(F_prime, F);
    nmod_mat_poly0_get_coef(constant_mat, F_prime, 0);

    rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank);
    nmod_mat_poly0_init_setIII(res_prime, res, sigma);

    for (ulong k = 1; k < sigma; k++)
    {
        nmod_mat_poly0_naive_mul_coef(constant_mat, res_prime, F_prime, k);

        nmod_mat_clear(A_k);
        rank = mbasis1_for_mbasis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_list_multiplication_blocks_full(res_prime, A_k, perm, rank);
    }

    nmod_poly_mat_set_from_mat_poly0(res, res_prime);

    nmod_mat_poly0_clear(F_prime);
    nmod_mat_poly0_clear(res_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}


void nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                          slong * shift,
                          const nmod_poly_mat_t pmat,
                          ulong order)
{
    nmod_mat_poly_t app, matp;
    nmod_mat_poly_init(matp, pmat->r, pmat->c, pmat->modulus);
    // TODO improve: set init
    nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, order);
    nmod_mat_poly_init(app, pmat->r, pmat->r, pmat->modulus);
    nmod_mat_poly_mbasis(app, shift, matp, order);
    // TODO improve: set init
    nmod_poly_mat_set_from_mat_poly(appbas, app);
    nmod_mat_poly_clear(matp);
    nmod_mat_poly_clear(app);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
