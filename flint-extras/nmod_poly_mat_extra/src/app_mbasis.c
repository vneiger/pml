#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_matpoly.h"
#include "nmod_poly_mat_utils.h"

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
    apply_perm_rows_to_poly_matrix(res, perm, rdim);

    /** work on R1 (top) **/
    nmod_poly_mat_window_init(R1, res, 0, 0, rank, cdim);
    nmod_poly_mat_set(R1_cp, R1);

    nmod_poly_mat_shift(R1, 1); //matpol function

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
    apply_perm_rows_to_poly_matrix(res, inv_perm, rdim);

    /** clear **/
    nmod_poly_clear(P);

    nmod_poly_mat_clear(A_poly);
    nmod_poly_mat_clear(R1_cp);
    nmod_poly_mat_clear(R2_cp);

    nmod_poly_mat_window_clear(R1);
    nmod_poly_mat_window_clear(R2);

    _perm_clear(inv_perm);
}

void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
             const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts)

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
    coefficient_matrix(constant_mat, F_prime, 0); //Compute F mod x

    rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank); //Compute P0
    for (ulong k = 1; k < sigma; k++)
    {
        // doing the operation x^(-k) F_prime mod x
        structured_multiplication_blocks(F_prime, A_k, perm, rank);
        coefficient_matrix(constant_mat, F_prime, k);

        nmod_mat_clear(A_k); /** clear A_k before compute a new A_k because
                              * Basis for M basis init A_k
                              */

        rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_multiplication_blocks(res, A_k, perm, rank);
    }

    /** clear **/
    nmod_poly_mat_clear(F_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
               const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;

    nmod_list_poly_mat_t F_prime, res_list_repr;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    /** init **/
    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_list_poly_mat_init_set(F_prime, F);
    nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0); //Compute F mod x

    rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank); //Compute P0

    for (ulong k = 1; k < sigma; k++)
    {
        nmod_list_poly_mat_init_set(res_list_repr, res);

        nmod_list_poly_mat_naive_mul_coef(constant_mat, res_list_repr, F_prime, k);

        nmod_mat_clear(A_k);
        rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_multiplication_blocks(res, A_k, perm, rank);

        nmod_list_poly_mat_clear(res_list_repr);
    }

    /** clear **/
    nmod_list_poly_mat_clear(F_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
                const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_list_poly_mat_t F_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    /** init **/
    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_list_poly_mat_init_setII(F_prime, F, sigma);

    nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0);
    rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank);
    //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

    for (ulong k = 1; k < sigma; k++)
    {
        structured_list_multiplication_blocks(F_prime, A_k, perm, rank, k - 1, sigma);
        nmod_list_poly_mat_get_coef(constant_mat, F_prime, k);

        nmod_mat_clear(A_k);
        rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_multiplication_blocks(res, A_k, perm, rank);
    }

    /** clear **/
    nmod_list_poly_mat_clear(F_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void M_basisIV(nmod_poly_mat_t res, int64_t *res_shifts,
               const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_list_poly_mat_t F_prime, res_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_list_poly_mat_init_setII(F_prime, F, sigma);
    nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0);

    rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank);
    nmod_list_poly_mat_init_setIII(res_prime, res, sigma);
    //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

    for (ulong k = 1; k < sigma; k++)
    {
        structured_list_multiplication_blocks(F_prime, A_k, perm, rank, k - 1, sigma);
        nmod_list_poly_mat_get_coef(constant_mat, F_prime, k);

        nmod_mat_clear(A_k);
        rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_list_multiplication_blocks_full(res_prime, A_k, perm, rank);
    }

    nmod_list_poly_mat_to_poly_mat(res, res_prime); // can be improved

    nmod_list_poly_mat_clear(F_prime);
    nmod_list_poly_mat_clear(res_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

void M_basisV(nmod_poly_mat_t res, int64_t *res_shifts,
              const nmod_poly_mat_t F, ulong sigma, const int64_t *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_list_poly_mat_t F_prime, res_prime;
    nmod_mat_t A_k, constant_mat;
    slong rank, *perm;

    perm = _perm_init(rdim);
    nmod_mat_init(constant_mat, rdim, cdim, prime);

    nmod_list_poly_mat_init_set(F_prime, F);
    nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0);

    rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

    nmod_poly_mat_one(res);
    structured_multiplication_blocks(res, A_k, perm, rank);
    nmod_list_poly_mat_init_setIII(res_prime, res, sigma);

    for (ulong k = 1; k < sigma; k++)
    {
        nmod_list_poly_mat_naive_mul_coef(constant_mat, res_prime, F_prime, k);

        nmod_mat_clear(A_k);
        rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

        structured_list_multiplication_blocks_full(res_prime, A_k, perm, rank);
    }

    nmod_list_poly_mat_to_poly_mat(res, res_prime);

    nmod_list_poly_mat_clear(F_prime);
    nmod_list_poly_mat_clear(res_prime);
    nmod_mat_clear(constant_mat);
    nmod_mat_clear(A_k);
    _perm_clear(perm);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
