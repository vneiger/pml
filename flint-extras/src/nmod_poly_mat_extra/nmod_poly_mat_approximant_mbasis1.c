#include <stdlib.h>
#include "nmod_mat_extra.h" // for left_nullspace, permutation
#include "nmod_poly_mat_extra.h" // for apply_perm
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_approximant.h"

/* type for sorting permutation */
typedef struct
{
    slong value;
    slong ord;
} slong_pair;

/* comparison function for quicksort */
static inline int compare(const void * a, const void * b)
{
    slong_pair aa = * (const slong_pair *) a;
    slong_pair bb = * (const slong_pair *) b;
    return aa.value - bb.value;
}

/** Creates a permutation from the sorting of a shift
 * After running this, perm is the list of integers which sorts
 * shift increasingly, i.e.
 * shift[perm[0]] <= shift[perm[1]] < ... < shift[perm[n-1]]
 *
 * \param perm permutation (list of integers)
 * \param shift list of integer to be sorted increasingly
 * \param n length of both perm and shift
 *
 * TODO make stable
 */
static inline void sort_and_create_perm(slong * perm, const slong * shift, slong n)
{
    slong_pair temp[n];
    for (slong i = 0; i < n; i++)
    {
        temp[i].value = shift[i];
        temp[i].ord = i;
    }

    qsort(temp, n, sizeof(slong_pair), compare);

    for (slong i = 0; i < n; i++)
        perm[i] = temp[i].ord;
}


void mbasis1(nmod_poly_mat_t appbas,
             slong *res_shifts,
             const nmod_mat_t mat,
             const slong *shifts)
{
    slong rdim = mat->r;
    slong i, j, alloc, temp[rdim],
    rank_mat, rank_kernel,
    *P, *P_inv, *comp, *comp_inv, *perm;
    nmod_mat_t K, mat_cp;
    nmod_poly_t One, constant;
    ulong prime = mat->mod.n;

    perm = _perm_init(rdim);
    sort_and_create_perm(perm, shifts, rdim);

    /* LU operation on pi*mat and extraction  L = [ [Lr, 0], [G, I] ] */

    nmod_mat_init_set(mat_cp, mat);

    nmod_mat_permute_rows(mat_cp, perm, NULL);

    P = _perm_init(rdim);

    rank_kernel = nmod_mat_left_nullspace_compact(K, P, mat_cp);
    rank_mat = rdim - rank_kernel;

    P_inv = _perm_init(rdim);
    for (i = 0; i < rank_kernel; i++)
        P_inv[rank_mat + i] = P[i];
    for (i = 0; i < rank_mat; i++)
        P_inv[i] = P[rank_kernel + i];

    _perm_inv(P_inv, P_inv, rdim);
    _perm_inv(perm, perm, rdim);

    /* Computation of the block matrix  [ [xIr, 0], [K, I_{rdim - rank}] ] */
    nmod_poly_init(One, prime);
    nmod_poly_set_coeff_ui(One, 0, 1); // 1

    for (i = rank_mat; i < rdim; i++)
        nmod_poly_set(nmod_poly_mat_entry(appbas, i, i), One);

    nmod_poly_shift_left(One, One, 1); // x
    for (i = 0; i < rank_mat; i++)
        nmod_poly_set(nmod_poly_mat_entry(appbas, i, i), One);

    nmod_poly_init(constant, prime);
    for (i = 0; i < rank_kernel; i++)
        for (j = 0; j < rank_mat; j++)
        {
            alloc = (slong) nmod_mat_get_entry(K, i, j);
            if (alloc != 0)
            {
                nmod_poly_set_coeff_ui(constant, 0, alloc);
                nmod_poly_set(nmod_poly_mat_entry(appbas, i + rank_mat, j), constant);
            }
        }
    /* Multiply by the permutations */
    comp = _perm_init(rdim);
    comp_inv = _perm_init(rdim);

    _perm_compose(comp, P_inv, perm, rdim);

    _perm_inv(comp_inv, comp, rdim);

    nmod_poly_mat_permute_columns(appbas, comp, NULL);
    nmod_poly_mat_permute_rows(appbas, comp_inv, NULL);

    /* Compute the new shift */
    apply_perm_to_vector(temp, shifts, comp, rdim);

    for (i = 0; i < rank_mat; i++)
        temp[i] += 1;

    apply_perm_to_vector(res_shifts, temp, comp_inv, rdim);

    nmod_mat_clear(mat_cp);
    nmod_mat_clear(K);
    nmod_poly_clear(constant);
    nmod_poly_clear(One);

    _perm_clear(perm);
    _perm_clear(P);
    _perm_clear(P_inv);
    _perm_clear(comp);
    _perm_clear(comp_inv);
}

slong mbasis1_for_mbasis(nmod_mat_t appbas,
                         slong *res_shifts,
                         slong *res_perm,
                         const nmod_mat_t mat,
                         const slong *shifts)
{
    slong rdim = mat->r;
    slong  i, rank_mat, rank_kernel, temp[rdim];

    slong *P, *P_inv, *comp_inv, *perm = _perm_init(rdim);
    nmod_mat_t K, mat_cp;

    sort_and_create_perm(perm, shifts, rdim);

    /* compute left kernel of perm*mat with rank profile */
    nmod_mat_init_set(mat_cp, mat);

    nmod_mat_permute_rows(mat_cp, perm, NULL);

    P = _perm_init(rdim);

    rank_kernel = nmod_mat_left_nullspace_compact(K, P, mat_cp);
    rank_mat = rdim - rank_kernel;

    P_inv = _perm_init(rdim);
    for (i = 0; i < rank_kernel; i++)
        P_inv[rank_mat + i] = P[i];
    for (i = 0; i < rank_mat; i++)
        P_inv[i] = P[rank_kernel + i];

    _perm_inv(perm, perm, rdim);
    _perm_inv(perm, perm, rdim);

    nmod_mat_init_set(appbas, K);

    /* Permutations */
    _perm_compose(res_perm, P_inv, perm, rdim);

    /* Compute the new shift */
    apply_perm_to_vector(temp, shifts, res_perm, rdim);

    for (i = 0; i < rank_mat; i++)
        temp[i] += 1;

    comp_inv = _perm_init(rdim);

    _perm_inv(comp_inv, res_perm, rdim);
    apply_perm_to_vector(res_shifts, temp, comp_inv, rdim);

    nmod_mat_clear(mat_cp);
    nmod_mat_clear(K);

    _perm_clear(P);
    _perm_clear(perm);
    _perm_clear(P_inv);
    _perm_clear(comp_inv);

    return rank_mat;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
