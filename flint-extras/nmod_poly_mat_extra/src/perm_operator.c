#include "nmod_poly_mat_approximant.h"

/* Type only for the function Basis */
typedef struct
{
    int64_t value;
    slong ord;
} int_tuple;


/* Parameter for quicksort */
static int compare(const void *a, const void *b)
{
    int_tuple a_prime = * (const int_tuple *) a;
    int_tuple b_prime = * (const int_tuple *) b;
    return a_prime.value - b_prime.value;
}

/**
 * \brief Will only be used in Basis.
 *
 *  Creates a permutation from the sorting of a shift
 *
 * \param perm, the result his length must be n
 * \param vec, the shift we want to sort increasly
 * \param n, length of perm and vec
 */
void sort_and_create_perm(slong *perm, const int64_t *vec, slong n)
{
    int_tuple temp[n];
    for (slong i = 0; i < n; i++)
    {
        temp[i].value = vec[i];
        temp[i].ord = i;
    }

    qsort(temp, n, sizeof(int_tuple), compare);
    for (slong i = 0; i < n; i++)
        perm[temp[i].ord] = i;

}

/**
 * \brief A function fund in the doc of nmod_mat_t
 *
 */
static void nmod_poly_mat_swap_cols(nmod_poly_mat_t mat, slong * perm,
                                    slong r, slong s)
{
    if (r != s)
    {
        slong t;
        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        for (t = 0; t < mat->r; t++)
            nmod_poly_swap(mat->rows[t] + r, mat->rows[t] + s);
    }
}

static void nmod_poly_mat_swap_rows(nmod_poly_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        nmod_poly_struct *u;
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}


void apply_perm_rows_to_matrix(nmod_mat_t mat, const slong * perm, slong rdim)
{
    slong *perm_copy = _perm_init(rdim);
    _perm_set(perm_copy, perm, rdim);
    for (slong i = 0; i < rdim; i++)
        for (slong j = 0; j < rdim; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
        }

    _perm_clear(perm_copy);
}

void apply_perm_cols_to_poly_matrix(nmod_poly_mat_t mat, const slong * perm, slong cdim)
{
    slong *perm_copy = _perm_init(cdim);
    _perm_inv(perm_copy, perm, cdim);
    for (slong i = 0; i < cdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_poly_mat_swap_cols(mat, perm_copy, i, perm_copy[i]);
        }
    }
    _perm_clear(perm_copy);
}

void apply_perm_rows_to_poly_matrix(nmod_poly_mat_t mat, const slong *perm, slong rdim)
{
    slong *perm_copy = _perm_init(rdim);
    _perm_set(perm_copy, perm, rdim);

    for (slong i = 0; i < rdim; i++)
        for (slong j = 0; j < rdim; j++)
        {
            if (i == perm_copy[i])
                break;
            nmod_poly_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
        }
    _perm_clear(perm_copy);
}



void apply_perm_to_vector(int64_t *res, const int64_t *initial_vect,
                          const slong *perm, slong length)
{
    for (slong i = 0; i < length; i++)
        res[perm[i]] = initial_vect[i];
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
