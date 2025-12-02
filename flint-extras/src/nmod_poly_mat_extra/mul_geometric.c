#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h" // for nmod_find_root

#include "nmod_poly_mat_multiply.h"

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input  (TODO make this consistent with existing functions)
 *  ASSUMPTION (not checked): existence of element of "large enough" order
 *  TODO -> fail flag when element not found
 *  FIXME -> version underscore with provided geometric progression?
 *           (if will be used with various degrees, may require recrafting interpolation for more flexible npoints)
 *  uses evaluation and interpolation at a geometric progression
 */
void nmod_poly_mat_mul_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong ellA, ellB, ellC;
    ulong i, j, ell, m, k, n;
    ulong p, w;
    nn_ptr val;
    nmod_t mod;
    nmod_geometric_progression_t F;

    m = A->r;
    k = A->c;
    n = B->c;
    p = A->modulus;

    if (m < 1 || n < 1 || k < 1)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_mul_geometric(T, A, B);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    // length = 0 iff matrix is zero
    ellA = nmod_poly_mat_max_length(A);
    ellB = nmod_poly_mat_max_length(B);

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    ellC = ellA + ellB - 1;  // length(C) = length(A) + length(B) - 1
    nmod_init(&mod, p);
    w = nmod_find_root(2*ellC, mod); /* TODO 2*ellC ok? */
    nmod_geometric_progression_init(F, w, ellC, mod);

    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);

#ifdef DIRTY_ALLOC_MATRIX
    /* TODO check if still has any interest after FLINT's new mat storage */
    // we alloc the memory for all matrices at once
    nn_ptr tmp = flint_malloc((m*k + k*n + m*n) * ellC * sizeof(ulong));
    nn_ptr ptr = tmp;
    for (i = 0; i < ellC; i++)
    {
        mod_A[i]->entries = ptr + i*m*k;
        mod_A[i]->stride = k;
        mod_A[i]->r = m;
        mod_A[i]->c = k;
        mod_A[i]->mod.n = mod.n;
        mod_A[i]->mod.norm = mod.norm;
        mod_A[i]->mod.ninv = mod.ninv;
    }
    ptr += ellC*m*k;

    for (i = 0; i < ellC; i++)
    {
        mod_B[i]->entries = ptr + i*k*n;
        mod_B[i]->stride = n;
        mod_B[i]->r = k;
        mod_B[i]->c = n;
        mod_B[i]->mod.n = mod.n;
        mod_B[i]->mod.norm = mod.norm;
        mod_B[i]->mod.ninv = mod.ninv;
    }
    ptr += ellC*k*n;

    for (i = 0; i < ellC; i++)
    {
        mod_C[i]->entries = ptr + i*m*n;
        mod_C[i]->stride = n;
        mod_C[i]->r = m;
        mod_C[i]->c = n;
        mod_C[i]->mod.n = mod.n;
        mod_C[i]->mod.norm = mod.norm;
        mod_C[i]->mod.ninv = mod.ninv;
    }
#else
    for (i = 0; i < ellC; i++)
    {
        nmod_mat_init(mod_A[i], m, k, p);
        nmod_mat_init(mod_B[i], k, n, p);
        nmod_mat_init(mod_C[i], m, n, p);
    }
#endif

    nmod_poly_struct * pol;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            pol = nmod_poly_mat_entry(A, i, j);
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, pol->coeffs, pol->length, F, ellC, mod);
            for (ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_A[ell], i, j) = val[ell];
        }
    }

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < n; j++)
        {
            pol = nmod_poly_mat_entry(B, i, j);
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, pol->coeffs, pol->length, F, ellC, mod);
            for (ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = val[ell];
        }
    }

    for (ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (ell = 0; ell < ellC; ell++)
                val[ell] = nmod_mat_entry(mod_C[ell], i, j);
            nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_mat_entry(C, i, j), val, F, ellC);
        }
    }

#ifdef DIRTY_ALLOC_MATRIX
    flint_free(tmp);
#else
    for (i = 0; i < ellC; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }
#endif

    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
    _nmod_vec_clear(val);
    nmod_geometric_progression_clear(F);
}
