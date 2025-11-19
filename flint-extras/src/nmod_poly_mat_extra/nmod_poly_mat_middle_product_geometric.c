#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_mat_multiply.h"

#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 4)

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  ASSUME: existence of primitive root
 *  uses geometric evaluation and interpolation
 */
void nmod_poly_mat_middle_product_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                            const ulong dA, const ulong dB)
{
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong ellA, ellB, ellC, order;
    ulong i, j, ell, m, k, n, u;
    long v;
    ulong p, w;
    nn_ptr val, val2;
    nmod_t mod;
    nmod_geometric_progression_t F;
    nmod_poly_t tmp_poly;

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
        nmod_poly_mat_middle_product_geometric(T, A, B, dA, dB);
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

    nmod_init(&mod, p);

    ellC = ellA + ellB - 1;  // length(C) = length(A) + length(B) - 1
    order = ellC;
    nmod_init(&mod, p);
    w = nmod_find_root(order, mod);
    nmod_geometric_progression_init(F, w, order, mod);

    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);
    val2 = _nmod_vec_init(ellC);
    nmod_poly_init2(tmp_poly, mod.n, dA+1);

#ifdef DIRTY_ALLOC_MATRIX
    // we alloc the memory for all matrices at once
    nn_ptr tmp = (nn_ptr) malloc((m*k + k*n + m*n) * ellC * sizeof(ulong));
    nn_ptr bak;

    bak = tmp;
    for (i = 0; i < ellC; i++)
    {
        mod_A[i]->entries = tmp + i*m*k;
        mod_A[i]->r = m;
        mod_A[i]->c = k;
        mod_A[i]->mod.n = mod.n;
        mod_A[i]->mod.norm = mod.norm;
        mod_A[i]->mod.ninv = mod.ninv;
    }
    tmp += ellC*m*k;

    for (i = 0; i < ellC; i++)
    {
        mod_B[i]->entries = tmp + i*k*n;
        mod_B[i]->r = k;
        mod_B[i]->c = n;
        mod_B[i]->mod.n = mod.n;
        mod_B[i]->mod.norm = mod.norm;
        mod_B[i]->mod.ninv = mod.ninv;
    }
    tmp += ellC*k*n;

    for (i = 0; i < ellC; i++)
    {
        mod_C[i]->entries = tmp + i*m*n;
        mod_C[i]->r = m;
        mod_C[i]->c = n;
        mod_C[i]->mod.n = mod.n;
        mod_C[i]->mod.norm = mod.norm;
        mod_C[i]->mod.ninv = mod.ninv;
    }
    tmp = bak;
#else
    for (i = 0; i < ellC; i++)
    {
        nmod_mat_init(mod_A[i], m, k, p);
        nmod_mat_init(mod_B[i], k, n, p);
        nmod_mat_init(mod_C[i], m, n, p);
    }
#endif

    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
        {
            ulong deg = nmod_poly_degree(nmod_poly_mat_entry(A, i, j));
            nn_ptr src = nmod_poly_mat_entry(A, i, j)->coeffs;
            nn_ptr dest = tmp_poly->coeffs;
            v = dA;
            for (u = 0; u <= deg; u++, v--)
                dest[v] = src[u];
            for (; v >= 0; v--)
                dest[v] = 0;
            tmp_poly->length = dA + 1;
            _nmod_poly_normalise(tmp_poly);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, tmp_poly->coeffs, tmp_poly->length, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_A[ell], i, j) = val[ell];
        }


    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
        {
            ulong deg = nmod_poly_degree(nmod_poly_mat_entry(B, i, j));
            nn_ptr src = nmod_poly_mat_entry(B, i, j)->coeffs;
            for (u = 0; u <= deg; u++)
                val2[u] = src[u];
            for (; u < ellC; u++)
                val2[u] = 0;
            nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(tmp_poly, val2, F, ellC);
            deg = nmod_poly_degree(tmp_poly);
            src = tmp_poly->coeffs;
            for (ell = 0; ell <= deg; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = src[ell];
            for (; ell < ellC; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = 0;
        }

    for (ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);


    nmod_poly_fit_length(tmp_poly, ellC);
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
            nn_ptr dest = tmp_poly->coeffs;
            for (ell = 0; ell < ellC; ell++)
                dest[ell] = nmod_mat_entry(mod_C[ell], i, j);
            tmp_poly->length = ellC;
            _nmod_poly_normalise(tmp_poly);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val2, tmp_poly->coeffs, tmp_poly->length, F, ellC);

            nmod_poly_realloc(nmod_poly_mat_entry(C, i, j), dB + 1);
            nmod_poly_mat_entry(C, i, j)->length = dB + 1;
            dest = nmod_poly_mat_entry(A, i, j)->coeffs;
            for (u = 0; u <= dB; u++)
                dest[u] = val2[u];
            _nmod_poly_normalise(nmod_poly_mat_entry(A, i, j));
        }


#ifdef DIRTY_ALLOC_MATRIX
    free(tmp);
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
    nmod_poly_clear(tmp_poly);
    _nmod_vec_clear(val2);
    _nmod_vec_clear(val);
    nmod_geometric_progression_clear(F);
}

#endif  /* (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 4) */
