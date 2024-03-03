#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#define DIRTY_ALLOC_MATRIX

#include "nmod_extra.h"
#include "nmod_poly_extra.h"
#include "nmod_poly_mat_multiply.h"

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
    mp_limb_t p, w;
    mp_ptr val, val2;
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
    nmod_geometric_progression_init_set(F, w, order, mod);


    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);
    val2 = _nmod_vec_init(ellC);
    nmod_poly_init2(tmp_poly, mod.n, dA+1);

    
#ifdef DIRTY_ALLOC_MATRIX
    // we alloc the memory for all matrices at once
    mp_ptr *tmp_rows = (mp_ptr *) malloc((m + k + m) * ellC * sizeof(mp_ptr));
    mp_ptr tmp = (mp_ptr) malloc((m*k + k*n + m*n) * ellC * sizeof(mp_limb_t));
    mp_ptr *bak_rows;
    mp_ptr bak;
    
    bak_rows = tmp_rows;
    j = 0;
    for (i = 0; i < m*ellC; i++)
    {
        tmp_rows[i] = tmp + j;
        j += k;
    }
    tmp_rows += m*ellC;
    
    for (i = 0; i < k*ellC; i++)
    {
        tmp_rows[i] = tmp + j;
        j += n;
    }
    tmp_rows += k*ellC;
    
    for (i = 0; i < m*ellC; i++)
    {
        tmp_rows[i] = tmp + j;
        j += n;
    }
    tmp_rows = bak_rows;

    bak_rows = tmp_rows;
    bak = tmp;
    for (i = 0; i < ellC; i++)
    {
        mod_A[i]->rows = tmp_rows + i*m;
        mod_A[i]->entries = tmp + i*m*k;
        mod_A[i]->r = m;
        mod_A[i]->c = k;
        mod_A[i]->mod.n = mod.n;
        mod_A[i]->mod.norm = mod.norm;
        mod_A[i]->mod.ninv = mod.ninv;
    }
    tmp_rows += ellC*m;
    tmp += ellC*m*k;
    
    for (i = 0; i < ellC; i++)
    {
        mod_B[i]->rows = tmp_rows + i*k;
        mod_B[i]->entries = tmp + i*k*n;
        mod_B[i]->r = k;
        mod_B[i]->c = n;
        mod_B[i]->mod.n = mod.n;
        mod_B[i]->mod.norm = mod.norm;
        mod_B[i]->mod.ninv = mod.ninv;
    }
    tmp_rows += ellC*k;
    tmp += ellC*k*n;
    
    for (i = 0; i < ellC; i++)
    {
        mod_C[i]->rows = tmp_rows + i*m;
        mod_C[i]->entries = tmp + i*m*n;
        mod_C[i]->r = m;
        mod_C[i]->c = n;
        mod_C[i]->mod.n = mod.n;
        mod_C[i]->mod.norm = mod.norm;
        mod_C[i]->mod.ninv = mod.ninv;
    }
    tmp_rows = bak_rows;
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
            ulong deg = nmod_poly_degree(&A->rows[i][j]);
            mp_ptr src = (A->rows[i] + j)->coeffs;
            mp_ptr dest = tmp_poly->coeffs;
            v = dA;
            for (u = 0; u <= deg; u++, v--)
                dest[v] = src[u];
            for (; v >= 0; v--)
                dest[v] = 0;
            tmp_poly->length = dA + 1;
            _nmod_poly_normalise(tmp_poly);
            
            nmod_geometric_progression_evaluate(val, tmp_poly, F);
            for (ell = 0; ell < ellC; ell++)
                mod_A[ell]->rows[i][j] = val[ell];
        }


    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
        {
            ulong deg = nmod_poly_degree(&B->rows[i][j]);
            mp_ptr src = (B->rows[i] + j)->coeffs;
            for (u = 0; u <= deg; u++)
                val2[u] = src[u];
            for (; u < ellC; u++)
                val2[u] = 0;
            nmod_geometric_progression_interpolate(tmp_poly, val2, F);
            deg = nmod_poly_degree(tmp_poly);
            src = tmp_poly->coeffs;
            for (ell = 0; ell <= deg; ell++)
                mod_B[ell]->rows[i][j] = src[ell];
            for (; ell < ellC; ell++)
                mod_B[ell]->rows[i][j] = 0;
        }

    for (ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);


    nmod_poly_fit_length(tmp_poly, ellC);
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
            mp_ptr dest = tmp_poly->coeffs;
            for (ell = 0; ell < ellC; ell++)
                dest[ell] = mod_C[ell]->rows[i][j];
            tmp_poly->length = ellC;
            _nmod_poly_normalise(tmp_poly);

            nmod_geometric_progression_evaluate(val2, tmp_poly, F);

            nmod_poly_realloc(C->rows[i] + j, dB + 1);
            (C->rows[i] + j)->length = dB + 1;
            dest = (C->rows[i] + j)->coeffs;
            for (u = 0; u <= dB; u++)
                dest[u] = val2[u];
            _nmod_poly_normalise(C->rows[i] + j);
        }
    
    
#ifdef DIRTY_ALLOC_MATRIX
    free(tmp_rows);
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
