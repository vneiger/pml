#include <flint/nmod_mat.h>

#include "nmod_poly_mat_multiply.h"

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  assume 2^(ceiling(log_2(lenA+lenB-1))) divides p-1
 *  uses tft multiplication
 */
void nmod_poly_mat_mul_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong ellA, ellB, ellC, order;
    ulong i, j, ell, m, k, n;
    mp_limb_t p, w;
    mp_ptr val;
    nmod_t mod;
    sd_fft_ctx_t Q;
    sd_fft_lctx_t QL;
    nmod_sd_fft_t F;


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
        nmod_poly_mat_mul_tft(T, A, B);
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
    order = n_clog(ellC, 2); // ceiling(log_2(ellC))

    nmod_init(&mod, p);

    sd_fft_ctx_init_prime(Q, p);
    w = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> order, mod); // w has order 2^order

    nmod_sd_fft_init_set(F, w, order, mod);
    sd_fft_lctx_init(QL, Q, order);

    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);

    for (i = 0; i < ellC; i++)
    {
        nmod_mat_init(mod_A[i], m, k, p);
        nmod_mat_init(mod_B[i], k, n, p);
        nmod_mat_init(mod_C[i], m, n, p);
    }


    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
        {
            nmod_sd_tft_evaluate(val, &A->rows[i][j], QL, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                mod_A[ell]->rows[i][j] = val[ell];
        }

    for (i = 0; i < k; i++)
        for (j = 0; j < m; j++)
        {
            nmod_sd_tft_evaluate(val, &B->rows[i][j], QL, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                mod_B[ell]->rows[i][j] = val[ell];
        }


    for (ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);

    
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
            for (ell = 0; ell < ellC; ell++)
                val[ell] = mod_C[ell]->rows[i][j];
            nmod_sd_tft_interpolate(&C->rows[i][j], val, QL, F, ellC);
        }
    

    for (i = 0; i < ellC; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }
    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
    _nmod_vec_clear(val);
    nmod_sd_fft_clear(F);
    sd_fft_lctx_clear(QL, Q);
    sd_fft_ctx_clear(Q);
}
