#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/fft_small.h>

#include "nmod_poly_extra.h"

#ifdef TIME_TFT
#include <time.h>
#endif

#include "nmod_poly_mat_multiply.h"

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: 2^(ceiling(log_2(lenA+lenB-1))) divides p-1 (assumption not checked)
 *  uses tft multiplication
 */
void nmod_poly_mat_mul_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong ellA, ellB, ellC, order;
    ulong i, j, ell, m, k, n;
    ulong p, w;
    nn_ptr val;
    nmod_t mod;
    sd_fft_ctx_t Q;
    sd_fft_lctx_t QL;
    nmod_sd_fft_t F;
    
#ifdef TIME_TFT
    double t = 0.0;
    clock_t tt;
    tt = clock();
#endif    

    // straightforward algorithm: tft-evaluate A and B, multiply the values, interpolate C
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
    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("\ninit: %lf\n", t);
    tt = clock();
#endif

    sd_fft_ctx_init_prime(Q, p);

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("ctx init: %lf\n", t);
    tt = clock();
#endif

    w = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> order, mod); // w has order 2^order
    nmod_sd_fft_init_set(F, w, order, mod);
    sd_fft_lctx_init(QL, Q, order);
    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);
    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("fft init: %lf\n", t);
    tt = clock();
#endif

#ifdef DIRTY_ALLOC_MATRIX
    // we alloc the memory for all matrices at once
    nn_ptr *tmp_rows = (nn_ptr *) malloc((m + k + m) * ellC * sizeof(nn_ptr));
    nn_ptr tmp = (nn_ptr) malloc((m*k + k*n + m*n) * ellC * sizeof(ulong));
    nn_ptr *bak_rows;
    nn_ptr bak;
    
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

    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("mem init: %lf\n", t);
    tt = clock();
#endif
    
    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
        {
            nmod_sd_tft_evaluate(val, &A->rows[i][j], QL, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                mod_A[ell]->rows[i][j] = val[ell];
        }

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("eval A: %lf\n", t);
    tt = clock();
#endif
    
    for (i = 0; i < k; i++)
        for (j = 0; j < m; j++)
        {
            nmod_sd_tft_evaluate(val, &B->rows[i][j], QL, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                mod_B[ell]->rows[i][j] = val[ell];
        }

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("eval B: %lf\n", t);
    tt = clock();
#endif
    
    for (ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);

    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("mul: %lf\n", t);
    tt = clock();
#endif
    
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
        {
            for (ell = 0; ell < ellC; ell++)
                val[ell] = mod_C[ell]->rows[i][j];
            nmod_sd_tft_interpolate(&C->rows[i][j], val, QL, F, ellC);
        }
    

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("interpolate: %lf\n", t);
    tt = clock();
#endif

    
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
    _nmod_vec_clear(val);
    nmod_sd_fft_clear(F);
    sd_fft_lctx_clear(QL, Q);
    sd_fft_ctx_clear(Q);

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("clean: %lf\n", t);
#endif

}
