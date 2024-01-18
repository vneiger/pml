#include <flint/nmod_mat.h>

#define DIRTY_ALLOC_MATRIX

#ifdef TIME_TFT
#include "time.h"
#endif

#include "nmod_poly_mat_multiply.h"

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1), assuming deg(A) <= dA and deg(B) <= dA + dB
 *  output can alias input
 *  assume 2^(ceiling(log_2(dA + dB + 1))) divides p-1
 *  uses tft multiplication
 */
void nmod_poly_mat_middle_product_tft(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                      const ulong dA, const ulong dB)
{
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong ellA, ellB, ellC, order;
    ulong i, j, ell, m, k, n, u;
    long v;
    mp_limb_t p, w;
    mp_ptr val, val2;
    nmod_t mod;
    sd_fft_ctx_t Q, Qt;
    sd_fft_lctx_t QL, QtL;
    nmod_sd_fft_t F;
    nmod_poly_t revA;
    
#ifdef TIME_TFT
    double t = 0.0;
    clock_t tt;
    tt = clock();
#endif    
    
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
        nmod_poly_mat_middle_product_tft(T, A, B, dA, dB);
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

    ellC = dA + dB + 1;  
    order = n_clog(ellC, 2); // ceiling(log_2(ellC))

    nmod_init(&mod, p);
    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("\ninit: %lf\n", t);
    tt = clock();
#endif

    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_init_inverse(Qt, Q);

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("ctx init: %lf\n", t);
    tt = clock();
#endif

    w = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> order, mod); // w has order 2^order
    nmod_sd_fft_init_set(F, w, order, mod);
    sd_fft_lctx_init(QL, Q, order);
    sd_fft_lctx_init(QtL, Qt, order);
    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);
    val2 = _nmod_vec_init(ellC);
    nmod_poly_init2(revA, mod.n, dA+1);

    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("fft init: %lf\n", t);
    tt = clock();
#endif

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

    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("mem init: %lf\n", t);
    tt = clock();
#endif
    
    
    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
        {
            ulong deg = nmod_poly_degree(&A->rows[i][j]);
            mp_ptr src = (A->rows[i] + j)->coeffs;
            mp_ptr dest = revA->coeffs;
            v = dA;
            for (u = 0; u <= deg; u++, v--)
                dest[v] = src[u];
            for (; v >= 0; v--)
                dest[v] = 0;
            revA->length = dA + 1;
            _nmod_poly_normalise(revA);
            
            nmod_sd_tft_evaluate(val, revA, QL, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                mod_A[ell]->rows[i][j] = val[ell];
        }


#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("eval A: %lf\n", t);
    tt = clock();
#endif
    
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
        {
            ulong deg = nmod_poly_degree(&B->rows[i][j]);
            mp_ptr src = (B->rows[i] + j)->coeffs;
            for (u = 0; u <= deg; u++)
                val2[u] = src[u];
            for (; u < ellC; u++)
                val2[u] = 0;
            nmod_sd_tft_interpolate_t(val, val2, QtL, F, ellC);
            for (ell = 0; ell < ellC; ell++)
                mod_B[ell]->rows[i][j] = val[ell];
        }

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("interpolate_t B: %lf\n", t);
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
            
            nmod_sd_tft_evaluate_t(val2, val, QtL, F, ellC);

            nmod_poly_realloc(C->rows[i] + j, dB + 1);
            (C->rows[i] + j)->length = dB + 1;
            mp_ptr dest = (C->rows[i] + j)->coeffs;
            for (u = 0; u <= dB; u++)
                dest[u] = val2[u];
            _nmod_poly_normalise(C->rows[i] + j);
        }
    
#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("evaluate_t: %lf\n", t);
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
    nmod_poly_clear(revA);
    _nmod_vec_clear(val2);
    _nmod_vec_clear(val);
    nmod_sd_fft_clear(F);
    sd_fft_lctx_clear(QtL, Qt);
    sd_fft_lctx_clear(QL, Q);
    sd_fft_ctx_clear(Qt);
    sd_fft_ctx_clear(Q);

#ifdef TIME_TFT
    t = (double)(clock()-tt) / CLOCKS_PER_SEC;
    printf("clean: %lf\n", t);
#endif

}
