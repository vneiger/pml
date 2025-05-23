#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_mat_multiply.h"

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  ASSUME: num columns of A < 2^30 and min(deg A, deg B) < 2^30 (assumption not checked)
 *  uses tft multiplication modulo 50 bits fft primes
 */
void nmod_poly_mat_mul_3_primes(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    ulong num_primes, num_bits, i, j, ell, m, k, n, len_A, len_B, len_C;
    ulong p, primes[4];
    nmod_multimod_CRT_t CRT;
    nn_ptr residues[4];
    nmod_poly_mat_t *mod_A, *mod_B, *mod_C;

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
        nmod_poly_mat_mul_3_primes(T, A, B);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    primes[0] = PRIME0;
    primes[1] = PRIME1;
    primes[2] = PRIME2;
    primes[3] = PRIME3;

    len_A = 0;
    len_B = 0;
   
    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
            len_A = FLINT_MAX(len_A, A->rows[i][j].length);
    
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            len_B = FLINT_MAX(len_B, B->rows[i][j].length);
    
    // slight overestimate
    num_bits = 2 * FLINT_BIT_COUNT(p) + FLINT_BIT_COUNT(k) + FLINT_BIT_COUNT(FLINT_MIN(len_A, len_B));
    // all entries of C fit in size len_C
    len_C = len_A + len_B - 1;

    // our assumption implies that num_bits < 4*49
    num_primes = 1;
    if (num_bits > 49)
        num_primes = 2;
    if (num_bits > 2*49)
        num_primes = 3;
    if (num_bits > 3*49)
        num_primes = 4;

    nmod_multimod_CRT_init(CRT, p, num_primes);

    mod_A = FLINT_ARRAY_ALLOC(num_primes, nmod_poly_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(num_primes, nmod_poly_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(num_primes, nmod_poly_mat_t);
    
    for (i = 0; i < num_primes; i++)
    {
        nmod_poly_mat_init(mod_A[i], m, k, primes[i]);
        nmod_poly_mat_init(mod_B[i], k, n, primes[i]);
        nmod_poly_mat_init(mod_C[i], m, n, primes[i]);
    }

    // computes A modulo all primes
    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
        {
            slong len;
            len = A->rows[i][j].length;
            for (ell = 0; ell < num_primes; ell++)
            {
                nmod_poly_realloc(&mod_A[ell]->rows[i][j], len);
                mod_A[ell]->rows[i][j].length = len;
                residues[ell] = mod_A[ell]->rows[i][j].coeffs;
            }
            nmod_multimod_CRT_reduce(residues, A->rows[i][j].coeffs, len, CRT);
            for (ell = 0; ell < num_primes; ell++)
                _nmod_poly_normalise(&mod_A[ell]->rows[i][j]);
        }


    // computes B modulo all primes
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
        {
            slong len;
            len = B->rows[i][j].length;
            for (ell = 0; ell < num_primes; ell++)
            {
                nmod_poly_realloc(&mod_B[ell]->rows[i][j], len);
                mod_B[ell]->rows[i][j].length = len;
                residues[ell] = mod_B[ell]->rows[i][j].coeffs;
            }
            nmod_multimod_CRT_reduce(residues, B->rows[i][j].coeffs, len, CRT);
            for (ell = 0; ell < num_primes; ell++)
                _nmod_poly_normalise(&mod_B[ell]->rows[i][j]);
        }


    // TFT multiplications
    for (ell = 0; ell < num_primes; ell++)
    {
        nmod_poly_mat_mul_tft(mod_C[ell], mod_A[ell], mod_B[ell]);
        residues[ell] = _nmod_vec_init(len_C);
    }


    // CRT to find C
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
            slong len;
            len = 0;
            for (ell = 0; ell < num_primes; ell++)
                len = FLINT_MAX(len, mod_C[ell]->rows[i][j].length);
            for (ell = 0; ell < num_primes; ell++)
            {
                ulong t, len_ell;
                nn_ptr res_ell, C_ell_ij;

                res_ell = residues[ell];
                C_ell_ij = mod_C[ell]->rows[i][j].coeffs;
                len_ell = mod_C[ell]->rows[i][j].length;
                for (t = 0; t < len_ell; t++)
                    res_ell[t] = C_ell_ij[t];
                for (; t < (ulong) len; t++)
                    res_ell[t] = 0;
            }
            nmod_poly_realloc(&C->rows[i][j], len);
            C->rows[i][j].length = len;
            nmod_multimod_CRT_CRT(C->rows[i][j].coeffs, residues, len, CRT);
            _nmod_poly_normalise(&C->rows[i][j]);
        }


    nmod_multimod_CRT_clear(CRT);
    for (i = 0; i < num_primes; i++)
    {
        nmod_poly_mat_clear(mod_A[i]);
        nmod_poly_mat_clear(mod_B[i]);
        nmod_poly_mat_clear(mod_C[i]);
        _nmod_vec_clear(residues[i]);
    }
    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
}


