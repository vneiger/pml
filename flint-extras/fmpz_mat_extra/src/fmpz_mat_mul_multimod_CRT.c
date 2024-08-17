#include <flint/nmod_mat.h>

#include "nmod_vec_extra.h"
#include "nmod_mat_extra.h"
#include "fmpz_extra.h"
#include "fmpz_mat_extra.h"

static void _fmpz_mat_mul_multimod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B, int sign, flint_bitcnt_t bits)
{
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong i, j, l, m, k, n, num_primes;
    flint_bitcnt_t primes_bits;
    fmpz_multimod_CRT_t CRT;
    nn_ptr primes, residues;
    fmpz_t half_top, C_nosign;
    
    m = A->r;
    k = A->c;
    n = B->c;

    if (m < 1 || n < 1 || k < 1)
    {
        fmpz_mat_zero(C);
        return;
    }

    bits += sign;

    primes_bits = 29;
    num_primes = (bits + primes_bits - 1) / primes_bits;
    primes = FLINT_ARRAY_ALLOC(num_primes, ulong);
    nmod_vec_primes(primes, num_primes, primes_bits + 1);

    mod_A = FLINT_ARRAY_ALLOC(num_primes, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(num_primes, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(num_primes, nmod_mat_t);
    
    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_init(mod_A[i], m, k, primes[i]);
        nmod_mat_init(mod_B[i], k, n, primes[i]);
        nmod_mat_init(mod_C[i], m, n, primes[i]);
    }
    
    fmpz_multimod_CRT_init(CRT, primes, num_primes);
    fmpz_init(half_top);

    if (sign == 1)
        fmpz_fdiv_q_2exp(half_top, CRT->product_primes, 1);

    residues = FLINT_ARRAY_ALLOC(num_primes, ulong);
    
    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
        {
            fmpz_multimod_CRT_reduce(residues, &A->rows[i][j], CRT);
            for (l = 0; l < num_primes; l++)
                mod_A[l]->rows[i][j] = residues[l];
        }

    
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
        {
            fmpz_multimod_CRT_reduce(residues, &B->rows[i][j], CRT);
            for (l = 0; l < num_primes; l++)
                mod_B[l]->rows[i][j] = residues[l];
        }


    for (i = 0; i < num_primes; i++)
        nmod_mat_mul_small_modulus(mod_C[i], mod_A[i], mod_B[i]);
    

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
            for (l = 0; l < num_primes; l++)
                residues[l] = mod_C[l]->rows[i][j];
            
            fmpz_multimod_CRT_CRT(&C->rows[i][j], residues, CRT);

            /* T > M/2 iff floor(M/2) < T */
            if (sign == 1 && fmpz_cmp(half_top, &C->rows[i][j]) < 0)
                fmpz_sub(&C->rows[i][j], &C->rows[i][j], CRT->product_primes);
        }

    flint_free(residues);
    fmpz_multimod_CRT_clear(CRT);

    for (i = 0; i < num_primes; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }

    fmpz_clear(C_nosign);
    fmpz_clear(half_top);
    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
    flint_free(primes);
}

/** ------------------------------------------------------------ */
/** a multimodular algorithm for matrix multiplication           */
/** based on flint's fmpz_mat_mul_multi_mod implementation       */
/** uses our fmpz_multimod_CRT nmod_mat_mul_small_modulus        */ 
/** ------------------------------------------------------------ */
void fmpz_mat_mul_multimod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong Abits, Bbits;
    int sign = 0;
    flint_bitcnt_t Cbits;

    Abits = fmpz_mat_max_bits(A);
    Bbits = fmpz_mat_max_bits(B);

    if (Abits < 0)
    {
        sign = 1;
        Abits = -Abits;
    }

    if (Bbits < 0)
    {
        sign = 1;
        Bbits = -Bbits;
    }

    Cbits = Abits + Bbits + FLINT_BIT_COUNT(A->c);

    _fmpz_mat_mul_multimod(C, A, B, sign, Cbits);
}

