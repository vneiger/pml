#include <stdlib.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/machine_vectors.h>

#include "nmod_vec_extra.h"

/** matrix multiplication using AVX2 instructions for moduli less than 2^30 */
void nmod_mat_mul_small_modulus(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    ulong i, j, m, k, n;
    nmod_mat_t BT;
    m = A->r;
    k = A->c;
    n = B->c;
    
    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, m, n, A->mod.n);
        nmod_mat_mul_small_modulus(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return;
    }

    nmod_mat_init(BT, n, k, B->mod.n);
    nmod_mat_transpose(BT, B);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            C->rows[i][j] = _nmod_vec_dot_small_modulus(A->rows[i], BT->rows[j], k, A->mod);
    
    nmod_mat_clear(BT);
}
