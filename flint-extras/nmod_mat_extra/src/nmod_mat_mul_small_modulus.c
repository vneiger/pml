#include <flint/nmod_mat.h>
#include <flint/machine_vectors.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"

/** matrix multiplication using AVX2 instructions for moduli less than 2^30 */
void nmod_mat_mul_small_modulus(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    ulong i, j, m, k, n;
    nmod_mat_t BT;
    double p, pinv;
    vec2d p2, pinv2;
    mp_limb_t res[2];
    mp_limb_t power_two;
    
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

    p = A->mod.n;
    pinv = 1 / (double) p;
    p2 = vec2d_set_d2(p, p);
    pinv2 = vec2d_set_d2(pinv, pinv);
    
    nmod_mat_init(BT, n, k, B->mod.n);
    nmod_mat_transpose(BT, B);

    NMOD_RED(power_two, 1L << 45, A->mod);
    
    for (i = 0; i + 1 < m; i += 2)
        for (j = 0; j < n; j++)
        {
            _nmod_vec_dot2_small_modulus(res, A->rows[i], A->rows[i+1], BT->rows[j], k, power_two, p2, pinv2);
            C->rows[i][j] = res[0];
            C->rows[i+1][j] = res[1];
        }

    for (; i < m; i ++)
        for (j = 0; j < n; j++)
            C->rows[i][j] = _nmod_vec_dot_small_modulus(A->rows[i], BT->rows[j], k, power_two, p, pinv);

    nmod_mat_clear(BT);
}
