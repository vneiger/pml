#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_vec_extra.h"

void nmod_mat_mul_newdot(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    // C aliases A or B: use a temporary
    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, A->r, B->c, A->mod.n);
        nmod_mat_mul_newdot(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return;
    }

    // number of bits of modulus
    const ulong nbits = FLINT_BIT_COUNT(A->mod.n);

    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // let's go
    //slong i = 0;
    //for (; i < A->r - 1; i += 2)
    //    for (slong j = 0; j < BT->r; j++)
    //    {
    //        _nmod_vec_dot2_small_modulus(res, A->rows[i], A->rows[i+1], BT->rows[j], A->c, power_two, p2, pinv2);
    //        C->rows[i][j] = res[0];
    //        C->rows[i+1][j] = res[1];
    //    }

    //for (; i < A->r; i++)
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
            C->rows[i][j] = nmod_vec_dot_product(A->rows[i], BT->rows[j], A->c, nbits, nbits, A->mod);

    nmod_mat_clear(BT);
}
