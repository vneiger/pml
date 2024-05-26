#include <flint/longlong.h> // for BIT_COUNT
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>  // for _nmod_vec_init

#include "nmod_vec_extra.h"  // for dot_product

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

    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // now let's compute
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
            C->rows[i][j] = nmod_vec_dot_product_v1(A->rows[i], BT->rows[j], A->c, A->mod);
            //C->rows[i][j] = nmod_vec_dot_product_v2(A->rows[i], BT->rows[j], A->c, A->mod, 1);
            //C->rows[i][j] = _nmod_vec_dot_product_1_avx2(A->rows[i], BT->rows[j], A->c, A->mod);
            //C->rows[i][j] = _nmod_vec_dot_product_1_avx512(A->rows[i], BT->rows[j], A->c, A->mod);

    nmod_mat_clear(BT);
}

void nmod_mat_mul_nmod_vec_newdot(nn_ptr v, const nmod_mat_t A, nn_srcptr u, ulong len)
{
    for (slong i = 0; i < A->r; i++)
        v[i] = nmod_vec_dot_product_v1(A->rows[i], u, len, A->mod);
}

