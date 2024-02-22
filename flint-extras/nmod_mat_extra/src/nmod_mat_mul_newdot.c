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

    // number of bits of modulus
    const ulong nbits = FLINT_BIT_COUNT(A->mod.n);

    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // now let's compute
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
            C->rows[i][j] = nmod_vec_dot_product(A->rows[i], BT->rows[j], A->c, nbits, nbits, A->mod);

    nmod_mat_clear(BT);
}

void nmod_mat_mul_nmod_vec_newdot(mp_ptr v, const nmod_mat_t A, mp_srcptr u, ulong len)
{
    // number of bits of modulus
    const ulong nbits = FLINT_BIT_COUNT(A->mod.n);

    for (slong i = 0; i < A->r; i++)
        v[i] = nmod_vec_dot_product(A->rows[i], u, len, nbits, nbits, A->mod);
}

