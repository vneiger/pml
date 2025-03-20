#include "nmod_vec_extra.h"
#include <flint/longlong.h> // for BIT_COUNT
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>  // for _nmod_vec_init

#ifdef __AVX2__
#define HAS_AVX2
#endif

#ifdef HAS_AVX2  // GV 

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

    const dot_params_t params = _nmod_vec_dot_params(A->c, A->mod);

    // now let's compute
    if (params.method == _DOT2_HALF)
    {
        for (slong i = 0; i < A->r; i++)
            for (slong j = 0; j < BT->r; j++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
                C->rows[i][j] = _nmod_vec_dot(A->rows[i], BT->rows[j], A->c, A->mod, params);
#else
                nmod_mat_entry(C, i, j) = _nmod_vec_dot2_half_avx(nmod_mat_entry_ptr(A, i, 0),
                                                nmod_mat_entry_ptr(BT, j, 0),
                                                A->c,
                                                A->mod);
#endif
    }
    else
    {
        for (slong i = 0; i < A->r; i++)
            for (slong j = 0; j < BT->r; j++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
                C->rows[i][j] = _nmod_vec_dot(A->rows[i], BT->rows[j], A->c, A->mod, params);
#else
                nmod_mat_entry(C, i, j) = _nmod_vec_dot(nmod_mat_entry_ptr(A, i, 0),
                                                nmod_mat_entry_ptr(BT, j, 0),
                                                A->c,
                                                A->mod,
                                                params);
#endif
    }

    nmod_mat_clear(BT);
}

void nmod_mat_mul_nmod_vec_newdot(nn_ptr v, const nmod_mat_t A, nn_srcptr u, ulong len)
{
    const dot_params_t params = _nmod_vec_dot_params(A->c, A->mod);
    if (params.method == _DOT2_HALF)
    {
        for (slong i = 0; i < A->r; i++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
            v[i] = _nmod_vec_dot(A->rows[i], u, len, A->mod, params);
#else
            v[i] = _nmod_vec_dot2_half_avx(nmod_mat_entry_ptr(A, i, 0), u, len, A->mod);
#endif
    }
    else
    {
        for (slong i = 0; i < A->r; i++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
            v[i] = _nmod_vec_dot(A->rows[i], u, len, A->mod, params);
#else
            v[i] = _nmod_vec_dot(nmod_mat_entry_ptr(A, i, 0), u, len, A->mod, params);
#endif
    }
}

#endif 
