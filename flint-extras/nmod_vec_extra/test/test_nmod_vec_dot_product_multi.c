#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

#include "nmod_vec_extra.h"
#include "nmod_mat_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot_product_multi(ulong len, ulong k, ulong bits1, ulong bits2, ulong n, flint_rand_t state)
{
    // moduli
    nmod_t mod, mod1, mod2;
    nmod_init(&mod, n);

    if (bits1 < FLINT_BITS) nmod_init(&mod1, (UWORD(1) << bits1) - 1);
    else nmod_init(&mod1, UWORD_MAX);
    if (bits2 < FLINT_BITS) nmod_init(&mod2, (UWORD(1) << bits2) - 1);
    else nmod_init(&mod2, UWORD_MAX);
    //printf("%lu, %lu, %lu\n", mod.n, mod1.n, mod2.n);

    // input vector
    nn_ptr u = _nmod_vec_init(len);
    nn_ptr ur = _nmod_vec_init(len);

    // input collection of vectors (= matrices)
    nn_ptr * v = flint_malloc(len * sizeof(nn_ptr));
    for (ulong i = 0; i < len; i++)
        v[i] = _nmod_vec_init(k);
    nmod_mat_t vr; nmod_mat_init(vr, len, k, n);

    // fill u / v at random, and reduce+copy into ur / vr
    _nmod_vec_rand(u, state, len, mod1);
    for (ulong i = 0; i < len; i++)
        ur[i] = u[i] % n;

    for (ulong i = 0; i < len; i++)
    {
        _nmod_vec_rand(v[i], state, k, mod2);
        for (ulong j = 0; j < k; j++)
            vr->rows[i][j] = v[i][j] % n;
    }

    // compute
    nn_ptr uv1 = _nmod_vec_init(k);
    nmod_vec_dot_product_multi(uv1, u, (nn_srcptr *) v, len, k, bits1, bits2, mod);

    nn_ptr uv2 = _nmod_vec_init(k);
    nmod_mat_nmod_vec_mul(uv2, ur, len, vr);

    assert (_nmod_vec_equal(uv1, uv2, k));

    // in case, also test mul_nmod_vec_newdot
    nmod_mat_t vrt;
    nmod_mat_init(vrt, vr->c, vr->r, vr->mod.n);
    nmod_mat_transpose(vrt, vr);
    nmod_mat_mul_nmod_vec_newdot(uv2, vrt, ur, len);
    assert (_nmod_vec_equal(uv1, uv2, k) && "mat_mul_nmod_vec_newdot");

    if (FLINT_BIT_COUNT(n) < 21)
    {
        _nmod_vec_dot_product_multi_1_v1_8(uv1, u, (nn_srcptr *) v, len, k, mod);
        assert (_nmod_vec_equal(uv1, uv2, k) && "v1_8");
        _nmod_vec_dot_product_multi_1_v8_8(uv1, u, (nn_srcptr *) v, len, k, mod);
        assert (_nmod_vec_equal(uv1, uv2, k) && "v8_8");
        _nmod_vec_dot_product_multi_1_v16_16(uv1, u, (nn_srcptr *) v, len, k, mod);
        assert (_nmod_vec_equal(uv1, uv2, k) && "v8_32");
        _nmod_vec_dot_product_multi_1_v16_16(uv1, u, (nn_srcptr *) v, len, k, mod);
        assert (_nmod_vec_equal(uv1, uv2, k) && "v16_16");
    }

    _nmod_vec_clear(u);
    for (ulong i = 0; i < len; i++)
        _nmod_vec_clear(v[i]);
    flint_free(v);
    _nmod_vec_clear(uv1);

    _nmod_vec_clear(ur);
    nmod_mat_clear(vr);
    _nmod_vec_clear(uv2);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("test running, various bitlengths, len and k from 1 to 300 (no error message means success)...\n");
    printf("len = ");
    for (slong len = 1; len < 300; len += 5)
    {
        if ((len-1) % 20 == 0)
        {
            printf("%ld..", len);
            fflush(stdout);
        }
        for (slong k = 1; k < 300; k += 5)
        {
            for (slong repeat = 0; repeat < 1; repeat++)
            {
                //check_nmod_vec_dot_product_multi(len, k, 1, 3, (UWORD(1) << 3) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 3, 3, (UWORD(1) << 3) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 5, 10, (UWORD(1) << 10) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 10, 10, (UWORD(1) << 10) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 10, 20, (UWORD(1) << 20) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 20, 20, (UWORD(1) << 20) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 15, 29, (UWORD(1) << 29) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 29, 29, (UWORD(1) << 29) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 15, 30, (UWORD(1) << 30) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 30, 30, (UWORD(1) << 30) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 15, 31, (UWORD(1) << 31) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 31, 31, (UWORD(1) << 31) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 16, 32, (UWORD(1) << 32) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 32, 32, (UWORD(1) << 32) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 20, 40, (UWORD(1) << 40) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 40, 40, (UWORD(1) << 40) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 25, 50, (UWORD(1) << 50) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 50, 50, (UWORD(1) << 50) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 30, 60, (UWORD(1) << 60) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 60, 60, (UWORD(1) << 60) + 1, state);
                //check_nmod_vec_dot_product_multi(len, k, 32, 64, UWORD_MAX, state);
                check_nmod_vec_dot_product_multi(len, k, 62, 62, UWORD_MAX, state);
                check_nmod_vec_dot_product_multi(len, k, 64, 64, UWORD_MAX, state);
            }
        }
    }
    printf("\n");

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
