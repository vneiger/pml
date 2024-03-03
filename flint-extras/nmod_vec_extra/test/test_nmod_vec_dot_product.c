#include <assert.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot_product(ulong len, ulong bits1, ulong bits2, ulong n, flint_rand_t state)
{
    mp_limb_t res1, res2;
    mp_ptr v1, v2, v1r, v2r;
    nmod_t mod, mod1, mod2;
    ulong i;

    nmod_init(&mod, n);
    if (bits1 < FLINT_BITS) nmod_init(&mod1, UWORD(1) << bits1);
    else nmod_init(&mod1, UWORD_MAX);
    if (bits2 < FLINT_BITS) nmod_init(&mod2, UWORD(1) << bits2);
    else nmod_init(&mod2, UWORD_MAX);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);
    v1r = _nmod_vec_init(len);
    v2r = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod1);
    _nmod_vec_rand(v2, state, len, mod2);

    for (i = 0; i < len; i++)
    {
        v1r[i] = v1[i] % n;
        v2r[i] = v2[i] % n;
    }

    res1 = _nmod_vec_dot(v1r, v2r, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res2 = nmod_vec_dot_product(v1, v2, len, bits1, bits2, mod);
    assert (res1 == res2);

    _nmod_vec_clear(v1r);
    _nmod_vec_clear(v2r);
    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_randinit(state);

    printf("test running, various bitlengths, len from 1 to ~1000 (no error message means success)...\n");
    for (slong len = 1; len < 1000; len += 1)
    {
        if (len % 20 == 0)
        {
            printf("%ld..", len);
            fflush(stdout);
        }
        for (slong repeat = 0; repeat < 10; repeat++)
        {
            check_nmod_vec_dot_product(len, 1, 3, (UWORD(1) << 3) + 1, state);
            check_nmod_vec_dot_product(len, 3, 3, (UWORD(1) << 3) + 1, state);
            check_nmod_vec_dot_product(len, 5, 10, (UWORD(1) << 10) + 1, state);
            check_nmod_vec_dot_product(len, 10, 10, (UWORD(1) << 10) + 1, state);
            check_nmod_vec_dot_product(len, 10, 20, (UWORD(1) << 20) + 1, state);
            check_nmod_vec_dot_product(len, 20, 20, (UWORD(1) << 20) + 1, state);
            check_nmod_vec_dot_product(len, 15, 29, (UWORD(1) << 29) + 1, state);
            check_nmod_vec_dot_product(len, 29, 29, (UWORD(1) << 29) + 1, state);
            check_nmod_vec_dot_product(len, 15, 30, (UWORD(1) << 30) + 1, state);
            check_nmod_vec_dot_product(len, 30, 30, (UWORD(1) << 30) + 1, state);
            check_nmod_vec_dot_product(len, 15, 31, (UWORD(1) << 31) + 1, state);
            check_nmod_vec_dot_product(len, 31, 31, (UWORD(1) << 31) + 1, state);
            check_nmod_vec_dot_product(len, 16, 32, (UWORD(1) << 32) + 1, state);
            check_nmod_vec_dot_product(len, 32, 32, (UWORD(1) << 32) + 1, state);
            check_nmod_vec_dot_product(len, 20, 40, (UWORD(1) << 40) + 1, state);
            check_nmod_vec_dot_product(len, 40, 40, (UWORD(1) << 40) + 1, state);
            check_nmod_vec_dot_product(len, 25, 50, (UWORD(1) << 50) + 1, state);
            check_nmod_vec_dot_product(len, 50, 50, (UWORD(1) << 50) + 1, state);
            check_nmod_vec_dot_product(len, 30, 60, (UWORD(1) << 60) + 1, state);
            check_nmod_vec_dot_product(len, 60, 60, (UWORD(1) << 60) + 1, state);
            check_nmod_vec_dot_product(len, 32, 64, UWORD_MAX, state);
            check_nmod_vec_dot_product(len, 64, 64, UWORD_MAX, state);
        }
    }
    printf("\n");

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
