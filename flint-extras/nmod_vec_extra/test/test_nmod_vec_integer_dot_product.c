#include <assert.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes an integer dot-product, length len, bitsize bit_len */
/* checks modulo a prime                                        */
/*--------------------------------------------------------------*/
void check_nmod_vec_integer_dot_product(slong len, mp_bitcnt_t bit_len, flint_rand_t state)
{
    mp_limb_t p, res1, res2;
    mp_ptr v1, v2, res;
    nmod_t mod;

    p = n_randprime(state, bit_len, 0);
    nmod_init(&mod, p);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_randtest(v1, state, len, mod);
    _nmod_vec_randtest(v2, state, len, mod);

    res = _nmod_vec_init(3);
    nmod_vec_integer_dot_product(res, v1, v2, len, bit_len, bit_len);
    NMOD_RED3(res1, res[2], res[1], res[0], mod);
    res2 = _nmod_vec_dot(v1, v2, len, mod, 3);

    assert (res1 == res2);

    _nmod_vec_clear(res);
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
    for (slong len = 1; len < 1000; len += 20)
    {
        printf("%ld..", len);
        fflush(stdout);
        for (slong repeat = 0; repeat < 100; repeat++)
        {
            check_nmod_vec_integer_dot_product(len, 3, state);
            check_nmod_vec_integer_dot_product(len, 10, state);
            check_nmod_vec_integer_dot_product(len, 20, state);
            check_nmod_vec_integer_dot_product(len, 29, state);
            check_nmod_vec_integer_dot_product(len, 30, state);
            check_nmod_vec_integer_dot_product(len, 31, state);
            check_nmod_vec_integer_dot_product(len, 32, state);
            check_nmod_vec_integer_dot_product(len, 40, state);
            check_nmod_vec_integer_dot_product(len, 50, state);
            check_nmod_vec_integer_dot_product(len, 60, state);
            check_nmod_vec_integer_dot_product(len, 64, state);
        }
    }
    printf("\n");

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
