#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>

#include "nmod_vec_extra.h"


void _nmod_vec_randtest_not_zero(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    slong i, sparseness;

    if (n_randint(state, 2))
    {
        for (i = 0; i < len; i++)
            vec[i] = 1 + (n_randtest(state) % (mod.n - 1));
    }
    else
    {
        sparseness = 1 + n_randint(state, FLINT_MAX(2, len));

        for (i = 0; i < len; i++)
        {
            if (n_randint(state, sparseness))
                vec[i] = 1;
            else
                vec[i] = 1 + (n_randtest(state) % (mod.n - 1));
        }
    }
}

void _nmod_vec_inv_naive(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
        res[k] = nmod_inv(vec[k], mod);
}


/*--------------------------------------------------------------*/
/* computes an integer dot-product, length len, bitsize bit_len */
/* checks modulo a prime                                        */
/*--------------------------------------------------------------*/
void check_nmod_vec_inv(slong len, flint_bitcnt_t bit_len, flint_rand_t state)
{
    ulong p;
    nn_ptr vec, res1, res2, res3;
    nmod_t mod;

    p = n_randprime(state, bit_len, 1);
    nmod_init(&mod, p);

    vec = _nmod_vec_init(len);
    _nmod_vec_randtest_not_zero(vec, state, len, mod);

    res1 = _nmod_vec_init(len);
    _nmod_vec_inv(res1, vec, len, mod);

    res2 = _nmod_vec_init(len);
    _nmod_vec_inv_naive(res2, vec, len, mod);

    res3 = _nmod_vec_init(len);
    _nmod_vec_inv2(res3, vec, len, mod);

    if (!_nmod_vec_equal(res1, res2, len))
        printf("res1 != res2\n");

    if (!_nmod_vec_equal(res3, res2, len))
        printf("res3 != res2\n");

    _nmod_vec_clear(res1);
    _nmod_vec_clear(res2);
    _nmod_vec_clear(res3);
    _nmod_vec_clear(vec);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("test running, various bitlengths, len from 1 to ~1000 (no error message means success)...\n");
    for (slong len = 1; len < 1000; len += 20)
    {
        printf("%ld..", len);
        fflush(stdout);
        for (slong repeat = 0; repeat < 100; repeat++)
        {
            check_nmod_vec_inv(len, 3, state);
            check_nmod_vec_inv(len, 10, state);
            check_nmod_vec_inv(len, 20, state);
            check_nmod_vec_inv(len, 29, state);
            check_nmod_vec_inv(len, 30, state);
            check_nmod_vec_inv(len, 31, state);
            check_nmod_vec_inv(len, 32, state);
            check_nmod_vec_inv(len, 40, state);
            check_nmod_vec_inv(len, 50, state);
            check_nmod_vec_inv(len, 60, state);
            check_nmod_vec_inv(len, 63, state);
        }
    }
    printf("\n");

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
