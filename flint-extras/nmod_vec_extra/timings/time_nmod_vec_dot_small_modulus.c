#include <time.h>
#include <assert.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n of nbits bits    */
/*--------------------------------------------------------------*/
void time_nmod_vec_dot_small_modulus(ulong len, ulong nbits, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    mp_ptr v1 = _nmod_vec_init(len);
    mp_ptr v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);

    // storing results
    mp_limb_t res1, res2, res3, res4;

    double t;
    clock_t tt;
    long nb_iter;

    t = 0.0;
    nb_iter = 0;
    while (t < 0.2)
    //while (t < 0.2 && nb_iter < 3)
    {
        tt = clock();
        res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
        res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
        res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
        res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
        res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.2)
    //while (t < 0.2 && nb_iter < 3)
    {
        tt = clock();
        res2 = nmod_vec_dot_product(v1, v2, len, nbits, nbits, mod);
        res2 = nmod_vec_dot_product(v1, v2, len, nbits, nbits, mod);
        res2 = nmod_vec_dot_product(v1, v2, len, nbits, nbits, mod);
        res2 = nmod_vec_dot_product(v1, v2, len, nbits, nbits, mod);
        res2 = nmod_vec_dot_product(v1, v2, len, nbits, nbits, mod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);


    t = 0.0;
    nb_iter = 0;
    while (t < 0.2)
    //while (t < 0.5 && nb_iter < 3)
    {
        tt = clock();
        res3 = _nmod_vec_dot_small_modulus(v1, v2, len, (1L << 45) % n, n, 1 / (double) n);
        res3 = _nmod_vec_dot_small_modulus(v1, v2, len, (1L << 45) % n, n, 1 / (double) n);
        res3 = _nmod_vec_dot_small_modulus(v1, v2, len, (1L << 45) % n, n, 1 / (double) n);
        res3 = _nmod_vec_dot_small_modulus(v1, v2, len, (1L << 45) % n, n, 1 / (double) n);
        res3 = _nmod_vec_dot_small_modulus(v1, v2, len, (1L << 45) % n, n, 1 / (double) n);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.2)
    //while (t < 0.5 && nb_iter < 3)
    {
        tt = clock();
        res4 = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        res4 = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        res4 = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        res4 = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        res4 = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    assert (res1 == res2);
    assert (res1 == res3);
    assert (res1 == res4);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_randinit(state);

    //printf("len\t4\t\t\t\t10\t\t\t\t20\t\t\t\t25\t\t\t\t29\t\t\t\t30\n");
    printf("len\t4\t\t\t\t25\t\t\t\t30\n");
    for (slong len = 1; len < 1000; len += 21)
    {
        printf("%ld\t", len);
        time_nmod_vec_dot_small_modulus(len, 4, (UWORD(1) << 3) + 1, state);
        //time_nmod_vec_dot_small_modulus(len, 10, (UWORD(1) << 9) + 1, state);
        //time_nmod_vec_dot_small_modulus(len, 20, (UWORD(1) << 19) + 1, state);
        time_nmod_vec_dot_small_modulus(len, 25, (UWORD(1) << 24) + 1, state);
        //time_nmod_vec_dot_small_modulus(len, 29, (UWORD(1) << 28) + 1, state);
        time_nmod_vec_dot_small_modulus(len, 30, (UWORD(1) << 29) + 1, state);
        printf("\n");
    }

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
