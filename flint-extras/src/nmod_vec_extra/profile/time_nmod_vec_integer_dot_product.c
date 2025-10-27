#include <flint/flint.h>
#include <time.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void time_nmod_vec_integer_dot_product(ulong len, ulong n)
{
    flint_rand_t state;
    nn_ptr res = flint_malloc(3 * sizeof(ulong));
    nn_ptr v1, v2;
    nmod_t mod;
    double t;
    clock_t tt;
    long nb_iter;

    flint_rand_init(state);
    nmod_init(&mod, n);

    ulong maxbits = FLINT_BIT_COUNT(n);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.2)
    //while (t < 0.5 && nb_iter<2)
    {
        tt = clock();
        nmod_vec_integer_dot_product(res, v1, v2, len, maxbits, maxbits);
        nmod_vec_integer_dot_product(res, v1, v2, len, maxbits, maxbits);
        nmod_vec_integer_dot_product(res, v1, v2, len, maxbits, maxbits);
        nmod_vec_integer_dot_product(res, v1, v2, len, maxbits, maxbits);
        nmod_vec_integer_dot_product(res, v1, v2, len, maxbits, maxbits);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%.1e\t", t);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
    flint_free(res);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    printf("len\t3\t10\t20\t30\t31\t32\t40\t50\t60\n");
    for (slong len = 1; len < 1000; len += 21)
    {
        printf("%ld\t", len);
        time_nmod_vec_integer_dot_product(len, (1L << 2) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 9) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 19) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 29) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 30) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 31) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 39) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 49) + 1);
        time_nmod_vec_integer_dot_product(len, (1L << 59) + 1);
        printf("\n");
    }

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
