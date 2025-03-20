#include <time.h>
#include <stdlib.h> //GV
#include <flint/flint.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/* creates and prints uniform random vector                   */
/*------------------------------------------------------------*/
int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("Usage %s [positive modulus] [vector length]\n",argv[0]);
        return 0;
    }

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    nmod_t mod;
    nmod_init(&mod, atol(argv[1]));

    const slong len = atol(argv[2]);

    nn_ptr vec_uni = _nmod_vec_init(len);
    nn_ptr vec_test = _nmod_vec_init(len);

    _nmod_vec_rand(vec_uni, state, len, mod);
    _nmod_vec_randtest(vec_test, state, len, mod);

    _nmod_vec_print_pretty(vec_uni, len, mod);
    _nmod_vec_print_pretty(vec_test, len, mod);

    _nmod_vec_clear(vec_uni);
    _nmod_vec_clear(vec_test);
    flint_rand_clear(state);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
