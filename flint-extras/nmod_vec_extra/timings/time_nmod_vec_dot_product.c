#include <assert.h>
#include <flint/flint.h>
#include <time.h>
#include <flint/nmod_vec.h>
#include <flint/nmod.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void time_nmod_vec_dot_product(ulong len, ulong maxbits1, ulong maxbits2, ulong n, flint_rand_t state)
{
    nmod_t mod, mod1, mod2;
    nmod_init(&mod, n);
    if (maxbits1 < FLINT_BITS) nmod_init(&mod1, UWORD(1) << maxbits1);
    else nmod_init(&mod1, UWORD_MAX);
    if (maxbits2 < FLINT_BITS) nmod_init(&mod2, UWORD(1) << maxbits2);
    else nmod_init(&mod2, UWORD_MAX);
    
    mp_ptr v1, v2;
    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    double t1, t2;
    clock_t tt;
    long nb_iter;

    mp_limb_t val1, val2;

    t1 = 0.0; nb_iter = 0;
    while (t1 < 0.2)
    {
        _nmod_vec_rand(v1, state, len, mod1);
        _nmod_vec_rand(v2, state, len, mod2);
        tt = clock();
        nmod_vec_dot_product(v1, v2, len, maxbits1, maxbits2, mod);
        nmod_vec_dot_product(v1, v2, len, maxbits1, maxbits2, mod);
        nmod_vec_dot_product(v1, v2, len, maxbits1, maxbits2, mod);
        nmod_vec_dot_product(v1, v2, len, maxbits1, maxbits2, mod);
        nmod_vec_dot_product(v1, v2, len, maxbits1, maxbits2, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    if (FLINT_BIT_COUNT(n) <= 31)
    {
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.2)
        {
            _nmod_vec_rand(v1, state, len, mod);
            _nmod_vec_rand(v2, state, len, mod);
            tt = clock();
            _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
            _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
            _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
            _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
            _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 5;
        }
        t2 /= nb_iter;
        printf("%.1e\t", t2);
    }

    t2 = 0.0; nb_iter = 0;
    while (t2 < 0.2)
    {
        _nmod_vec_rand(v1, state, len, mod);
        _nmod_vec_rand(v2, state, len, mod);
        tt = clock();
        _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t2 /= nb_iter;
    printf("%.1e\t", t2);


    // VERSUS FLINT:
    //int nlimbs = _nmod_vec_dot_bound_limbs(len, mod);
    //t2 = 0.0;
    //nb_iter = 0;
    //while (t2 < 0.2)
    ////while (t2 < 0.5 && nb_iter<2)
    //{
    //    tt = clock();
    //    val2 = _nmod_vec_dot(v1, v2, len, mod, nlimbs);
    //    val2 = _nmod_vec_dot(v1, v2, len, mod, nlimbs);
    //    val2 = _nmod_vec_dot(v1, v2, len, mod, nlimbs);
    //    val2 = _nmod_vec_dot(v1, v2, len, mod, nlimbs);
    //    val2 = _nmod_vec_dot(v1, v2, len, mod, nlimbs);
    //    t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    nb_iter += 5;
    //}
    ////t = 1000 * t;
    //t2 /= nb_iter;
    ////printf("%.1e\t", t);

    //printf("%.1e\t", t1/t2);

    val1 = nmod_vec_dot_product(v1, v2, len, maxbits1, maxbits2, mod);
    if (FLINT_BIT_COUNT(n) <= 31)
    {
        val2 = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        assert (val1 == val2 && "2_split16");
    }
    if (FLINT_BIT_COUNT(n) <= 52)
    {
        val2 = _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        if (val1 != val2) {fflush(stdout); printf("\n");}
        assert (val1 == val2 && "2_split26");
    }

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

    printf("len\t3\t3\t10\t10\t20\t20\t29\t29\t30\t30\t31\t31\t32\t32\t40\t40\t50\t50\t60\t60\t64\t64\n");
    //for (slong len = 100; len < 1000; len += 21)
    for (slong len = 800; len < 4000; len *= 2)
    {
        printf("%ld\t", len);
        //time_nmod_vec_dot_product(len, 1, 3, (UWORD(1) << 3) + 1, state);
        //time_nmod_vec_dot_product(len, 3, 3, (UWORD(1) << 3) + 1, state);
        //time_nmod_vec_dot_product(len, 5, 10, (UWORD(1) << 10) + 1, state);
        //time_nmod_vec_dot_product(len, 10, 10, (UWORD(1) << 10) + 1, state);
        //time_nmod_vec_dot_product(len, 10, 20, (UWORD(1) << 20) + 1, state);
        time_nmod_vec_dot_product(len, 20, 20, (UWORD(1) << 20) - 1, state);
        printf("\t\t");
        //time_nmod_vec_dot_product(len, 15, 29, (UWORD(1) << 29) + 1, state);
        //time_nmod_vec_dot_product(len, 28, 28, (UWORD(1) << 27) + 1, state);
        //printf("\t\t");
        //time_nmod_vec_dot_product(len, 15, 30, (UWORD(1) << 30) + 1, state);
        //time_nmod_vec_dot_product(len, 30, 30, (UWORD(1) << 30) + 1, state);
        //time_nmod_vec_dot_product(len, 15, 31, (UWORD(1) << 31) + 1, state);
        time_nmod_vec_dot_product(len, 31, 31, (UWORD(1) << 31) - 1, state);
        printf("\t\t");
        //time_nmod_vec_dot_product(len, 16, 32, (UWORD(1) << 32) + 1, state);
        //time_nmod_vec_dot_product(len, 32, 32, (UWORD(1) << 32) + 1, state);
        //time_nmod_vec_dot_product(len, 20, 40, (UWORD(1) << 40) + 1, state);
        time_nmod_vec_dot_product(len, 40, 40, (UWORD(1) << 40) + 1, state);
        //printf("\n");
        //time_nmod_vec_dot_product(len, 25, 50, (UWORD(1) << 50) + 1, state);
        //time_nmod_vec_dot_product(len, 50, 50, (UWORD(1) << 50) - 1, state);
        printf("\n");
        //time_nmod_vec_dot_product(len, 30, 60, (UWORD(1) << 60) + 1, state);
        //time_nmod_vec_dot_product(len, 60, 60, (UWORD(1) << 60) + 1, state);
        //time_nmod_vec_dot_product(len, 32, 64, UWORD_MAX, state);
        //time_nmod_vec_dot_product(len, 64, 64, UWORD_MAX, state);
        //printf("\n");
    }

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
