#include <time.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_extra.h"  // for rand
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void time_nmod_vec_dot_product_multi(ulong len, ulong k, ulong maxbits1, ulong maxbits2, ulong n, flint_rand_t state)
{
    // moduli
    nmod_t mod, mod1, mod2;
    nmod_init(&mod, n);

    if (maxbits1 < FLINT_BITS) nmod_init(&mod1, UWORD(1) << maxbits1);
    else nmod_init(&mod1, UWORD_MAX);
    if (maxbits2 < FLINT_BITS) nmod_init(&mod2, UWORD(1) << maxbits2);
    else nmod_init(&mod2, UWORD_MAX);

    // input vector
    mp_ptr u = _nmod_vec_init(len);

    // input collection of vectors (= matrices)
    mp_ptr * v = flint_malloc(len * sizeof(mp_ptr));
    for (ulong i = 0; i < len; i++)
        v[i] = _nmod_vec_init(k);

    mp_ptr ur = _nmod_vec_init(len);
    nmod_mat_t vmat; nmod_mat_init(vmat, len, k, n);

    // fill entries at random
    _nmod_vec_rand(u, state, len, mod1);
    for (ulong i = 0; i < len; i++)
        _nmod_vec_rand(v[i], state, k, mod2);

    for (ulong i = 0; i < len; i++)
        ur[i] = u[i] % n;
    nmod_mat_rand(vmat, state);

    double t1, t2;
    clock_t tt;
    long nb_iter;

    t1 = 0.0;
    nb_iter = 0;
    while (t1 < 0.2)
    //while (t1 < 0.5 && nb_iter<2)
    {
        mp_ptr uv = _nmod_vec_init(k);
        tt = clock();
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _nmod_vec_clear(uv);
        nb_iter += 5;
    }
    //t = 1000 * t;
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    t2 = 0.0;
    nb_iter = 0;
    while (t2 < 0.2)
    //while (t2 < 0.5 && nb_iter<2)
    {
        mp_ptr uv = _nmod_vec_init(k);
        tt = clock();
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _nmod_vec_clear(uv);
        nb_iter += 5;
    }
    //t = 1000 * t;
    t2 /= nb_iter;
    printf("%.1e\t", t2);

    _nmod_vec_clear(u);
    for (ulong i = 0; i < len; i++)
        _nmod_vec_clear(v[i]);
    flint_free(v);
    nmod_mat_clear(vmat);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_randinit(state);

    // launch full suite
    if (argc == 1)
    {
        const ulong mods[11] =
        {
                (UWORD(1) << 3) + 1,
                (UWORD(1) << 10) + 1,
                (UWORD(1) << 20) + 1,
                (UWORD(1) << 25) + 1,
                (UWORD(1) << 29) + 1,
                (UWORD(1) << 30) + 1,
                (UWORD(1) << 31) + 1,
                (UWORD(1) << 40) + 1,
                (UWORD(1) << 50) + 1,
                (UWORD(1) << 60) + 1,
                UWORD_MAX
        };
        const ulong nbits[11] = { 4, 11, 21, 26, 30, 31, 32, 41, 51, 61, 64 };

        const ulong dims[24] = {1,2,3,4,5,8,12,16,20,25,30,40,50,75,100,150,200,300,400,500,750,1000,1500,2000};

        printf("nbits\tbits1\tbits2\tlen\tk\tmulti\tflint\ttrspdot\n");
        for (slong i = 0; i < 11; i++)
        {
            //for (slong r = 0; r < 24; r++)
            //{
            //    for (slong c = 0; c < 24; c++)
            //    {
            //        printf("%ld\t%ld\t%ld\t%ld\t%ld\t", nbits[i]/2, nbits[i], nbits[i], dims[r], dims[c]);
            //        time_nmod_vec_dot_product_multi(dims[r], dims[c], nbits[i]/2, nbits[i], mods[i], state);
            //        printf("\n");
            //    }
            //}
            for (slong r = 0; r < 24; r++)
            {
                for (slong c = 0; c < 24; c++)
                {
                    printf("%ld\t%ld\t%ld\t%ld\t%ld\t", nbits[i], nbits[i], nbits[i], dims[r], dims[c]);
                    time_nmod_vec_dot_product_multi(dims[r], dims[c], nbits[i], nbits[i], mods[i], state);
                    printf("\n");
                }
            }
        }
    }

    // fixed number of bits
    else if (argc == 2)
    {
        printf("nbits\tbits1\tbits2\tlen\tk\tmulti\tflint\ttrspdot\n");
        const ulong nbit = atoi(argv[1]);
        const ulong dims[24] = {1,2,3,4,5,8,12,16,20,25,30,40,50,75,100,150,200,300,400,500,750,1000,1500,2000};
        for (slong r = 0; r < 24; r++)
        {
            for (slong c = 0; c < 24; c++)
            {
                printf("%ld\t%ld\t%ld\t%ld\t%ld\t", nbit, nbit, nbit, dims[r], dims[c]);
                time_nmod_vec_dot_product_multi(dims[r], dims[c], nbit, nbit, (UWORD(1) << (nbit-1))+1, state);
                printf("\n");
            }
        }
    }

    // fixed number of bits and dimensions
    else if (argc == 4)
    {
        printf("nbits\tbits1\tbits2\tlen\tk\tmulti\tflint\ttrspdot\n");
        const ulong nbit = atoi(argv[1]);
        const ulong len = atoi(argv[2]);
        const ulong k = atoi(argv[3]);
        printf("%ld\t%ld\t%ld\t%ld\t%ld\t", nbit, nbit, nbit, len, k);
        time_nmod_vec_dot_product_multi(len, k, nbit, nbit, (UWORD(1) << (nbit-1))+1, state);
        printf("\n");
    }

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
