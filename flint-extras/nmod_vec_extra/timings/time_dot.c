#include <stdlib.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/machine_vectors.h>

#include "dot.h"

// small values for testing before launching test:
//#define TIME_THRES 0.02
//#define NB_ITER 100
// full test:
#define TIME_THRES 0.1
#define NB_ITER 5000
//#define SMALL_SUITE 1

// utility
static inline
void _nmod_vec_rand(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = n_randint(state, mod.n);
}

// uniform random
static inline
void nmod_mat_rand(nmod_mat_t mat, flint_rand_t state)
{
    _nmod_vec_rand(mat->entries, state, mat->r * mat->c, mat->mod);
}


/*--------------------------------------------------------*/
/* timing of dot-expr against current flint               */
/*--------------------------------------------------------*/

ulong time_dotexpr_vs_flint_plain(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    nn_srcptr v1i = v1;
    nn_srcptr v2i = v2;

    { // TEST
        ulong res_new;
        ulong j;
        res_new = _nmod_vec_newdot(v1, v2, len, mod, params);
        ulong res_flint;
        v1i = v1; v2i = v2;
        FLINT_NMOD_VEC_DOT(res_flint, j, len, v1i[j], v2i[j], mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tstart, tend;
    long nb_iter;

    t1 = 0.0; t2 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        ulong j;

        tstart = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_newdot(v1, v2, len, mod, params);
        tend = clock();
        t1 += (double)(tend-tstart) / CLOCKS_PER_SEC;

        tstart = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            FLINT_NMOD_VEC_DOT(buf, j, len, v1i[j], v2i[j], mod, n_limbs);
            res += buf;
        }
        tend = clock();
        t2 += (double)(tend-tstart) / CLOCKS_PER_SEC;

        nb_iter += NB_ITER;
    }

    printf("%.1e\t%.1e\t%.1e\t", t1/nb_iter, t2/nb_iter, t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_dotexpr_vs_flint_rev(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    nn_srcptr v1i, v2i;
    v1i = v1; v2i = v2;

    { // TEST
        ulong res_new;
        ulong j;
        res_new = _nmod_vec_newdot_rev(v1, v2, len, mod, params);
        ulong res_flint;
        FLINT_NMOD_VEC_DOT(res_flint, j, len, v1i[j], v2i[len - 1 - j], mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tstart, tend;
    long nb_iter;

    t2 = 0.0; t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        tstart = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_newdot_rev(v1, v2, len, mod, params);
        tend = clock();
        t1 += (double)(tend-tstart) / CLOCKS_PER_SEC;

        ulong buf;
        ulong j;
        tstart = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            v1i = v1; v2i = v2;
            FLINT_NMOD_VEC_DOT(buf, j, len, v1i[j], v2i[len - 1 - j], mod, n_limbs);
            res += buf;
        }
        tend = clock();
        t2 += (double)(tend-tstart) / CLOCKS_PER_SEC;

        nb_iter += NB_ITER;
    }

    printf("%.1e\t%.1e\t%.1e\t", t1/nb_iter, t2/nb_iter, t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_dotexpr_vs_flint_rev2(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    nn_srcptr v1i, v2i;
    v1i = v1; v2i = v2;

    { // TEST
        ulong res_new;
        ulong j;
        _NMOD_VEC_DOT_NEW(res_new, j, len/9, v1[len - 1 - 9*j], v2[len - 1 - 9*j], mod, params);
        ulong res_flint;
        FLINT_NMOD_VEC_DOT(res_flint, j, len/9, v1i[len - 1 - 9*j], v2i[len - 1 - 9*j], mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tstart, tend;
    long nb_iter;

    t2 = 0.0; t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        ulong j;

        tstart = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            _NMOD_VEC_DOT_NEW(buf, j, len/9, v1[len - 1 - 9*j], v2[len - 1 - 9*j], mod, params);
            res += buf;
        }
        tend = clock();
        t1 += (double)(tend-tstart) / CLOCKS_PER_SEC;

        tstart = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            v1i = v1; v2i = v2;
            NMOD_VEC_DOT(buf, j, len/9, v1i[len - 1 - 9*j], v2i[len - 1 - 9*j], mod, n_limbs);
            res += buf;
        }
        tend = clock();
        t2 += (double)(tend-tstart) / CLOCKS_PER_SEC;

        nb_iter += NB_ITER;
    }

    printf("%.1e\t%.1e\t%.1e\t", t1/nb_iter, t2/nb_iter, t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_dotexpr_vs_flint_matmul(ulong len, ulong n, flint_rand_t state)
{
    if (len > 2000)
    {
        printf("inf\tinf\t1.0\t");
        return 0;
    }
    else
    {
        { // TEST
            nmod_mat_t mat1, mat2, mat, tmp;
            nmod_mat_init(mat1, len, len, n);
            nmod_mat_rand(mat1, state);
            nmod_mat_init(mat2, len, len, n);
            nmod_mat_rand(mat2, state);
            nmod_mat_init(mat, len, len, n);
            nmod_mat_init(tmp, len, len, n);

            nmod_mat_mul_flint(mat, mat1, mat2);
            nmod_mat_mul_newdot(tmp, mat1, mat2);
            if (!nmod_mat_equal(mat, tmp))
                printf("\n\n\nMATMUL ERROR!\n\n\n");

            nmod_mat_clear(mat1);
            nmod_mat_clear(mat2);
            nmod_mat_clear(mat);
            nmod_mat_clear(tmp);
        }

        double t1, t2;
        clock_t tt;
        long nb_iter;

        const ulong NB_MAT_ITER = FLINT_MAX(1, NB_ITER / len);

        t1 = 0.0; t2 = 0.0; nb_iter = 0;
        while (t1 < TIME_THRES)
        {
            nmod_mat_t mat1, mat2, mat;
            nmod_mat_init(mat1, len, len, n);
            nmod_mat_rand(mat1, state);
            nmod_mat_init(mat2, len, len, n);
            nmod_mat_rand(mat2, state);
            nmod_mat_init(mat, len, len, n);

            tt = clock();
            for (ulong i = 0; i < NB_MAT_ITER; i++)
                nmod_mat_mul_newdot(mat, mat1, mat2);
            t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;

            tt = clock();
            for (ulong i = 0; i < NB_MAT_ITER; i++)
                nmod_mat_mul_flint(mat, mat1, mat2);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;

            nb_iter += NB_MAT_ITER;

            nmod_mat_clear(mat1);
            nmod_mat_clear(mat2);
            nmod_mat_clear(mat);
        }

        printf("%.1e\t%.1e\t%.1e\t", t1/nb_iter, t2/nb_iter, t2/t1);
        return 0;
    }
}


/*--------------------------------------------------------------*/
/* main                                                         */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

#if SMALL_SUITE
    const slong nlens = 4;
    const ulong lens[] = {2, 20, 200, 2000};

    const slong nbits = 10;
    const ulong bits[] = {20, 28, 30, 32, 33, 50, 60, 62, 63, 64};
#else
    const slong nlens = 13;
    const ulong lens[] = {5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 100000, 1000000};

    const slong nbits = 19;
    const ulong bits[] = {17, 20, 23, 26, 29, 30, 31, 32, 33, 40, 50, 55, 57, 59, 60, 61, 62, 63, 64};
#endif

    const slong nfuns = 3;
    typedef ulong (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_dotexpr_vs_flint_plain,                   // 0
        time_dotexpr_vs_flint_rev,                     // 1
        time_dotexpr_vs_flint_rev2,                    // 2
        time_dotexpr_vs_flint_matmul,                  // 3
    };

    printf("warmup...");
    {
        ulong res = 0;
        for (long i = 0; i < NB_ITER; i++)
        {
            const slong len = 10000;
            nmod_t mod;
            nmod_init(&mod, (1<<20) + 1);
            const dot_params_t params = _nmod_vec_dot_params(len, mod);

            nn_ptr v1;
            v1 = _nmod_vec_init(len);
            _nmod_vec_rand(v1, state, len, mod);
            nn_ptr v2;
            v2 = _nmod_vec_init(len);
            _nmod_vec_rand(v2, state, len, mod);

            ulong buf;
            ulong j;
            _NMOD_VEC_DOT_NEW(buf, j, len, v1[j], v2[j], mod, params);
            res += buf;
        }
        printf(" done (res = %ld)\n", res);
    }

    if (argc == 1)
    {
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];

            printf("bit/len");
            for (slong i = 0; i < nlens; i++)
                printf("\t%ld\t\t", lens[i]);
            printf("\n");

            for (slong j = 0; j < nbits; j++)
            {
                const slong b = bits[j];

#ifndef SMALL_SUITE
                printf("%ldmin\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
                printf("\n");
#endif

                printf("%ldmid\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
                printf("\n");

#ifndef SMALL_SUITE
                printf("%ldmax\t", b);
                if (b < 64)
                    for (slong i = 0; i < nlens; i++)
                        tfun(lens[i], (UWORD(1) << b) - 1, state);
                else
                    for (slong i = 0; i < nlens; i++)
                        tfun(lens[i], UWORD_MAX, state);
                printf("\n");
#endif
            }
        }
    }
    else if (argc == 2)  // function is given
    {
        const timefun tfun = funs[atoi(argv[1])];

        printf("bit/len");
        for (slong i = 0; i < nlens; i++)
            printf("\t%ld\t\t", lens[i]);
        printf("\n");

        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];

#ifndef SMALL_SUITE
            printf("%ldmin\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
            printf("\n");
#endif

            printf("%ldmid\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
            printf("\n");

#ifndef SMALL_SUITE
            printf("%ldmax\t", b);
            if (b < 64)
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << b) - 1, state);
            else
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], UWORD_MAX, state);
            printf("\n");
#endif
        }
    }
    else if (argc == 3)  // function + nbits given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);

        printf("bit/len");
        for (slong i = 0; i < nlens; i++)
            printf("\t%ld\t\t", lens[i]);
        printf("\n");

#ifndef SMALL_SUITE
        printf("%ldmin\t", b);
        for (slong i = 0; i < nlens; i++)
            tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
        printf("\n");
#endif

        printf("%ldmid\t", b);
        for (slong i = 0; i < nlens; i++)
            tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");

#ifndef SMALL_SUITE
        printf("%ldmax\t", b);
        if (b < 64)
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << b) - 1, state);
        else
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], UWORD_MAX, state);
        printf("\n");
#endif
    }
    else if (argc == 4)  // function + nbits + len given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);
        const slong len = atoi(argv[3]);

        printf("bit/len");
        printf("\t%ld", len);
        printf("\n");

#ifndef SMALL_SUITE
        printf("%ldmin\t", b);
        tfun(len, (UWORD(1) << (b-1)) + 1, state);
        printf("\n");
#endif

        printf("%ldmid\t", b);
        tfun(len, (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");

#ifndef SMALL_SUITE
        printf("%ldmax\t", b);
        if (b < 64)
            tfun(len, (UWORD(1) << b) - 1, state);
        else
            tfun(len, UWORD_MAX, state);
        printf("\n");
#endif
    }

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
