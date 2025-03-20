#include <flint/flint.h>
#include <time.h>
#include <flint/nmod_vec.h>
#include <flint/nmod.h>

#include "nmod_vec_extra.h"


//// small values for testing before launching test:
//#define TIME_THRES 0.002
//#define NB_ITER 10
#define TIME_THRES 0.2
#define NB_ITER 2500


/*--------------------------------------------------------------*/
/* flint's DOT                                                  */
/*--------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_flint_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res += _nmod_vec_dot(v1, v2, len, mod, params);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot(v1, v2, len, mod, params);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}


ulong time_nmod_vec_dot_product_flint_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v1[i], state, len, mod);
    }
    nn_ptr v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v2[i], state, len, mod);
    }
    ulong res = 0;

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res += _nmod_vec_dot(v1[i], v2[i], len, mod, params);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot(v1[i], v2[i], len, mod, params);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        _nmod_vec_clear(v1[i]);
        _nmod_vec_clear(v2[i]);
    }

    return res;
}


/*------------------------------------------------------------*/
/* FROM HERE, EXPERIMENTAL                                    */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* split, experimental                                        */
/*------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_split16_avx_v1_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v1[i], state, len, mod);
    }
    nn_ptr v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v2[i], state, len, mod);
    }
    ulong res[NB_ITER];

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v1(v1[0], v2[0], len, mod);
        ulong res_correct = _nmod_vec_dot(v1[0], v2[0], len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res[i] += _nmod_vec_dot_product_2_split16_avx_v1(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_2_split16_avx_v1(v1[i], v2[i], len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        _nmod_vec_clear(v1[i]);
        _nmod_vec_clear(v2[i]);
    }

    return 0;
}

ulong time_nmod_vec_dot_product_split16_avx_v1_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);

    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v1(v1, v2, len, mod);
        ulong res_correct = _nmod_vec_dot(v1, v2, len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res += _nmod_vec_dot_product_2_split16_avx_v1(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_2_split16_avx_v1(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_nmod_vec_dot_product_split16_avx_v2_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    uint * v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = flint_malloc(len * sizeof(uint));
        for (ulong j = 0; j < len; j++)
            v1[i][j] = (uint) n_randint(state, mod.n);
    }
    uint * v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = flint_malloc(len * sizeof(uint));
        for (ulong j = 0; j < len; j++)
            v2[i][j] = (uint) n_randint(state, mod.n);
    }
    ulong res[NB_ITER];

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v2(v1[0], v2[0], len, mod);
        nn_ptr v1long = _nmod_vec_init(len);
        nn_ptr v2long = _nmod_vec_init(len);
        for (ulong j = 0; j < len; j++)
        {
            v1long[j] = v1[0][j];
            v2long[j] = v2[0][j];
        }
        ulong res_correct = _nmod_vec_dot(v1long, v2long, len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res[i] += _nmod_vec_dot_product_2_split16_avx_v2(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_2_split16_avx_v2(v1[i], v2[i], len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        flint_free(v1[i]);
        flint_free(v2[i]);
    }

    return 0;
}

ulong time_nmod_vec_dot_product_split16_avx_v2_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    uint * v1 = flint_malloc(len * sizeof(uint));
    for (ulong j = 0; j < len; j++)
        v1[j] = (uint) n_randint(state, mod.n);
    uint * v2 = flint_malloc(len * sizeof(uint));
    for (ulong j = 0; j < len; j++)
        v2[j] = (uint) n_randint(state, mod.n);

    ulong res = 0;

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v2(v1, v2, len, mod);
        nn_ptr v1long = _nmod_vec_init(len);
        nn_ptr v2long = _nmod_vec_init(len);
        for (ulong j = 0; j < len; j++)
        {
            v1long[j] = v1[j];
            v2long[j] = v2[j];
        }
        ulong res_correct = _nmod_vec_dot(v1long, v2long, len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res += _nmod_vec_dot_product_2_split16_avx_v2(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_2_split16_avx_v2(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    flint_free(v1);
    flint_free(v2);

    return res;
}

ulong time_nmod_vec_dot_product_split16_avx_v3_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v1[i], state, len, mod);
    }
    nn_ptr v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v2[i], state, len, mod);
    }
    ulong res[NB_ITER];

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v3(v1[0], v2[0], len, mod);
        ulong res_correct = _nmod_vec_dot(v1[0], v2[0], len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res[i] += _nmod_vec_dot_product_2_split16_avx_v3(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_2_split16_avx_v3(v1[i], v2[i], len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        _nmod_vec_clear(v1[i]);
        _nmod_vec_clear(v2[i]);
    }

    return 0;
}

ulong time_nmod_vec_dot_product_split16_avx_v3_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);

    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v3(v1, v2, len, mod);
        ulong res_correct = _nmod_vec_dot(v1, v2, len, mod, params);
        if (res_split != res_correct)
        //if (0)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res += _nmod_vec_dot_product_2_split16_avx_v3(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_2_split16_avx_v3(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_nmod_vec_dot_product_split16_avx_v4_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    uint * v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = flint_malloc(len * sizeof(uint));
        for (ulong j = 0; j < len; j++)
            v1[i][j] = (uint) n_randint(state, mod.n);
    }
    uint * v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = flint_malloc(len * sizeof(uint));
        for (ulong j = 0; j < len; j++)
            v2[i][j] = (uint) n_randint(state, mod.n);
    }
    ulong res[NB_ITER];

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v4(v1[0], v2[0], len, mod);
        nn_ptr v1long = _nmod_vec_init(len);
        nn_ptr v2long = _nmod_vec_init(len);
        for (ulong j = 0; j < len; j++)
        {
            v1long[j] = v1[0][j];
            v2long[j] = v2[0][j];
        }
        ulong res_correct = _nmod_vec_dot(v1long, v2long, len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res[i] += _nmod_vec_dot_product_2_split16_avx_v4(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_2_split16_avx_v4(v1[i], v2[i], len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        flint_free(v1[i]);
        flint_free(v2[i]);
    }

    return 0;
}

ulong time_nmod_vec_dot_product_split16_avx_v4_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    uint * v1 = flint_malloc(len * sizeof(uint));
    for (ulong j = 0; j < len; j++)
        v1[j] = (uint) n_randint(state, mod.n);
    uint * v2 = flint_malloc(len * sizeof(uint));
    for (ulong j = 0; j < len; j++)
        v2[j] = (uint) n_randint(state, mod.n);

    ulong res = 0;

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16_avx_v4(v1, v2, len, mod);
        nn_ptr v1long = _nmod_vec_init(len);
        nn_ptr v2long = _nmod_vec_init(len);
        for (ulong j = 0; j < len; j++)
        {
            v1long[j] = v1[j];
            v2long[j] = v2[j];
        }
        ulong res_correct = _nmod_vec_dot(v1long, v2long, len, mod, params);
        //if (res_split != res_correct)
        if (0)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res += _nmod_vec_dot_product_2_split16_avx_v4(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_2_split16_avx_v4(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    flint_free(v1);
    flint_free(v2);

    return res;
}



ulong time_nmod_vec_dot_product_split26_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v1[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v1[i], state, len, mod);
    }
    nn_ptr v2[NB_ITER];
    for (slong i = 0; i < NB_ITER; i++)
    {
        v2[i] = _nmod_vec_init(len);
        _nmod_vec_rand(v2[i], state, len, mod);
    }
    ulong res[NB_ITER];

    { // TEST
        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        ulong res_split = _nmod_vec_dot_product_split26(v1[0], v2[0], len, mod);
        ulong res_correct = _nmod_vec_dot(v1[0], v2[0], len, mod, params);
        if (res_split != res_correct)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        for (slong i = 0; i < NB_ITER; i++) // warmup
            res[i] += _nmod_vec_dot_product_split26_avx(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_split26_avx(v1[i], v2[i], len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        _nmod_vec_clear(v1[i]);
        _nmod_vec_clear(v2[i]);
    }

    return 0;
}

ulong time_nmod_vec_dot_product_split26_cf(ulong len, ulong n, flint_rand_t state)
{
nmod_t mod;
nmod_init(&mod, n);

nn_ptr v1;
v1 = _nmod_vec_init(len);
_nmod_vec_rand(v1, state, len, mod);

nn_ptr v2;
v2 = _nmod_vec_init(len);
_nmod_vec_rand(v2, state, len, mod);

ulong res = 0;

{ // TEST
    const dot_params_t params = _nmod_vec_dot_params(len, mod);
    ulong res_split = _nmod_vec_dot_product_split26(v1, v2, len, mod);
    ulong res_correct = _nmod_vec_dot(v1, v2, len, mod, params);
    if (res_split != res_correct)
    {
        printf("\nDOT PRODUCT ERROR!\n");
        return 0;
    }
}

double t1;
clock_t tt;
long nb_iter;

t1 = 0.0; nb_iter = 0;
while (t1 < TIME_THRES)
{
    for (slong i = 0; i < NB_ITER; i++) // warmup
        res += _nmod_vec_dot_product_split26_avx(v1, v2, len, mod);

    tt = clock();
    for (slong i = 0; i < NB_ITER; i++)
        res += _nmod_vec_dot_product_split26_avx(v1, v2, len, mod);
    t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
    nb_iter += NB_ITER;
}
t1 /= nb_iter;
printf("%.1e\t", t1);

_nmod_vec_clear(v1);
_nmod_vec_clear(v2);

return res;
}


/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125);

    //const slong nlens = 33;
    //const slong lens[] = {2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 30, 40, 50, 65, 80, 100, 130, 160, 200, 260, 320, 400, 500, 600, 700, 800, 900, 1000, 2000, 4000, 8000, 16000};

    const slong nlens = 12;
    const slong lens[] = {50, 100, 200, 400, 600, 800, 1000, 2000, 4000, 8000, 16000, 50000};

    const slong nbits = 19;
    const slong bits[] = {17, 20, 23, 26, 29, 30, 31, 32, 33, 40, 50, 55, 57, 59, 60, 61, 62, 63, 64};

    const slong nfuns = 6;
    typedef ulong (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_nmod_vec_dot_product_flint_cf,      // 0
        time_nmod_vec_dot_product_flint_cu,      // 1
        time_nmod_vec_dot_product_split26_cf,    // 2
        time_nmod_vec_dot_product_split26_cu,    // 3
        time_nmod_vec_dot_product_split16_avx_v1_cf,    // 4
        time_nmod_vec_dot_product_split16_avx_v1_cu,    // 5
        time_nmod_vec_dot_product_split16_avx_v2_cf,    // 6
        time_nmod_vec_dot_product_split16_avx_v2_cu,    // 7
        time_nmod_vec_dot_product_split16_avx_v3_cf,    // 8
        time_nmod_vec_dot_product_split16_avx_v3_cu,    // 9
        time_nmod_vec_dot_product_split16_avx_v4_cf,    // 10
        time_nmod_vec_dot_product_split16_avx_v4_cu,    // 11
    };

    if (argc == 1)
    {
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];

            printf("bit/len");
            for (slong i = 0; i < nlens; i++)
                printf("\t%ld", lens[i]);
            printf("\n");

            for (slong j = 0; j < nbits; j++)
            {
                const slong b = bits[j];

                printf("%ldmid\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
                printf("\n");
            }
        }
    }
    else if (argc == 2)  // function is given
    {
        const timefun tfun = funs[atoi(argv[1])];

        printf("bit/len");
        for (slong i = 0; i < nlens; i++)
            printf("\t%ld", lens[i]);
        printf("\n");

        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];

            printf("%ldmid\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
            printf("\n");
        }
    }
    else if (argc == 3)  // function + nbits given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);

        printf("bit/len");
        for (slong i = 0; i < nlens; i++)
            printf("\t%ld", lens[i]);
        printf("\n");

        printf("%ldmid\t", b);
        for (slong i = 0; i < nlens; i++)
            tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");
    }
    else if (argc == 4)  // function + nbits + len given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);
        const slong len = atoi(argv[3]);

        printf("bit/len");
        printf("\t%ld", len);
        printf("\n");

        printf("%ldmid\t", b);
        tfun(len, (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
