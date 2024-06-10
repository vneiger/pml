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

ulong time_nmod_vec_dot_product_flint_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

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
        ulong buf;
        slong j;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += _nmod_vec_dot(v1[i], v2[i], len, mod, n_limbs);
        {
            nn_srcptr v1i = v1[i];
            nn_srcptr v2i = v2[i];
            NMOD_VEC_DOT(buf, j, len, v1i[j], v2i[j], mod, n_limbs);
            res += buf;
        }


        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += _nmod_vec_dot(v1[i], v2[i], len, mod, n_limbs);
        {
            nn_srcptr v1i = v1[i];
            nn_srcptr v2i = v2[i];
            NMOD_VEC_DOT(buf, j, len, v1i[j], v2i[j], mod, n_limbs);
            res += buf;
        }
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

ulong time_nmod_vec_dot_product_flint_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);

    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        slong j;
        nn_srcptr v1s = v1;
        nn_srcptr v2s = v2;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += _nmod_vec_dot(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT(buf, j, len, v1s[j], v2s[j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += _nmod_vec_dot(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT(buf, j, len, v1s[j], v2s[j], mod, n_limbs);
            res += buf;
        }
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
/* general function                                             */
/*--------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

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

    { // TEST
        ulong ii;
        ulong res_new;
        nn_srcptr v1i = v1[0];
        nn_srcptr v2i = v2[0];
        NMOD_VEC_DOT_PRODUCT(res_new, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
        ulong res_flint = _nmod_vec_dot(v1[0], v2[0], len, mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1;
    clock_t tt;
    long nb_iter;

    nn_srcptr v1i, v2i;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }
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

ulong time_nmod_vec_dot_product_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);

    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    { // TEST
        ulong res_new;
        ulong ii;
        NMOD_VEC_DOT_PRODUCT(res_new, ii, len, v1[ii], v2[ii], mod, n_limbs);
        //res_new = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        ulong res_flint = _nmod_vec_dot(v1, v2, len, mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1[ii], v2[ii], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1[ii], v2[ii], mod, n_limbs);
            res += buf;
        }
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

/*------------------------------------------------------------*/
/* AVX2, supports more lengths / moduli                       */
/*------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_mod32_avx2_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    if (n > 1518500249)
        return 0;

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

    uint power2 = (1L << DOT_SP_NB) % n;

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_avx = _nmod_vec_dot_mod32_avx2(v1[0], v2[0], len, mod, power2);
        ulong res_lng = nmod_vec_dot_product(v1[0], v2[0], len, mod, n_limbs);
        if (res_avx != res_lng)
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
            res += _nmod_vec_dot_mod32_avx2(v1[i], v2[i], len, mod, power2);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_mod32_avx2(v1[i], v2[i], len, mod, power2);
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

ulong time_nmod_vec_dot_product_mod32_avx2_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    if (n > 1518500249)
        return 0;
    //printf("%ld\t", n_limbs);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    ulong power2 = (1L << DOT_SP_NB) % n;

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_avx = _nmod_vec_dot_mod32_avx2(v1, v2, len, mod, power2);
        ulong res_lng = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        if (res_avx != res_lng)
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
            res += _nmod_vec_dot_mod32_avx2(v1, v2, len, mod, power2);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_mod32_avx2(v1, v2, len, mod, power2);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}


/*------------------------------------------------------------*/
/* small mod, supports more lengths than nlimbs=1 versions    */
/*------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_mod32_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    if (n > 1518500249)
        return 0;

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

    uint power2 = (1L << DOT_SP_NB) % n;

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_mod32 = _nmod_vec_dot_mod32(v1[0], v2[0], len, mod, power2);
        ulong res_lng = nmod_vec_dot_product(v1[0], v2[0], len, mod, n_limbs);
        if (res_mod32 != res_lng)
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
            res += _nmod_vec_dot_mod32(v1[i], v2[i], len, mod, power2);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_mod32(v1[i], v2[i], len, mod, power2);
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

ulong time_nmod_vec_dot_product_mod32_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    if (n > 1518500249)
        return 0;
    //printf("%ld\t", n_limbs);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    ulong power2 = (1L << DOT_SP_NB) % n;

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_mod32 = _nmod_vec_dot_mod32(v1, v2, len, mod, power2);
        ulong res_lng = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        if (res_mod32 != res_lng)
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
            res += _nmod_vec_dot_mod32(v1, v2, len, mod, power2);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_mod32(v1, v2, len, mod, power2);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

/*------------------------------------------------------------*/
/* timing of global func against current flint                */
/*------------------------------------------------------------*/

ulong time_vs_current_flint_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

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

    { // TEST
        ulong ii;
        ulong res_new;
        nn_srcptr v1i = v1[0];
        nn_srcptr v2i = v2[0];
        NMOD_VEC_DOT_PRODUCT(res_new, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
        ulong res_flint = _nmod_vec_dot(v1[0], v2[0], len, mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tt;
    long nb_iter;

    nn_srcptr v1i, v2i;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;

    t2 = 0.0; nb_iter = 0;
    while (t2 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1[i], v2[i], len, mod, n_limbs);
        {
            v1i = v1[i]; v2i = v2[i];
            ulong ii;
            NMOD_VEC_DOT(buf, ii, len, v1i[ii], v2i[ii], mod, n_limbs);
            res += buf;
        }
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t2 /= nb_iter;

    printf("%.1e\t", t2/t1);

    for (slong i = 0; i < NB_ITER; i++)
    {
        _nmod_vec_clear(v1[i]);
        _nmod_vec_clear(v2[i]);
    }

    return res;
}

ulong time_vs_current_flint_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);

    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    { // TEST
        ulong res_new;
        ulong ii;
        NMOD_VEC_DOT_PRODUCT(res_new, ii, len, v1[ii], v2[ii], mod, n_limbs);
        //res_new = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        ulong res_flint = _nmod_vec_dot(v1, v2, len, mod, n_limbs);
        if (res_new != res_flint)
        {
            printf("\nDOT PRODUCT ERROR!\n");
            return 0;
        }
    }

    ulong res = 0;

    double t1, t2;
    clock_t tt;
    long nb_iter;

    t1 = 0.0; nb_iter = 0;
    while (t1 < TIME_THRES)
    {
        ulong buf;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1[ii], v2[ii], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            ulong ii;
            NMOD_VEC_DOT_PRODUCT(buf, ii, len, v1[ii], v2[ii], mod, n_limbs);
            res += buf;
        }
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;

    t2 = 0.0; nb_iter = 0;
    while (t2 < TIME_THRES)
    {
        ulong buf;
        slong j;
        nn_srcptr v1s = v1;
        nn_srcptr v2s = v2;
        for (slong i = 0; i < NB_ITER; i++) // warmup
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT(buf, j, len, v1s[j], v2s[j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            //res += nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        {
            NMOD_VEC_DOT(buf, j, len, v1s[j], v2s[j], mod, n_limbs);
            res += buf;
        }
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t2 /= nb_iter;

    printf("%.1e\t", t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}





/*------------------------------------------------------------*/
/* FROM HERE, EXPERIMENTAL                                    */
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* AVX2                                                       */
/*------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_avx2_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    if (n_limbs > 1)
        return 0;
    //printf("%ld\t", n_limbs);

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

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_avx = _nmod_vec_dot_product_1_avx2(v1[0], v2[0], len, mod);
        ulong res_lng = nmod_vec_dot_product(v1[0], v2[0], len, mod, n_limbs);
        if (res_avx != res_lng)
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
            res += _nmod_vec_dot_product_1_avx2(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_1_avx2(v1[i], v2[i], len, mod);
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

ulong time_nmod_vec_dot_product_avx2_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    if (n_limbs > 1)
        return 0;

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_avx = _nmod_vec_dot_product_1_avx2(v1, v2, len, mod);
        ulong res_lng = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        if (res_avx != res_lng)
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
            res += _nmod_vec_dot_product_1_avx2(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_1_avx2(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}



/*------------------------------------------------------------*/
/* AVX512i                                                    */
/*------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_avx512_cu(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    if (n_limbs > 1)
        return 0;

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

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_avx = _nmod_vec_dot_product_1_avx512(v1[0], v2[0], len, mod);
        ulong res_lng = nmod_vec_dot_product(v1[0], v2[0], len, mod, n_limbs);
        if (res_avx != res_lng)
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
            res += _nmod_vec_dot_product_1_avx512(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_1_avx512(v1[i], v2[i], len, mod);
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

ulong time_nmod_vec_dot_product_avx512_cf(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
    if (n_limbs > 1)
        return 0;

    nn_ptr v1;
    v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2;
    v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    ulong res = 0;

    { // TEST
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_avx = _nmod_vec_dot_product_1_avx512(v1, v2, len, mod);
        ulong res_lng = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        if (res_avx != res_lng)
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
            res += _nmod_vec_dot_product_1_avx512(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_1_avx512(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

/*------------------------------------------------------------*/
/* split, experimental                                        */
/*------------------------------------------------------------*/

ulong time_nmod_vec_dot_product_split16_cu(ulong len, ulong n, flint_rand_t state)
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
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16(v1[0], v2[0], len, mod);
        ulong res_correct = nmod_vec_dot_product(v1[0], v2[0], len, mod, n_limbs);
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
            res[i] += _nmod_vec_dot_product_2_split16(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_2_split16(v1[i], v2[i], len, mod);
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

ulong time_nmod_vec_dot_product_split16_cf(ulong len, ulong n, flint_rand_t state)
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
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        ulong res_correctj = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        if (res_split != res_correctj)
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
            res += _nmod_vec_dot_product_2_split16(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_2_split16(v1, v2, len, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t1 /= nb_iter;
    printf("%.1e\t", t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

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
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split26(v1[0], v2[0], len, mod);
        ulong res_correct = nmod_vec_dot_product(v1[0], v2[0], len, mod, n_limbs);
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
            res[i] += _nmod_vec_dot_product_2_split26(v1[i], v2[i], len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res[i] += _nmod_vec_dot_product_2_split26(v1[i], v2[i], len, mod);
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
        const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);
        ulong res_split = _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
        ulong res_correctj = nmod_vec_dot_product(v1, v2, len, mod, n_limbs);
        if (res_split != res_correctj)
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
            res += _nmod_vec_dot_product_2_split26(v1, v2, len, mod);

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
            res += _nmod_vec_dot_product_2_split26(v1, v2, len, mod);
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

    const slong nlens = 33;
    const slong lens[] = {2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 30, 40, 50, 65, 80, 100, 130, 160, 200, 260, 320, 400, 500, 600, 700, 800, 900, 1000, 2000, 4000, 8000, 16000};

    const slong nbits = 20;
    const slong bits[] = {3, 10, 20, 23, 26, 29, 30, 31, 32, 33, 40, 50, 55, 57, 59, 60, 61, 62, 63, 64};

    const slong nfuns = 10;
    typedef ulong (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_nmod_vec_dot_product_cf,            // 0
        time_nmod_vec_dot_product_mod32_cf,      // 1
        time_nmod_vec_dot_product_mod32_avx2_cf, // 2
        time_nmod_vec_dot_product_flint_cf,      // 3
        time_nmod_vec_dot_product_cu,            // 4
        time_nmod_vec_dot_product_mod32_cu,      // 5
        time_nmod_vec_dot_product_mod32_avx2_cu, // 6
        time_nmod_vec_dot_product_flint_cu,      // 7
        time_vs_current_flint_cf,                // 8
        time_vs_current_flint_cu,                // 9
    };

    if (argc == 1)
    {
        for (slong ifun = 8; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];

            printf("bit/len");
            for (slong i = 0; i < nlens; i++)
                printf("\t%ld", lens[i]);
            printf("\n");

            for (slong j = 0; j < nbits; j++)
            {
                const slong b = bits[j];

                printf("%ldmin\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
                printf("\n");

                printf("%ldmid\t", b);
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
                printf("\n");

                printf("%ldmax\t", b);
                if (b < 64)
                    for (slong i = 0; i < nlens; i++)
                        tfun(lens[i], (UWORD(1) << b) - 1, state);
                else
                    for (slong i = 0; i < nlens; i++)
                        tfun(lens[i], UWORD_MAX, state);
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

            printf("%ldmin\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
            printf("\n");

            printf("%ldmid\t", b);
            for (slong i = 0; i < nlens; i++)
                tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
            printf("\n");

            printf("%ldmax\t", b);
            if (b < 64)
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], (UWORD(1) << b) - 1, state);
            else
                for (slong i = 0; i < nlens; i++)
                    tfun(lens[i], UWORD_MAX, state);
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

        //printf("%ldmin\t", b);
        //for (slong i = 0; i < nlens; i++)
        //    tfun(lens[i], (UWORD(1) << (b-1)) + 1, state);
        //printf("\n");

        printf("%ldmid\t", b);
        for (slong i = 0; i < nlens; i++)
            tfun(lens[i], (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");

        //printf("%ldmax\t", b);
        //if (b < 64)
        //    for (slong i = 0; i < nlens; i++)
        //        tfun(lens[i], (UWORD(1) << b) - 1, state);
        //else
        //    for (slong i = 0; i < nlens; i++)
        //        tfun(lens[i], UWORD_MAX, state);
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

        //printf("%ldmin\t", b);
        //tfun(len, (UWORD(1) << (b-1)) + 1, state);
        //printf("\n");

        printf("%ldmid\t", b);
        tfun(len, (UWORD(1) << (b-1)) + (UWORD(1) << (b-2)), state);
        printf("\n");

        //printf("%ldmax\t", b);
        //if (b < 64)
        //    tfun(len, (UWORD(1) << b) - 1, state);
        //else
        //    tfun(len, UWORD_MAX, state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
