#include <stdlib.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/machine_vectors.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>

// small values for testing before launching test:
//#define TIME_THRES 0.002
//#define NB_ITER 10
//#define NB_MAT_ITER 1
// full test:
#define TIME_THRES 0.2
#define NB_ITER 2500
#define NB_MAT_ITER 3
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


// new general dot macro
#define DOT_SP_NB 56
#define _NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)        \
do                                                                   \
{                                                                    \
    res = UWORD(0);                                                  \
                                                                     \
    if (nlimbs == 1)                                                 \
    {                                                                \
        for (i = 0; i < len; i++)                                    \
            res += (expr1) * (expr2);                                \
                                                                     \
        NMOD_RED(res, res, mod);                                     \
    }                                                                \
    else if (mod.n <= 1515531528 && len <= 134744072)                \
    {                                                                \
        const ulong low_bits = (1L << DOT_SP_NB) - 1;                \
        const uint red_pow = (1L << DOT_SP_NB) % mod.n;              \
        ulong dp_lo = 0;                                             \
        uint dp_hi = 0;                                              \
                                                                     \
        for (i = 0; i+7 < len; )                                     \
        {                                                            \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
            dp_lo += (expr1) * (expr2);                              \
            i++;                                                     \
                                                                     \
            dp_hi += dp_lo >> DOT_SP_NB;                             \
            dp_lo &= low_bits;                                       \
        }                                                            \
                                                                     \
        for ( ; i < len; i++)                                        \
            dp_lo += (expr1) * (expr2);                              \
                                                                     \
        NMOD_RED(res, ((ulong)red_pow * dp_hi) + dp_lo, mod);        \
    }                                                                \
    else if (nlimbs == 2)                                            \
    {                                                                \
        ulong s0, s1;                                                \
        ulong u0 = UWORD(0);                                         \
        ulong u1 = UWORD(0);                                         \
                                                                     \
        for (i = 0; i+7 < len; )                                     \
        {                                                            \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
            i++;                                                     \
        }                                                            \
        for (; i < len; i++)                                         \
        {                                                            \
            umul_ppmm(s1, s0, (expr1), (expr2));                     \
            add_ssaaaa(u1, u0, u1, u0, s1, s0);                      \
        }                                                            \
                                                                     \
        NMOD2_RED2(res, u1, u0, mod);                                \
    }                                                                \
                                                                     \
    else if (nlimbs == 3)                                            \
    {                                                                \
        ulong s0, s1;                                                \
        ulong t2 = UWORD(0);                                         \
        ulong t1 = UWORD(0);                                         \
        ulong t0 = UWORD(0);                                         \
                                                                     \
        i = 0;                                                       \
        /* we can accumulate 8 terms if n == mod.n is such that */   \
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */   \
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */   \
        if (mod.n <= 6521908912666391107L)                           \
        {                                                            \
            slong u0, u1;                                            \
            for (; i+7 < len; )                                      \
            {                                                        \
                umul_ppmm(u1, u0, (expr1), (expr2));                 \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
                i++;                                                 \
                add_sssaaaaaa(t2, t1, t0,                            \
                              t2, t1, t0,                            \
                              UWORD(0), u1, u0);                     \
            }                                                        \
                                                                     \
            u0 = UWORD(0);                                           \
            u1 = UWORD(0);                                           \
            for (; i < len; i++)                                     \
            {                                                        \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);                  \
            }                                                        \
                                                                     \
            add_sssaaaaaa(t2, t1, t0,                                \
                          t2, t1, t0,                                \
                          UWORD(0), u1, u0);                         \
        }                                                            \
        else                                                         \
        {                                                            \
            for (; i < len; i++)                                     \
            {                                                        \
                umul_ppmm(s1, s0, (expr1), (expr2));                 \
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, s1, s0);    \
            }                                                        \
        }                                                            \
                                                                     \
        NMOD_RED(t2, t2, mod);                                       \
        NMOD_RED3(res, t2, t1, t0, mod);                             \
    }                                                                \
} while(0);

ulong
_nmod_vec_flintdot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs)
{
    ulong res;
    slong i;
    NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i], mod, nlimbs);
    return res;
}

// matmul, C does not alias A or B
void nmod_mat_mul_flint(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // now let's compute
    const int nlimbs = _nmod_vec_dot_bound_limbs(A->c, A->mod);
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
            C->rows[i][j] = _nmod_vec_flintdot(A->rows[i], BT->rows[j], A->c, A->mod, nlimbs);

    nmod_mat_clear(BT);
}

ulong
_nmod_vec_newdot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs)
{
    ulong res;
    slong i;
    _NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i], mod, nlimbs);
    return res;
}

void nmod_mat_mul_newdot(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // now let's compute
    const int nlimbs = _nmod_vec_dot_bound_limbs(A->c, A->mod);
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
            C->rows[i][j] = _nmod_vec_newdot(A->rows[i], BT->rows[j], A->c, A->mod, nlimbs);

    nmod_mat_clear(BT);
}

/*--------------------------------------------------------*/
/* timing of dot-expr against current flint               */
/*--------------------------------------------------------*/

ulong time_dotexpr_vs_flint_plain(ulong len, ulong n, flint_rand_t state)
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

    nn_srcptr v1i, v2i;

    { // TEST
        ulong res_new;
        ulong j;
        _NMOD_VEC_DOT(res_new, j, len, v1[j], v2[j], mod, n_limbs);
        ulong res_flint;
        v1i = v1; v2i = v2;
        NMOD_VEC_DOT(res_flint, j, len, v1i[j], v2i[j], mod, n_limbs);
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
        ulong j;
        for (slong i = 0; i < NB_ITER; i++) // warmup
        {
            _NMOD_VEC_DOT(buf, j, len, v1[j], v2[j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            _NMOD_VEC_DOT(buf, j, len, v1[j], v2[j], mod, n_limbs);
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
        ulong j;
        for (slong i = 0; i < NB_ITER; i++) // warmup
        {
            v1i = v1; v2i = v2;
            NMOD_VEC_DOT(buf, j, len, v1i[j], v2i[j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            v1i = v1; v2i = v2;
            NMOD_VEC_DOT(buf, j, len, v1i[j], v2i[j], mod, n_limbs);
            res += buf;
        }
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t2 /= nb_iter;

    printf("%.1e\t%.1e\t%.1e\t", t1, t2, t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_dotexpr_vs_flint_rev2(ulong len, ulong n, flint_rand_t state)
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

    nn_srcptr v1i, v2i;

    { // TEST
        ulong res_new;
        ulong j;
        _NMOD_VEC_DOT(res_new, j, len/2, v1[len - 1 - 2*j], v2[len - 1 - 2*j], mod, n_limbs);
        ulong res_flint;
        v1i = v1; v2i = v2;
        NMOD_VEC_DOT(res_flint, j, len/2, v1i[len - 1 - 2*j], v2i[len - 1 - 2*j], mod, n_limbs);
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
        ulong j;
        for (slong i = 0; i < NB_ITER; i++) // warmup
        {
            _NMOD_VEC_DOT(buf, j, len/2, v1[len - 1 - 2*j], v2[len - 1 - 2*j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            _NMOD_VEC_DOT(buf, j, len/2, v1[len - 1 - 2*j], v2[len - 1 - 2*j], mod, n_limbs);
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
        ulong j;
        for (slong i = 0; i < NB_ITER; i++) // warmup
        {
            v1i = v1; v2i = v2;
            NMOD_VEC_DOT(buf, j, len/2, v1i[len - 1 - 2*j], v2i[len - 1 - 2*j], mod, n_limbs);
            res += buf;
        }

        tt = clock();
        for (slong i = 0; i < NB_ITER; i++)
        {
            v1i = v1; v2i = v2;
            NMOD_VEC_DOT(buf, j, len/2, v1i[len - 1 - 2*j], v2i[len - 1 - 2*j], mod, n_limbs);
            res += buf;
        }
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += NB_ITER;
    }
    t2 /= nb_iter;

    printf("%.1e\t%.1e\t%.1e\t", t1, t2, t2/t1);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return res;
}

ulong time_dotexpr_vs_flint_matmul(ulong len, ulong n, flint_rand_t state)
{
    if (len > 1000)
    {
        printf("inf\tinf\t1.0\t");
        return 0;
    }
    else
    {
        nmod_mat_t mat1, mat2, mat, acc;
        nmod_mat_init(mat1, len, len, n);
        nmod_mat_rand(mat1, state);
        nmod_mat_init(mat2, len, len, n);
        nmod_mat_rand(mat2, state);
        nmod_mat_init(mat, len, len, n);
        nmod_mat_init(acc, len, len, n);

        { // TEST
            nmod_mat_mul_flint(mat, mat1, mat2);
            nmod_mat_mul_newdot(acc, mat1, mat2);
            if (!nmod_mat_equal(mat, acc))
            {
                printf("\nMATMUL ERROR!\n");
                return 0;
            }
        }

        double t1, t2;
        clock_t tt;
        long nb_iter;

        t1 = 0.0; nb_iter = 0;
        while (t1 < TIME_THRES)
        {
            for (slong i = 0; i < NB_MAT_ITER; i++) // warmup
            {
                nmod_mat_mul_newdot(mat, mat1, mat2);
                nmod_mat_add(acc, acc, mat);
            }

            tt = clock();
            for (slong i = 0; i < NB_MAT_ITER; i++)
            {
                nmod_mat_mul_newdot(mat, mat1, mat2);
                nmod_mat_add(acc, acc, mat);
            }
            t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += NB_MAT_ITER;
        }
        t1 /= nb_iter;

        t2 = 0.0; nb_iter = 0;
        while (t2 < TIME_THRES)
        {
            for (slong i = 0; i < NB_MAT_ITER; i++) // warmup
            {
                nmod_mat_mul_flint(mat, mat1, mat2);
                nmod_mat_add(acc, acc, mat);
            }

            tt = clock();
            for (slong i = 0; i < NB_MAT_ITER; i++)
            {
                nmod_mat_mul_flint(mat, mat1, mat2);
                nmod_mat_add(acc, acc, mat);
            }
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += NB_MAT_ITER;
        }
        t2 /= nb_iter;

        printf("%.1e\t%.1e\t%.1e\t", t1, t2, t2/t1);

        nmod_mat_clear(mat1);
        nmod_mat_clear(mat2);
        nmod_mat_clear(mat);
        nmod_mat_clear(acc);
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
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125);

#if SMALL_SUITE
    const slong nlens = 4;
    const ulong lens[] = {2, 20, 200, 2000};

    const slong nbits = 10;
    const ulong bits[] = {20, 28, 30, 32, 33, 50, 60, 62, 63, 64};
#else
    const slong nlens = 12;
    const ulong lens[] = {2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000};

    const slong nbits = 19;
    const ulong bits[] = {17, 20, 23, 26, 29, 30, 31, 32, 33, 40, 50, 55, 57, 59, 60, 61, 62, 63, 64};
#endif

    const slong nfuns = 3;
    typedef ulong (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_dotexpr_vs_flint_plain,                // 0
        time_dotexpr_vs_flint_rev2,                 // 1
        time_dotexpr_vs_flint_matmul,                 // 2
        //time_dot_dotrev_vs_flint_cf,                // ??
        //time_dot_dotrev_vs_flint_cu,                // ??
    };

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
