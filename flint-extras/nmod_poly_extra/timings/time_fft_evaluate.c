#include <flint/profiler.h>
#include <flint/nmod_vec.h>
#include <flint/fft_small.h>
#include "nmod_poly_fft.h"

#define VERSIONS 1

#define num_primes 5
#define num_strides 1

typedef struct
{
   ulong prime;
   ulong depth;
   ulong maxdepth;
   ulong stride;
} info_t;

#define SAMPLE(fun, _variant)                                                    \
void sample_##fun##_variant(void * arg, ulong count)                             \
{                                                                                \
    info_t * info = (info_t *) arg;                                              \
    const ulong p = info->prime;                                                 \
    const ulong depth = info->depth;                                             \
    const ulong maxdepth = info->maxdepth;                                       \
    const ulong stride = info->stride;                                           \
                                                                                 \
    const ulong len = stride * (UWORD(1) << depth);                              \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));                \
                                                                                 \
    /* modulus, roots of unity */                                                \
    nmod_t mod;                                                                  \
    nmod_init(&mod, p);                                                          \
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxdepth, mod); \
    ulong w = nmod_pow_ui(w0, 1UL<<(maxdepth - depth), mod);                     \
    n_fft_ctx_t F;                                                               \
    n_fft_ctx_init2_root(F, w, depth, depth, p);                                 \
                                                                                 \
    FLINT_TEST_INIT(state);                                                      \
                                                                                 \
    ulong * coeffs = _nmod_vec_init(len);                                        \
    _nmod_vec_randtest(coeffs, state, len, mod);                                 \
                                                                                 \
    for (ulong i = 0; i < count; i++)                                            \
    {                                                                            \
        prof_start();                                                            \
        for (ulong j = 0; j < rep; j++)                                          \
            _n_fft_##fun##_variant(coeffs, len, depth, F);                       \
        prof_stop();                                                             \
    }                                                                            \
                                                                                 \
    n_fft_ctx_clear(F);                                                          \
    FLINT_TEST_CLEAR(state);                                                     \
}                                                                                \

#define SAMPLEOLD(fun)                                                           \
void sample_##fun(void * arg, ulong count)                                       \
{                                                                                \
    info_t * info = (info_t *) arg;                                              \
    const ulong p = info->prime;                                                 \
    const ulong depth = info->depth;                                             \
    const ulong maxdepth = info->maxdepth;                                       \
                                                                                 \
    const ulong len = UWORD(1) << depth;                                         \
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));                \
                                                                                 \
    /* modulus, roots of unity */                                                \
    nmod_t mod;                                                                  \
    nmod_init(&mod, p);                                                          \
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxdepth, mod); \
    ulong w = nmod_pow_ui(w0, 1UL<<(maxdepth - depth), mod);                     \
    n_fft_old_ctx_t F;                                                           \
    n_fft_old_ctx_init_set(F, w, depth, p);                                      \
                                                                                 \
    FLINT_TEST_INIT(state);                                                      \
                                                                                 \
    ulong * coeffs = _nmod_vec_init(len);                                        \
    _nmod_vec_randtest(coeffs, state, len, mod);                                 \
                                                                                 \
    for (ulong i = 0; i < count; i++)                                            \
    {                                                                            \
        prof_start();                                                            \
        for (ulong j = 0; j < rep; j++)                                          \
            _n_fft_##fun(coeffs, len, depth, F);                                 \
        prof_stop();                                                             \
    }                                                                            \
                                                                                 \
    n_fft_old_ctx_clear(F);                                                      \
    FLINT_TEST_CLEAR(state);                                                     \
}                                                                                \


SAMPLEOLD(old_dif_rec2_lazy)
SAMPLEOLD(old_dif_iter2_lazy)
SAMPLEOLD(old_dif_rec4_lazy)
SAMPLEOLD(old_dif_rec8_lazy)
SAMPLE(red_rec2_lazy, )
SAMPLE(red_iter2_lazy, )
SAMPLE(red_rec4_lazy, )
SAMPLE(red_rec2_lazy, _stride)

void sample_sd_fft(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong depth = info->depth;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    sd_fft_ctx_t Q;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_fit_depth(Q, depth);

    ulong sz = sd_fft_ctx_data_size(depth)*sizeof(double);

    FLINT_TEST_INIT(state);

    nmod_t mod;
    nmod_init(&mod, p);
    ulong * coeffs = _nmod_vec_init(len);
    _nmod_vec_randtest(coeffs, state, len, mod);

    double* data = flint_aligned_alloc(4096, n_round_up(sz, 4096));
    for (ulong i = 0; i < len; i++)
        data[i] = coeffs[i];

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
            sd_fft_trunc(Q, data, depth, len, len);
        prof_stop();
    }

    sd_fft_ctx_clear(Q);
    FLINT_TEST_CLEAR(state);
}

/*------------------------------------------------------------*/
/* computes init for FFT for several bit lengths and depths   */
/*------------------------------------------------------------*/
void time_evaluate()
{
    flint_printf("- depth is log(fft length)\n");
    flint_printf("- timing init FFT tables + DIF evaluate for several bit lengths and depths\n");
    //flint_printf("depth\tsd_fft\trec2\trec4\titer2\t   ||\tsd_fft\trec2\trec4\titer2\toldr2\toldr4\toldr8\toldi2\t\n");
    flint_printf("depth\tsd_fft\trec2\trec4\titer2\tr2s1\tr2s2\tr2s4\tr2s9\tr2s16\tr2s100\t\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_depths[num_primes] = { 18, 25, 25, 25, 25 };

    ulong strides[num_strides] = {1, 2, 4, 9, 16, 100};
    ulong rep_stride[num_strides] = {0, 0, 0, 0, 0, 0};

    for (ulong k = 3; k < 4; k++)
    {
        for (ulong depth = 3; depth <= max_depths[k]; depth++)
        {
            printf("%ld\t", depth);

            info_t info;
            info.prime = primes[k];
            info.maxdepth = max_depths[k];
            info.depth = depth;
            info.stride = 1;

            const ulong len = UWORD(1) << depth;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            double min[15];
            double max;

            prof_repeat(min+0, &max, sample_sd_fft, (void *) &info);
            prof_repeat(min+1, &max, sample_red_rec2_lazy, (void *) &info);
            prof_repeat(min+2, &max, sample_red_rec4_lazy, (void *) &info);
            prof_repeat(min+3, &max, sample_red_iter2_lazy, (void *) &info);
            //prof_repeat(min+4, &max, sample_old_dif_rec2_lazy, (void *) &info);
            //prof_repeat(min+5, &max, sample_old_dif_rec4_lazy, (void *) &info);
            //prof_repeat(min+6, &max, sample_old_dif_iter2_lazy, (void *) &info);
            //prof_repeat(min+7, &max, sample_old_dif_rec8_lazy, (void *) &info);

            for (ulong i = 0; i < num_strides; i++)
            {
                info.stride = strides[i];
                rep_stride[i] = FLINT_MAX(1, FLINT_MIN(1000, 1000000 / (len * strides[i])));
                prof_repeat(min+8+i, &max, sample_red_rec2_lazy_stride, (void *) &info);
            }

            flint_printf("%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\n",
                    min[0]/(double)1000000/rep,
                    min[1]/(double)1000000/rep,
                    min[2]/(double)1000000/rep,
                    min[3]/(double)1000000/rep,
                    //min[4]/(double)1000000/rep,
                    //min[5]/(double)1000000/rep,
                    //min[6]/(double)1000000/rep,
                    //min[7]/(double)1000000/rep,
                    min[8]/(double)1000000/rep_stride[0],
                    min[9]/(double)1000000/rep_stride[1],
                    min[10]/(double)1000000/rep_stride[2],
                    min[11]/(double)1000000/rep_stride[3],
                    min[12]/(double)1000000/rep_stride[4],
                    min[13]/(double)1000000/rep_stride[5]
                    );
        }
    }
}


/*------------------------------------------------------------*/
/* main just calls time_init_set()                            */
/*------------------------------------------------------------*/
int main()
{
    time_evaluate();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
