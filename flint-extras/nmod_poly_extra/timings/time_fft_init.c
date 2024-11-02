#include "flint/nmod.h"
#include "flint/profiler.h"

#include "nmod_poly_fft.h"
#include <flint/flint.h>

#define num_primes 5

typedef struct
{
   ulong prime;
   ulong depth;
   ulong maxdepth;
} info_t;

void sample_old_init_set(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong depth = info->depth;
    const ulong maxdepth = info->maxdepth;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    // modulus, roots of unity
    nmod_t mod;
    nmod_init(&mod, p);
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxdepth, mod);
    ulong w = nmod_pow_ui(w0, 1UL<<(maxdepth - depth), mod);

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
        {
            n_fft_old_ctx_t F;
            n_fft_old_ctx_init_set(F, w, depth, p);
            n_fft_old_ctx_clear(F);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

void sample_init2_root(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong p = info->prime;
    ulong depth = info->depth;
    ulong maxdepth = info->maxdepth;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    // modulus, roots of unity
    nmod_t mod;
    nmod_init(&mod, p);
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxdepth, mod);
    ulong w = nmod_pow_ui(w0, 1UL<<(maxdepth - depth), mod);

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
        {
            n_fft_ctx_t F;
            n_fft_ctx_init2_root(F, w, depth, depth, p);
            n_fft_ctx_clear(F);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/*-----------------------------------------------------------------*/
/* initialize context for FFT for several bit lengths and depths   */
/*-----------------------------------------------------------------*/
void time_fft_init(ulong * primes, ulong * max_depths)
{
    for (ulong k = 4; k < num_primes; k++)
    {
        for (ulong depth = 3; depth <= max_depths[k]; depth++)
        {
            printf("%ld\t", depth);

            info_t info;
            info.prime = primes[k];
            info.maxdepth = max_depths[k];
            info.depth = depth;

            const ulong len = UWORD(1) << depth;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            double min[10];
            double max;

            prof_repeat(min+0, &max, sample_init2_root, (void *) &info);
            prof_repeat(min+1, &max, sample_old_init_set, (void *) &info);

            flint_printf("\t%.1e|%.1e\t||\t%.1e|%.1e",
                    min[0]/(double)FLINT_CLOCK_SCALE_FACTOR/len/rep,
                    min[1]/(double)FLINT_CLOCK_SCALE_FACTOR/len/rep,
                    min[0]/(double)1000000/rep,
                    min[1]/(double)1000000/rep
                    );
            flint_printf("\n");
        }
    }

}

/*------------------------------------------------------------*/
/* main just calls time_init_set()                            */
/*------------------------------------------------------------*/
int main()
{
    printf("- depth is log(fft length)\n");
    printf("- timing init FFT context at this depth\n");
    printf("depth\t\tred init new\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_depths[num_primes] = { 18, 25, 25, 25, 25 };

    time_fft_init(primes, max_depths);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
