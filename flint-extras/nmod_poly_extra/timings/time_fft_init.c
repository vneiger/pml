#include "flint/nmod.h"
#include "flint/profiler.h"

#include "nmod_poly_fft.h"
#include <flint/flint.h>

#define num_primes 5

typedef struct
{
   ulong prime;
   ulong order;
   ulong maxorder;
} info_t;

void sample_init_set(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    const ulong p = info->prime;
    const ulong order = info->order;
    const ulong maxorder = info->maxorder;

    const ulong len = UWORD(1) << order;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    // modulus, roots of unity
    nmod_t mod;
    nmod_init(&mod, p);
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxorder, mod);
    ulong w = nmod_pow_ui(w0, 1UL<<(maxorder - order), mod);

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
        {
            nmod_fft_ctx_t F;
            nmod_fft_ctx_init_set(F, w, order, p);
            nmod_fft_ctx_clear(F);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}


void sample_init_set_red(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong p = info->prime;
    ulong order = info->order;
    ulong maxorder = info->maxorder;

    const ulong len = UWORD(1) << order;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    // modulus, roots of unity
    nmod_t mod;
    nmod_init(&mod, p);
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxorder, mod);
    ulong w = nmod_pow_ui(w0, 1UL<<(maxorder - order), mod);

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
        {
            nmod_fft_ctx_t F;
            nmod_fft_ctx_init_set_red(F, w, order, p);
            nmod_fft_ctx_clear_red(F);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

void sample_init_set_new(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong p = info->prime;
    ulong order = info->order;
    ulong maxorder = info->maxorder;

    const ulong len = UWORD(1) << order;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    // modulus, roots of unity
    nmod_t mod;
    nmod_init(&mod, p);
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> maxorder, mod);
    ulong w = nmod_pow_ui(w0, 1UL<<(maxorder - order), mod);

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
        {
            nmod_fft_ctx_t F;
            nmod_fft_ctx_init_set_new(F, w, order, p);
            nmod_fft_ctx_clear_new(F);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/*-----------------------------------------------------------------*/
/* initialize context for FFT for several bit lengths and orders   */
/*-----------------------------------------------------------------*/
void time_fft_init(ulong * primes, ulong * max_orders)
{
    for (ulong k = 4; k < num_primes; k++)
    {
        for (ulong order = 3; order <= max_orders[k]; order++)
        {
            printf("%ld\t", order);

            info_t info;
            info.prime = primes[k];
            info.maxorder = max_orders[k];
            info.order = order;

            const ulong len = UWORD(1) << order;
            const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

            double min[10];
            double max;

            prof_repeat(min+0, &max, sample_init_set_red, (void *) &info);
            prof_repeat(min+1, &max, sample_init_set, (void *) &info);
            prof_repeat(min+2, &max, sample_init_set_new, (void *) &info);

            flint_printf("\t%.1e|%.1e|%.1e",
                    min[0]/(double)FLINT_CLOCK_SCALE_FACTOR/len/rep,
                    min[1]/(double)FLINT_CLOCK_SCALE_FACTOR/len/rep,
                    min[2]/(double)FLINT_CLOCK_SCALE_FACTOR/len/rep 
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
    printf("- order is log(fft length)\n");
    printf("- timing init FFT context at this order\n");
    printf("order\tinit     initred  new\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_orders[num_primes] = { 18, 25, 25, 25, 25 };

    time_fft_init(primes, max_orders);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
