#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include <stdlib.h>

static void warmup(flint_rand_t state)
{
    nmod_mat_t mat;
    nmod_mat_init(mat, 2000, 2000, 65537);
    nmod_mat_randfull(mat, state);
    nmod_mat_det(mat);
    nmod_mat_clear(mat);
}

// main function
void benchmark_generic_hnf_withtrans(slong rdim, slong cdim, slong len,
                                     ulong prime, flint_rand_t state)
{
    // create random matrix
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, rdim, cdim, prime);
    nmod_poly_mat_rand(mat, state, len);

    // init copy of mat
    nmod_poly_mat_t hnf;
    nmod_poly_mat_init(hnf, mat->r, mat->c, mat->modulus);

    // init unimodular transformation tsf such that hnf = tsf * mat
    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, mat->r, mat->r, mat->modulus);

    // initialize pivind list and rank profiles
    slong max_rank = FLINT_MIN(mat->r, mat->c);
    slong * pivind = flint_malloc(max_rank * sizeof(slong));

    // for timing
    long thres = 1000; // 1000ms = 1s
    long nb_iter = 0;
    long t_uref = 0; // time for uref, in ms
    long t_norm = 0; // time for normalization, in ms
    timeit_t timer;

    flint_printf("%ld\t%ld\t%ld\t", rdim, cdim, len);

    { // maxdeg, atomic
        nb_iter = 0; t_uref = 0; t_norm = 0;
        while (t_uref+t_norm < thres)
        {
            nmod_poly_mat_set(hnf, mat);
            nmod_poly_mat_one(tsf);
            timeit_start(timer);
            slong rk = nmod_poly_mat_uref_maxdeg_atomic(hnf, tsf, pivind);
            timeit_stop(timer);
            t_uref += timer->wall;
            timeit_start(timer);
            _normalize_uref(hnf, tsf, pivind, rk);
            timeit_stop(timer);
            t_norm += timer->wall;
            nb_iter++;
        }
        flint_printf("%0.1e\t%0.1e\t%0.1e\t",
                     (double)t_uref/nb_iter,
                     (double)t_norm/nb_iter,
                     (double)(t_uref+t_norm) / nb_iter);
    }

    { // revlex, xgcd
        nb_iter = 0; t_uref = 0; t_norm = 0;
        while (t_uref+t_norm < thres)
        {
            nmod_poly_mat_set(hnf, mat);
            nmod_poly_mat_one(tsf);
            timeit_start(timer);
            slong rk = nmod_poly_mat_uref_revlex_xgcd(hnf, tsf, pivind, NULL);
            timeit_stop(timer);
            t_uref += timer->wall;
            timeit_start(timer);
            _normalize_uref(hnf, tsf, pivind, rk);
            timeit_stop(timer);
            t_norm += timer->wall;
            nb_iter++;
        }
        flint_printf("%0.1e\t%0.1e\t%0.1e\t",
                     (double)t_uref/nb_iter,
                     (double)t_norm/nb_iter,
                     (double)(t_uref+t_norm) / nb_iter);
    }

    { // lex, xgcd
        nb_iter = 0; t_uref = 0; t_norm = 0;
        while (t_uref+t_norm < thres)
        {
            nmod_poly_mat_set(hnf, mat);
            nmod_poly_mat_one(tsf);
            timeit_start(timer);
            slong rk = nmod_poly_mat_uref_lex_xgcd(hnf, tsf, pivind, NULL);
            timeit_stop(timer);
            t_uref += timer->wall;
            timeit_start(timer);
            _normalize_uref(hnf, tsf, pivind, rk);
            timeit_stop(timer);
            t_norm += timer->wall;
            nb_iter++;
        }
        flint_printf("%0.1e\t%0.1e\t%0.1e\t",
                     (double)t_uref/nb_iter,
                     (double)t_norm/nb_iter,
                     (double)(t_uref+t_norm) / nb_iter);
    }

    { // Kannan-Bachem's algorithm
        nb_iter = 0; t_uref = 0;
        while (t_uref < thres)
        {
            nmod_poly_mat_set(hnf, mat);
            nmod_poly_mat_one(tsf);
            timeit_start(timer);
            nmod_poly_mat_hnf_ur_revlex_xgcd_delayed_zero(hnf, tsf, pivind, NULL);
            timeit_stop(timer);
            t_uref += timer->wall;
            nb_iter++;
        }
        flint_printf("%0.1e\t", (double)t_uref/nb_iter);
    }

    { // Mulder-Storjohann's algorithm
        nb_iter = 0; t_uref = 0; t_norm = 0;
        slong rk = 0;
        while (t_uref+t_norm < thres)
        {
            nmod_poly_mat_set(hnf, mat);
            nmod_poly_mat_one(tsf);
            timeit_start(timer);
            rk = nmod_poly_mat_uref_matrixgcd_iter(hnf, tsf, pivind, NULL, NULL);
            timeit_stop(timer);
            t_uref += timer->wall;
            if (rk > 0 || nmod_poly_mat_is_zero(mat)) // avoid non-generic early exit
            {
                timeit_start(timer);
                _normalize_uref(hnf, tsf, pivind, rk);
                timeit_stop(timer);
                t_norm += timer->wall;
            }
            nb_iter++;
        }
        if (rk <= 0 && !nmod_poly_mat_is_zero(mat)) // early exit, non-generic
            flint_printf("%0.1e\tXXX\tXXX\t", (double)t_uref/nb_iter);
        else
            flint_printf("%0.1e\t%0.1e\t%0.1e\t",
                        (double)t_uref/nb_iter,
                        (double)t_norm/nb_iter,
                        (double)(t_uref+t_norm) / nb_iter);
    }

    flint_printf("\n");
    nmod_poly_mat_clear(hnf);
    nmod_poly_mat_clear(tsf);
    flint_free(pivind);
}

/** Launches a series of benchmarks for a given prime size.
 *
 *  Launches benchmark for many matrix dimensions and lengths for a given
 *  bitlength for the prime defining the field of coefficients.
 *
 * \param nbits bitlength of the prime modulus
 * \param state flint's random generator
 * \return void
 */
void benchmark_nbits(ulong nbits, flint_rand_t state)
{
    flint_randinit(state);
    const ulong prime = n_randprime(state, nbits, 0);

    slong rdims[] = { 2, 4, 8, 16, 32, 64, 128 };
    slong lens[] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};

    printf("------------------------\n");
    printf("Bench HNF (times in ms)\n");
    printf("nbits=%ld, prime=%ld\n", nbits, prime);
    if (nbits < 15)
        printf("Warning: for 'small' primes, risk of non-generic CRP; Mulders-Storjohann's HNF may early exit before the computation is complete\n");
    printf("------------------------\n");
    printf("~~~~WARMUP~~~~\n");
    warmup(state);
    printf("~~~~WARMUP DONE~~~~\n");
    printf("------------------------\n");
    printf("rdim\tcdim\tlen\t");
    printf("Ros-ef\tRos-nm\tRos-tot\t");
    printf("rlx-ef\trlx-nm\trlx-tot\t");
    printf("lx-ef\tlx-nm\tlx-tot\t");
    printf("KB-tot\t");
    printf("mx-ef\tmx-nm\tmx-tot\t");
    printf("\n");
    for (size_t i = 0; i < sizeof(rdims) / sizeof(rdims[0]); ++i)
    {
        long rdim = rdims[i];
        for (size_t j = 0; j < sizeof(lens) / sizeof(lens[0]); ++j)
        {
            long len = lens[j];
            benchmark_generic_hnf_withtrans(rdim, 1, len, prime, state);
            benchmark_generic_hnf_withtrans(1, rdim, len, prime, state);
            benchmark_generic_hnf_withtrans(rdim, rdim/2, len, prime, state);
            benchmark_generic_hnf_withtrans(rdim/2, rdim, len, prime, state);
            benchmark_generic_hnf_withtrans(rdim, rdim, len, prime, state);
        }
    }

    flint_randclear(state);
}

/** Launches one benchmark for fixed parameters.
 *
 *  Launches benchmark for given dimensions, lengths, and modulus bitlength.
 *
 * \param nbits bitlength of the prime modulus
 * \param rdim row dimension of input matrix
 * \param cdim column dimension of input matrix
 * \param len length of polynomials in input matrix
 * \param state flint's random generator
 * \return void
 */
void benchmark_nbits_dim_deg(ulong nbits, ulong rdim, ulong cdim, ulong len, flint_rand_t state)
{
    const ulong prime = n_randprime(state, nbits, 0);

    printf("Bench HNF (times in ms):\n");
    printf("nbits=%ld, prime=%ld, rdim=%ld, cdim=%ld, len=%ld\n",nbits, prime, rdim, cdim, len);
    if (nbits < 15)
        printf("Warning: for 'small' primes, risk of non-generic CRP; Mulders-Storjohann's HNF may early exit before the computation is complete\n");
    printf("~~~~WARMUP~~~~\n");
    warmup(state);
    printf("~~~~WARMUP DONE~~~~\n");
    printf("rdim\tcdim\tlen\t");
    printf("Ros-ef\tRos-nm\tRos-tot\t");
    printf("rlx-ef\trlx-nm\trlx-tot\t");
    printf("lx-ef\tlx-nm\tlx-tot\t");
    printf("KB-tot\t");
    printf("mx-ef\tmx-nm\tmx-tot\t");
    printf("\n");
    benchmark_generic_hnf_withtrans(rdim, cdim, len, prime, state);
}

int main(int argc, char *argv[])
{
    setlinebuf(stdout);
    srand(time(NULL));
    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, rand(), rand());

    if (argc!=2 && argc!=5)
    {
        printf("Usage: %s nbits OR %s nbits rdim cdim len\n",argv[0],argv[0]);
        return 1;
    }

    if (argc==2)
    {
        warmup(state);
        ulong nbits = atoi(argv[1]);
        benchmark_nbits(nbits, state);
        return 0;
    }

    if (argc==5)
    {
        warmup(state);
        ulong nbits = atoi(argv[1]);
        ulong rdim = atoi(argv[2]);
        ulong cdim = atoi(argv[3]);
        ulong deg = atoi(argv[4]);
        benchmark_nbits_dim_deg(nbits,rdim,cdim,deg,state);
    }

    flint_randclear(state);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
