#include <flint/nmod_poly_mat.h>
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_approximant.h"
#include "time.h"

// main function
void benchmark_pmbasis(slong rdim, slong cdim, slong sigma, slong len,
                       ulong prime, flint_rand_t state)
{
    // create random matrices
    nmod_poly_mat_t mat;

    nmod_poly_mat_init(mat, rdim, cdim, prime);
    nmod_poly_mat_randtest(mat, state, len + 1);

    nmod_poly_mat_t res_mbasis;
    slong shifts[rdim], res_shifts[rdim];

    for (slong i = 0; i < rdim; i++)
        shifts[i] = 0;

    nmod_poly_mat_init(res_mbasis, rdim, rdim, prime);

    double thres = 1.0;

    // parameters for measuring time
    double t = 0.0;
    clock_t tt;
    long nb_iter = 0;

    // let's go
    t = 0.0;
    nb_iter = 0;

    {
        t = 0.0;
        nb_iter = 0;

        while (t<thres)
        {
            nmod_poly_mat_zero(res_mbasis);
            tt = clock();
            slong * cshift = flint_calloc(mat->r, sizeof(slong));
            nmod_poly_mat_pmbasis(res_mbasis, cshift, mat, sigma);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t /= nb_iter;
        printf("%s\t%ld\t%ld\t%ld\t%ld\t%f\n", "NEW", rdim, cdim, sigma, len, t);
    }
}

/** Launches a series of benchmarks for a given prime size.
 *
 *  Launches benchmark for many matrix dimensions and orders for a given
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
    slong sigmas[] = { 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};

    printf("Bench pmbasis\n");
    printf("nbits=%ld, prime=%ld\n", nbits, prime);
    printf("n\trdim\tcdim\torder\tdeg\ttime\n");
    printf("~~~~WARMUP~~~~\n");
    benchmark_pmbasis(4, 2, 128, 128, 997, state);
    printf("~~~~WARMUP DONE~~~~\n");
    for (size_t i = 0; i < sizeof(rdims) / sizeof(rdims[0]); ++i)
    {
        long rdim = rdims[i];

        for (size_t j = 0; j < sizeof(sigmas) / sizeof(sigmas[0]); ++j)
        {
            long ord = sigmas[j];
            benchmark_pmbasis(rdim, rdim / 2, ord, ord, prime, state);
        }
    }

    flint_randclear(state);
}

/** Launches one benchmark for fixed parameters.
 *
 *  Launches benchmark for given dimensions, orders, and modulus bitlength.
 *
 * \param nbits bitlength of the prime modulus
 * \param rdim row dimension of input matrix
 * \param cdim column dimension of input matrix
 * \param state flint's random generator
 * \return void
 */
void benchmark_nbits_dim_deg(ulong nbits, ulong rdim, ulong cdim, ulong deg, flint_rand_t state)
{
    const ulong prime = n_randprime(state, nbits, 0);

    printf("Bench pmbasis:\n");
    printf("nbits=%ld, prime=%ld, rdim=%ld, cdim=%ld, deg=%ld\n",nbits, prime, rdim, cdim, deg);
    printf("~~~~WARMUP~~~~\n");
    benchmark_pmbasis(4, 2, 128, 128, 997, state);
    printf("~~~~WARMUP DONE~~~~\n");
    printf("n\trdim\tcdim\torder\tdeg\ttime\n");
    benchmark_pmbasis(rdim, cdim, deg, deg, prime, state);
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
        printf("Usage: %s nbits OR %s nbits rdim cdim order\n",argv[0],argv[0]);
        return 1;
    }

    if (argc==2)
    {
        ulong nbits = atoi(argv[1]);
        benchmark_nbits(nbits, state);
        return 0;
    }

    if (argc==5)
    {
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
