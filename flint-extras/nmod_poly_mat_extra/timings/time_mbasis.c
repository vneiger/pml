#include <stdlib.h>
#include <time.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_utils.h" // for rand
#include "nmod_poly_mat_approximant.h"

#define NUMBER_MBASIS 5

// testing different variants of mbasis implementation
static void (*m_basis[NUMBER_MBASIS])(nmod_poly_mat_t, long * ,
                                       const nmod_poly_mat_t, ulong, const long*) =
{
    mbasis, mbasisII, mbasisIII, mbasisIV, mbasisV
    //mbasisIV
};

static char* nb_mbasis[NUMBER_MBASIS] =  {"I", "II", "III", "IV", "V"};
//static char* nb_mbasis[NUMBER_MBASIS] =  {"IV"};

// main function
void benchmark_mbasis(slong rdim, slong cdim, slong sigma, slong len,
                       ulong prime, flint_rand_t state)
{
    // create random matrices
    nmod_poly_mat_t mat;

    nmod_poly_mat_init(mat, rdim, cdim, prime);
    nmod_poly_mat_rand(mat, state, len + 1);

    nmod_poly_mat_t res_mbasis;
    slong shifts[rdim], res_shifts[rdim];

    for (slong i = 0; i < rdim; i++)
        shifts[i] = 0;

    double thres = 0.2;

    // parameters for measuring time
    double t = 0.0;
    clock_t tt;
    long nb_iter = 0;

    // let's go
    for (int i = 0; i < NUMBER_MBASIS; i++)
    {
        t = 0.0;
        nb_iter = 0;

        while (t<thres)
        {
            nmod_poly_mat_init(res_mbasis, rdim, rdim, prime);
            tt = clock();
            m_basis[i](res_mbasis, res_shifts, mat, sigma, shifts);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
            nmod_poly_mat_clear(res_mbasis);
        }
        t /= nb_iter;
        printf("%s\t%ld\t%ld\t%ld\t%ld\t%f\n", nb_mbasis[i], rdim, cdim, sigma, len, t);
    }

    {
        t = 0.0;
        nb_iter = 0;

        while (t<thres)
        {
            nmod_poly_mat_init(res_mbasis, rdim, rdim, prime);
            tt = clock();
            slong * cshift = flint_calloc(mat->r, sizeof(slong));
            nmod_poly_mat_mbasis(res_mbasis, cshift, mat, sigma);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
            nmod_poly_mat_clear(res_mbasis);
        }
        t /= nb_iter;
        printf("%s\t%ld\t%ld\t%ld\t%ld\t%f\n", "NEW", rdim, cdim, sigma, len, t);
    }
    printf("\n");

    nmod_poly_mat_clear(mat);
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
    flint_rand_init(state);
    const ulong prime = n_randprime(state, nbits, 0);

    slong rdims[] = { 2, 4, 8, 16, 32, 64, 128, 256, 512};
    slong orders[] = { 4, 8, 16, 32, 64};

    printf("Bench mbasis\n");
    printf("nbits=%ld, prime=%ld\n", nbits, prime);
    printf("n\trdim\tcdim\torder\tdeg\ttime\n");
    printf("~~~~WARMUP~~~~\n");
    benchmark_mbasis(4, 2, 128, 128, 997, state);
    printf("~~~~WARMUP DONE~~~~\n");
    for (size_t i = 0; i < sizeof(rdims) / sizeof(rdims[0]); ++i)
    {
        long rdim = rdims[i];
        for (size_t j = 0; j < sizeof(orders) / sizeof(orders[0]); ++j)
        {
            long ord = orders[j];
            benchmark_mbasis(rdim, 1, ord, ord, prime, state);
            benchmark_mbasis(1, rdim, ord, ord, prime, state);
            benchmark_mbasis(rdim, rdim / 2, ord, ord, prime, state);
            benchmark_mbasis(rdim / 2, rdim, ord, ord, prime, state);
            benchmark_mbasis(rdim, rdim-1, ord, ord, prime, state);
            benchmark_mbasis(rdim-1, rdim, ord, ord, prime, state);
        }
    }

    flint_rand_clear(state);
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

    printf("Bench mbasis:\n");
    printf("nbits=%ld, prime=%ld, rdim=%ld, cdim=%ld, deg=%ld\n",nbits, prime, rdim, cdim, deg);
    printf("~~~~WARMUP~~~~\n");
    benchmark_mbasis(4, 2, 128, 128, 997, state);
    printf("~~~~WARMUP DONE~~~~\n");
    printf("n\trdim\tcdim\torder\tdeg\ttime\n");
    benchmark_mbasis(rdim, cdim, deg, deg, prime, state);
}


int main(int argc, char *argv[])
{
    setlinebuf(stdout);

    srand(time(NULL));
    flint_rand_t state;
    flint_rand_init(state);
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

    flint_rand_clear(state);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
