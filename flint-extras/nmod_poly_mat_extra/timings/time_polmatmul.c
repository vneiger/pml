/** \file time_polmatmul.c
 *
 * This file benchmarks several polynomial matrix multiplication
 * implementations from FLINT and FLINT-based PML.
 *
 */

#include <time.h>
#include <stdlib.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"


/** Bench of polmatmul for given dimension, degree, modulus.
 *
 *  This function takes the parameters (matrix dimensions, degree bound,
 *  and prime modulus), creates two random matrices a and b, and measures
 *  the time it takes to compute c = a*b.
 *
 * \param rdim row dimension
 * \param idim inner dimension
 * \param cdim column dimension
 * \param deg degree bound (strict)
 * \param prime prime modulus
 * \param fftprime 0 if not an fft prime
 * \param state random state for generating random matrix entries
 * \return void
 */
void benchmark_polmatmul(long rdim, long idim, long cdim, long deg, ulong prime, int fftprime, flint_rand_t state)
{
    // create random matrices
    nmod_poly_mat_t a, b, c1, c2;

    flint_randinit(state);
    nmod_poly_mat_init(a, rdim, idim, prime);
    nmod_poly_mat_rand(a, state, deg);

    nmod_poly_mat_init(b, idim, cdim, prime);
    nmod_poly_mat_rand(b, state, deg);
    flint_randclear(state);

    printf("%ld\t%ld\t%ld\t%ld",rdim,idim,cdim,deg);

    // parameters for measuring time
    double t = 0.0;
    clock_t tt;
    long nb_iter = 0;

    // warmup
    t = 0.0;
    nb_iter = 0;
    while (t<0.5)
    {
            tt = clock();
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
    }

    if (1)
    { // FLINT native
        nmod_poly_mat_init(c1, rdim, cdim, prime);
        t = 0.0; nb_iter = 0;
        while (t<0.5)
        {
            tt = clock();
            nmod_poly_mat_mul(c1, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
            }
        t /= nb_iter;
        printf("\t%f",t);
    }
    else
        printf("\t-\t");

    if (1)
    { // PML 3-primes
        nmod_poly_mat_init(c2, rdim, cdim, prime);
        t = 0.0; nb_iter = 0;
        while (t<0.5)
        {
            tt = clock();
            nmod_poly_mat_mul_3_primes(c2, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
            t /= nb_iter;
            printf("\t%f", t);
            if (!nmod_poly_mat_equal(c1, c2))
            {
                printf("\terror with 3 primes\n");
                return;
            }
    }
    else
        printf("\t-\t");

    // should test if the field is large enough
    if (1)
    { // PML geometric
        nmod_poly_mat_init(c2, rdim, cdim, prime);
        t = 0.0; nb_iter = 0;
        while (t<0.5)
        {
            tt = clock();
            nmod_poly_mat_mul_geometric(c2, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t /= nb_iter;
        printf("\t%f", t);
        if (!nmod_poly_mat_equal(c1, c2))
        {
            printf("\terror with geometric\n");
            return;
        }
    }
    else
        printf("\t-\t");

    // should test if the field is large enough
    if (deg < 100)
    { // PML VDM1
        nmod_poly_mat_init(c2, rdim, cdim, prime);
        t = 0.0; nb_iter = 0;
        while (t<0.5)
        {
            tt = clock();
            nmod_poly_mat_mul_vdm1(c2, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t /= nb_iter;
        printf("\t%f", t);
        if (!nmod_poly_mat_equal(c1, c2))
        {
            printf("\terror with vdm1\n");
            return;
        }
    }
    else
        printf("\t-\t");

    // should test if the field is large enough
    if (deg < 100)
    { // PML VDM2
        nmod_poly_mat_init(c2, rdim, cdim, prime);
        t = 0.0; nb_iter = 0;
        while (t<0.5)
        {
            tt = clock();
            nmod_poly_mat_mul_vdm2(c2, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t /= nb_iter;
        printf("\t%f", t);
        if (!nmod_poly_mat_equal(c1, c2))
        {
            printf("\terror with vdm2\n");
            return;
        }
    }
    else
        printf("\t-\t");
        
    if (rdim < 20 && cdim < 20)
    { // PML WAKSMAN
        nmod_poly_mat_init(c2, rdim, cdim, prime);
        t = 0.0; nb_iter = 0;
        while (t<0.5)
        {
            tt = clock();
            nmod_poly_mat_mul_waksman(c2, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t /= nb_iter;
        printf("\t%f", t);
        if (!nmod_poly_mat_equal(c1, c2))
        {
            printf("\terror with waksman\n");
            return;
        }
    }
    else
        printf("\t-");
    

    nmod_poly_mat_clear(c1);
    nmod_poly_mat_clear(c2);
            
    /* if (fftprime) */
    /* { // PML TFT */
    /*     nmod_poly_mat_t c; */
    /*     nmod_poly_mat_init(c, rdim, cdim, prime); */
    /*     t = 0.0; nb_iter = 0; */
    /*     while (t<0.5) */
    /*     { */
    /*         tt = clock(); */
    /*         nmod_poly_mat_mul_tft(c, a, b); */
    /*         t += (double)(clock()-tt) / CLOCKS_PER_SEC; */
    /*         ++nb_iter; */
    /*     } */
    /*     t /= nb_iter; */
    /*     printf("\t%f",t); */
    /*     nmod_poly_mat_clear(c); */
    /* } */
    nmod_poly_mat_clear(a);
    nmod_poly_mat_clear(b);
    printf("\n");
}

/** Launches a series of benchmarks for a given prime size.
 *
 *  Launches benchmark for many matrix dimensions and degrees for a given
 *  bitlength for the prime defining the field of coefficients.
 *
 * \param nbits bitlength of the prime modulus
 * \return void
 */
void benchmark_nbits(ulong nbits)
{
    flint_rand_t state;
    flint_randinit(state);
    const ulong prime = n_randprime(state, nbits, 0);
    flint_randclear(state);

    //long dims[] = { 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, };
    long dims[] = { 2, 4, 10, 20, 32, 75, 180 };
    //long degs[] =
    //{
    //    8, 10, 12, 14,
    //    16, 20, 24, 28,
    //    32, 40, 48, 56,
    //    64, 80, 96, 112,
    //    128, 160, 192, 224,
    //    256, 320, 384, 448,
    //    512, 640, 768, 896,
    //    1024, 1280, 1536, 1792,
    //    2048, 2560, 3072, 3584,
    //    4096, 5120, 6144, 7168,
    //    8192, 10240, 12288, 14336,
    //    16384, 20480, 24576, 28672,
    //    32768, 40960, 49152, 57344
    //};
    long degs[] =
    {
        8,
        32,
        128,
        512,
        2048,
        8192,
    };

    printf("Bench square polynomial matrix multiplication\n");
    printf("nbits=%ld, prime=%ld\n",nbits,prime);
    printf("rdim\tidim\tcdim\tdeg\tflint\t\tpml\n");
    for (size_t i=0; i<sizeof(dims)/sizeof(dims[0]); ++i)
    {
        long dim = dims[i];
        //long maxdeg = 50000;
        //if (dim > 16)
        //    maxdeg = 10000;
        //if (dim > 128)
        //    maxdeg = 1600;
        //if (dim > 256)
        //    maxdeg = 800;
        for (size_t j=0; j<sizeof(degs)/sizeof(degs[0]); ++j)
        {
            long d = degs[j];
            //if (d < maxdeg)
            if (dim*dim*dim*d < 30000000000)
                benchmark_polmatmul(dim,dim,dim,d,prime,0,state);
        }
    }
}

/** Launches single benchmark for specified parameters.
 *
 *  Launches one benchmark for the given square matrix dimensions, matrix
 *  degree, and bitlength for the prime defining the field of coefficients.
 *
 * \param nbits bitlength of the prime modulus
 * \param dim square matrix dimension
 * \param deg matrix degree (strict bound)
 * \return void
 */
void benchmark_nbits_dim_deg(ulong nbits, ulong dim, ulong deg)
{
    flint_rand_t state;
    flint_randinit(state);
    const ulong prime = n_randprime(state, nbits, 0);
    flint_randclear(state);

    printf("Bench square polynomial matrix multiplication:\n");
    printf("nbits=%ld, dim=%ld, deg=%ld\n",nbits,dim,deg);
    printf("nbits=%ld, prime=%ld, dim=%ld, deg=%ld\n",nbits,prime,dim,deg);
    printf("rdim\tidim\tcdim\tdeg\tflint\t\tpml(3 primes)\tpml(geometric)\tpml(vdm1)\tpml(vdm2)\tpml(waksman)\n");
    benchmark_polmatmul(dim,dim,dim,deg,prime,0,state);
}

int main(int argc, char *argv[])
{
    if (argc!=2 && argc!=4)
    {
        printf("Usage: %s nbits OR %s nbits sz deg\n",argv[0],argv[0]);
        return 1;
    }

    if (argc==2)
    {
        ulong nbits = atoi(argv[1]);
        benchmark_nbits(nbits);
        return 0;
    }

    if (argc==4)
    {
        ulong nbits = atoi(argv[1]);
        ulong dim = atoi(argv[2]);
        ulong deg = atoi(argv[3]);
        benchmark_nbits_dim_deg(nbits,dim,deg);
    }

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
