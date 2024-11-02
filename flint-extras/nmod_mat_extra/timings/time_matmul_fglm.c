#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_mat_extra.h"

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
void time_nmod_mat_mul(ulong nbits, ulong dim1, ulong dim2, ulong block, ulong iter, ulong n, flint_rand_t state)
{
    printf("%lu\t%lu\t%lu\t%lu\t%lu\t", nbits, dim1, dim2, block, iter);

    nmod_mat_t mat, vec, res;
    nmod_t mod;

    nmod_init(&mod, n);

    nmod_mat_init(mat, dim1, dim2, mod.n);
    nmod_mat_init(vec, dim2, block, mod.n);
    nmod_mat_init(res, dim1, block, mod.n);

    nmod_mat_rand(mat, state);
    nmod_mat_rand(vec, state);
    nmod_mat_rand(res, state);

    double t;
    clock_t tt;
    long nb_iter;

    //t = 0.0;
    //nb_iter = 0;
    //while (t < 0.5)
    //{
    //    tt = clock();
    //    for (ulong it = 0; it < iter; it++)
    //        nmod_mat_mul(vec, mat, vec);
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    nb_iter += 1;
    //}
    //t /= nb_iter;
    //printf("%.1e\t", t);

    //t = 0.0;
    //nb_iter = 0;
    //while (t < 0.5)
    //{
    //    tt = clock();
    //    for (ulong it = 0; it < iter; it++)
    //        nmod_mat_mul_blas(res, mat, vec);
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    nb_iter += 1;
    //}
    //t /= nb_iter;
    //printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        for (ulong it = 0; it < iter; it++)
            nmod_mat_mul_newdot(res, mat, vec);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%.1e\t", t);

    if (n < (UWORD(1) << 30))
    {
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            for (ulong it = 0; it < iter; it++)
                nmod_mat_mul_small_modulus(res, mat, vec);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 1;
        }
        t /= nb_iter;
        printf("%.1e\t", t);
    }
    else
        printf("xxxxx\t");

    printf("\n");
    fflush(stdout);

    nmod_mat_clear(mat);
    nmod_mat_clear(vec);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("Block-Wied sequence generation over nmod:\n");
    //printf("  - mul: flint's nmod_mat_mul\n");
    printf("  - blas: flint's BLAS-based mul\n");
    printf("  - newdot: pml's mul using new nmod_vec_dot_product\n");
    printf("  - small_mod: pml's mul for small moduli using AVX-based nmod_vec_dot_product\n\n");
    printf("Args: nbits, matrix dim1, matrix dim2, block size, nb sequence terms\n\n");
    printf("      nbits == 0 --> test several nbits\n\n");
    printf("   (dim1 == nrows; dim2 == ncols; idea is dim1 == number of nontrivial rows)\n\n");

    if (argc == 6)
    {
        long nbits = atoi(argv[1]);
        long dim1 = atoi(argv[2]);
        long dim2 = atoi(argv[3]);
        long block = atoi(argv[4]);
        long iter = atoi(argv[5]);

        //printf("nbits\tdim1\tdim2\tblock\titer\tmul\tblas\tnewdot\tsmall_mod\n");
        printf("nbits\tdim1\tdim2\tblock\titer\tblas\tnewdot\tsmall_mod\n");
        if (nbits > 0)
        {
            time_nmod_mat_mul(nbits, dim1, dim2, block, iter, (1L << atoi(argv[1])) - 1, state);
        }
        else
        {
            const ulong mods[11] =
            {
                (1L << 3) + 1,
                (1L << 10) + 1,
                (1L << 20) + 1,
                (1L << 25) + 1,
                (1L << 29) + 1,
                (1L << 30) + 1,
                (1L << 31) + 1,
                (1L << 40) + 1,
                (1L << 50) + 1,
                (1L << 60) + 1,
                UWORD_MAX
            };
            const ulong nbits[11] = { 4, 11, 21, 26, 30, 31, 32, 41, 51, 61, 64 };

            for (slong i = 0; i < 11; i++)
            {
                time_nmod_mat_mul(nbits[i], dim1, dim2, block, iter, mods[i], state);
            }
        }
    }

    flint_rand_clear(state);
    return 0;
}
