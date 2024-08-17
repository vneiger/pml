#include <time.h>
#include <stdlib.h>
#include <flint/nmod_types.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

/*--------------------------------------------------------------*/
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
void time_nmod_poly_mat_mul(ulong m, ulong n, ulong p, ulong deg)
{
    flint_rand_t state;
    nmod_poly_mat_t A, B, C;
    double t;
    clock_t tt;
    long nb_iter;
    ulong modulus;

    modulus = 1108307720798209;
    flint_rand_init(state);

    nmod_poly_mat_init(A, m, n, modulus);
    nmod_poly_mat_init(B, n, p, modulus);
    nmod_poly_mat_init(C, m, p, modulus);

    nmod_poly_mat_rand(A, state, deg);
    nmod_poly_mat_rand(B, state, deg);
    nmod_poly_mat_rand(C, state, deg);

    printf("%lu\t%lu\t%lu\t%lu\t", m, n, p, deg);

    //t = 0.0;
    //nb_iter = 0;
    //while (t < 0.5)
    //{
    //    tt = clock();
    //    nmod_poly_mat_mul(C, A, B);
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    nb_iter += 1;
    //}
    //t /= nb_iter;
    //printf("%4g\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_poly_mat_mul_tft(C, A, B);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4g\t", t);
    

    printf("\n");
    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ulong i;
    flint_set_num_threads(1);
    
    if (argc==1)
    {
        for (i = 1; i < 300; i += 40)
            time_nmod_poly_mat_mul(i, i, i, 200);

        for (i = 1; i < 300; i += 40)
            time_nmod_poly_mat_mul(i, i, i, 2000);

        for (i = 1; i < 100; i += 20)
            time_nmod_poly_mat_mul(i, i, i, 20000);
    }

    else if (argc==5)
    {
        long rdim = atoi(argv[1]);
        long idim = atoi(argv[2]);
        long cdim = atoi(argv[3]);
        long deg = atoi(argv[4]);
        time_nmod_poly_mat_mul(rdim, idim, cdim, deg);
    }

    return 0;
}
