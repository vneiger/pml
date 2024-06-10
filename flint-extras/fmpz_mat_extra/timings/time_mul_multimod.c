#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "fmpz_mat_extra.h"

/*--------------------------------------------------------------*/
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
void time_fmpz_mat_mul(ulong m, ulong n, ulong p, ulong n_bits)
{
    flint_rand_t state;
    fmpz_mat_t A, B, C;
    double t;
    clock_t tt;
    long nb_iter;

    flint_rand_init(state);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, n, p);
    fmpz_mat_init(C, m, p);

    fmpz_mat_randbits(A, state, n_bits);
    fmpz_mat_randbits(B, state, n_bits);
    fmpz_mat_randbits(C, state, n_bits);

    printf("%lu\t%lu\t%lu\t%lu\t", m, n, p, n_bits);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_mat_mul(C, A, B);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4g\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_mat_mul_multimod(C, A, B);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4g\t", t);


    printf("\n");
    fmpz_mat_clear(C);
    fmpz_mat_clear(B);
    fmpz_mat_clear(A);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ulong i;
    flint_set_num_threads(1);
    
    printf("#m\tn\tp\tn_bits\tmul\t\tmultimod\n");
    for (i = 1; i < 300; i += 40)
        time_fmpz_mat_mul(i, i, i, 200);

    for (i = 1; i < 300; i += 40)
        time_fmpz_mat_mul(i, i, i, 2000);

    for (i = 1; i < 100; i += 20)
        time_fmpz_mat_mul(i, i, i, 20000);

    for (i = 1; i < 10; i += 2)
        time_fmpz_mat_mul(i, i, i, 200000);

    return 0;
}
