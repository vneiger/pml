#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz_mat.h>

#include "fmpz_mat_extra.h"


/*--------------------------------------------------------------*/
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
void test_fmpz_mat_mul(ulong m, ulong n, ulong p, ulong n_bits)
{
    flint_rand_t state;
    fmpz_mat_t A, B, C1, C2;

    flint_rand_init(state);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, n, p);
    fmpz_mat_init(C1, m, p);
    fmpz_mat_init(C2, m, p);

    fmpz_mat_randbits(A, state, n_bits);
    fmpz_mat_randbits(B, state, n_bits);
    fmpz_mat_randbits(C1, state, n_bits);
    fmpz_mat_randbits(C2, state, n_bits);

    fmpz_mat_mul(C1, A, B);
    fmpz_mat_mul_waksman(C2, A, B);

    if (!fmpz_mat_equal(C1, C2))
    {
        printf("error at m=%lu n=%lu p=%lu n_bits=%lu\n", m, n, p, n_bits);
        exit(-1);
    }
    
    fmpz_mat_clear(C1);
    fmpz_mat_clear(C2);
    fmpz_mat_clear(B);
    fmpz_mat_clear(A);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls tets                                              */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ulong i;
    flint_set_num_threads(1);

    for (i = 1; i < 300; i += 40)
        test_fmpz_mat_mul(i, i, i, 200);

    for (i = 1; i < 300; i += 40)
        test_fmpz_mat_mul(i, i, i, 2000);

    for (i = 1; i < 100; i += 20)
        test_fmpz_mat_mul(i, i, i, 20000);

    for (i = 1; i < 10; i += 2)
        test_fmpz_mat_mul(i, i, i, 200000);

    return 0;
}
