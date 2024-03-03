#include <stdlib.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h" // for rand
#include "nmod_poly_mat_multiply.h"


/*--------------------------------------------------------------*/
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
void test_nmod_poly_mat_mul(ulong m, ulong n, ulong p, ulong deg)
{
    flint_rand_t state;
    nmod_poly_mat_t A, B, C1, C2;
    mp_limb_t modulus;

    modulus = 1108307720798209;
    flint_randinit(state);

    nmod_poly_mat_init(A, m, n, modulus);
    nmod_poly_mat_init(B, n, p, modulus);
    nmod_poly_mat_init(C1, m, p, modulus);
    nmod_poly_mat_init(C2, m, p, modulus);

    nmod_poly_mat_rand(A, state, deg);
    nmod_poly_mat_rand(B, state, deg);
    nmod_poly_mat_rand(C1, state, deg);
    nmod_poly_mat_rand(C2, state, deg);

    nmod_poly_mat_mul(C1, A, B);
    nmod_poly_mat_mul_geometric(C2, A, B);

    if (!nmod_poly_mat_equal(C1, C2))
    {
        printf("error at m=%lu n=%lu p=%lu deg=%lu\n", m, n, p, deg);
        exit(-1);
    }
    
    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls tets                                              */
/*--------------------------------------------------------------*/
int main()
{
    ulong i;
    flint_set_num_threads(1);

    for (i = 1; i < 300; i += 1)
        test_nmod_poly_mat_mul(i, i, i, 2);

    for (i = 1; i < 300; i += 1)
        test_nmod_poly_mat_mul(i, i, i, 200);

    for (i = 1; i < 300; i += 5)
        test_nmod_poly_mat_mul(i, i, i, 2000);

    for (i = 1; i < 100; i += 20)
        test_nmod_poly_mat_mul(i, i, i, 20000);

    return 0;
}
