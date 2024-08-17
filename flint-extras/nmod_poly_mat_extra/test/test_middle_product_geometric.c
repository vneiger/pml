#include <stdlib.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h" // for rand
#include "nmod_poly_mat_multiply.h"


/*--------------------------------------------------------------*/
/* middle product using different implementations               */
/*--------------------------------------------------------------*/
void test_nmod_poly_mat_middle_product(ulong m, ulong n, ulong p, ulong deg)
{
    flint_rand_t state;
    nmod_poly_mat_t A, B, C1, C2;
    ulong modulus;

    modulus = 1108307720798209;
    flint_rand_init(state);

    nmod_poly_mat_init(A, m, n, modulus);
    nmod_poly_mat_init(B, n, p, modulus);
    nmod_poly_mat_init(C1, m, p, modulus);
    nmod_poly_mat_init(C2, m, p, modulus);

    nmod_poly_mat_rand(A, state, deg);
    nmod_poly_mat_rand(B, state, 2 * deg - 1);
    nmod_poly_mat_rand(C1, state, deg);
    nmod_poly_mat_rand(C2, state, deg);

    nmod_poly_mat_middle_product_naive(C1, A, B, deg-1, deg-1);
    nmod_poly_mat_middle_product_geometric(C2, A, B, deg-1, deg-1);

    if (!nmod_poly_mat_equal(C1, C2))
    {
        printf("error at m=%lu n=%lu p=%lu deg=%lu\n", m, n, p, deg-1);
        nmod_poly_mat_print(A, "x");
        nmod_poly_mat_print(B, "x");
        nmod_poly_mat_print(C1, "x");
        nmod_poly_mat_print(C2, "x");
        exit(-1);
    }
    
    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls tets                                              */
/*--------------------------------------------------------------*/
int main()
{
    ulong i;
    flint_set_num_threads(1);

    for (i = 1; i < 100; i += 1)
        test_nmod_poly_mat_middle_product(1, 1, 1, i);

    for (i = 1; i < 300; i += 5)
        test_nmod_poly_mat_middle_product(2, 2, 2, i);

    for (i = 1; i < 100; i += 20)
        test_nmod_poly_mat_middle_product(i, i, i, 300);

    return 0;
}
