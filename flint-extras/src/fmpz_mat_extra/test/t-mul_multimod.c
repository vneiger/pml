#include <flint/fmpz_mat.h>
#include <flint/test_helpers.h>

#include "fmpz_mat_extra.h"


/*--------------------------------------------------------------*/
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
int test_fmpz_mat_mul(ulong m, ulong n, ulong p, ulong n_bits, flint_rand_t state)
{
    fmpz_mat_t A, B, C1, C2;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, n, p);
    fmpz_mat_init(C1, m, p);
    fmpz_mat_init(C2, m, p);

    fmpz_mat_randbits(A, state, n_bits);
    fmpz_mat_randbits(B, state, n_bits);
    fmpz_mat_randbits(C1, state, n_bits);
    fmpz_mat_randbits(C2, state, n_bits);

    fmpz_mat_mul(C1, A, B);
    fmpz_mat_mul_multimod(C2, A, B);

    int ret = fmpz_mat_equal(C1, C2);

    fmpz_mat_clear(C1);
    fmpz_mat_clear(C2);
    fmpz_mat_clear(B);
    fmpz_mat_clear(A);

    return ret;
}

TEST_FUNCTION_START(fmpz_mat_mul_multimod, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong n_bits = (1 + n_randint(state, 100)) * 200;
        ulong m = n_randint(state, 300);
        ulong n = n_randint(state, 300);
        ulong p = n_randint(state, 300);

        result = test_fmpz_mat_mul(m, n, p, n_bits, state);;

        if (!result)
            TEST_FUNCTION_FAIL(
                    "m = %wu, n = %wu, p = %wu\n"
                    "n_bits = %wu\n",
                    m, n, p, n_bits);
    }

    TEST_FUNCTION_END(state);
}

