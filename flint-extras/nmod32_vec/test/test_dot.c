#include "flint/test_helpers.h"
#include "flint/nmod.h"
#include "flint/nmod_vec.h"

#include "nmod32_vec.h"
#include <flint/ulong_extras.h>

// utility (nmod vec uniform random)
static inline
void _nmod32_vec_rand(n32_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = (uint)n_randint(state, mod.n);
}

TEST_FUNCTION_START(nmod32_vec_dot, state)
{
    flint_rand_set_seed(state, time(NULL), time(NULL)+12984125L);

    int i;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        // acc8 == 1 <=> room for accumulating 8 terms
        const int acc8 = (i < 2500);

        slong len = n_randint(state, 1000) + 1;

        const long pow2_61_sqrt = UWORD(1518500249);  // floor(2**30.5)
        const long pow2_31 = UWORD(1) << 31;
        ulong m = 0;
        if (acc8)
        {
            while (m == 0)
                m = n_randtest_not_zero(state) % pow2_61_sqrt;
        }
        else
        {
            while (m == 0)
                m = n_randtest_not_zero(state) % pow2_31;
        }

        nmod_t mod;
        nmod_init(&mod, m);

        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        // force computation of pow2_precomp, as flint's dot might not need it (and then will not compute it)
        ulong pow2_precomp;
        NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

        n32_ptr x = _nmod32_vec_init(len);
        n32_ptr y = _nmod32_vec_init(len);
        nn_ptr xx = _nmod_vec_init(len);
        nn_ptr yy = _nmod_vec_init(len);

        _nmod32_vec_rand(x, state, len, mod);
        _nmod32_vec_rand(y, state, len, mod);
        for (slong k = 0; k < len; k++)
        {
            xx[k] = x[k];
            yy[k] = y[k];
        }
        ulong correct = _nmod_vec_dot(xx, yy, len, mod, params);

        {  // dot_split
            ulong res = _nmod32_vec_dot_split(x, y, len, mod, pow2_precomp);

            if (res != correct)
            {
                flint_printf("%ld\n", i);
                TEST_FUNCTION_FAIL("dot_split, m = %wu, len = %wd\n", m, len);
            }
        }

        {  // dot_split_avx2
            ulong res = _nmod32_vec_dot_split_avx2(x, y, len, mod, pow2_precomp);

            if (res != correct)
            {
                flint_printf("%ld\n", i);
                TEST_FUNCTION_FAIL("dot_split_avx2, m = %wu, len = %wd\n", m, len);
            }
        }

        {  // dot_split_avx512
            ulong res = _nmod32_vec_dot_split_avx512(x, y, len, mod, pow2_precomp);

            if (res != correct)
            {
                flint_printf("%ld\n", i);
                TEST_FUNCTION_FAIL("dot_split_avx512, m = %wu, len = %wd\n", m, len);
            }
        }

        // seems to fail for large primes
        if (acc8)
        {  // dot_msolve_avx2
            ulong res = _nmod32_vec_dot_msolve_avx2(x, y, len, mod.n);

            if (res != correct)
            {
                flint_printf("%ld\n", i);
                TEST_FUNCTION_FAIL("dot_msolve_avx2, m = %wu, len = %wd\n", m, len);
            }
        }

        _nmod32_vec_clear(x);
        _nmod32_vec_clear(y);
        _nmod_vec_clear(xx);
        _nmod_vec_clear(yy);
    }

    TEST_FUNCTION_END(state);
}

int main()
{
    test_nmod32_vec_dot();
    return 0;
}
