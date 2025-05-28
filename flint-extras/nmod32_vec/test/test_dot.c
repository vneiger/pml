#include "flint/test_helpers.h"
#include "flint/nmod.h"
#include "flint/nmod_vec.h"

#include "nmod32_vec.h"

// utility (nmod vec uniform random)
static inline
void _nmod32_vec_rand(n32_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = (uint)n_randint(state, mod.n);
}

TEST_FUNCTION_START(nmod_vec_dot, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        ulong m = 0;
        ulong res;
        n32_ptr x, y;
        nn_ptr xx, yy;
        ulong correct;

        len = n_randint(state, 1000) + 1;
        while (m == 0)
            m = n_randtest_not_zero(state) % (UWORD(1) << 30);

        nmod_init(&mod, m);

        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        // force computation of pow2_precomp, as flint's dot might not need it (and then will not compute it)
        ulong pow2_precomp;
        NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

        x = _nmod32_vec_init(len);
        y = _nmod32_vec_init(len);
        xx = _nmod_vec_init(len);
        yy = _nmod_vec_init(len);

        _nmod32_vec_rand(x, state, len, mod);
        _nmod32_vec_rand(y, state, len, mod);
        for (slong k = 0; k < len; k++)
        {
            xx[k] = x[k];
            yy[k] = y[k];
        }
        correct = _nmod_vec_dot(xx, yy, len, mod, params);

        {  // dot2_half_avx2
            res = _nmod32_vec_dot2_split_avx2(x, y, len, mod, pow2_precomp);

            flint_printf("%ld\n", i);
            if ((ulong)res != correct)
            {
                TEST_FUNCTION_FAIL("i = %wu\n"
                                   "m = %wu\n"
                                   "len = %wd\n",
                                   i, m, len);
            }
        }

        {  // dot2_half_avx512
            res = _nmod32_vec_dot2_split_avx512(x, y, len, mod, pow2_precomp);

            flint_printf("%ld\n", i);
            if ((ulong)res != correct)
            {
                TEST_FUNCTION_FAIL("i = %wu\n"
                                   "m = %wu\n"
                                   "len = %wd\n",
                                   i, m, len);
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
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+12984125L);

    test_nmod_vec_dot();
    return 0;
}
