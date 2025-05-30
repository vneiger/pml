#include <flint/ulong_extras.h>
#include <flint/test_helpers.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

#include "nmod32_vec.h"

// utility (nmod vec uniform random)
static inline
void _nmod32_vec_rand(n32_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = (uint)n_randint(state, mod.n);
}

TEST_FUNCTION_START(nmod32_vec_mdot, state)
{
    flint_rand_set_seed(state, time(NULL), time(NULL)+12984125L);

    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        // acc8 == 1 <=> room for accumulating 8 terms
        const int acc8 = (i < 500);

        slong nrows = n_randint(state, 100) + 1;
        slong len = n_randint(state, 1000) + 1;
        const slong stride = len;

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

        n32_ptr vec = _nmod32_vec_init(len);
        n32_ptr mat = _nmod32_vec_init(nrows*len);
        nn_ptr vec64 = _nmod_vec_init(len);
        nmod_mat_t mat64;
        nmod_mat_init(mat64, nrows, len, mod.n);

        _nmod32_vec_rand(vec, state, len, mod);
        _nmod32_vec_rand(mat, state, nrows*len, mod);
        for (slong k = 0; k < len; k++)
            vec64[k] = vec[k];
        for (slong r = 0; r < nrows; r++)
            for (slong k = 0; k < len; k++)
                nmod_mat_entry(mat64, r, k) = mat[r*stride + k];

        nn_ptr correct = _nmod_vec_init(nrows);
        nmod_mat_mul_nmod_vec(correct, mat64, vec64, len);

        {  // dot_split
            n32_ptr res = _nmod32_vec_init(nrows);
            n32_ptr res_avx2 = _nmod32_vec_init(nrows);
            _nmod32_vec_mdot_split(res, mat, vec, nrows, len, len, mod);
            _nmod32_vec_mdot_split_avx2(res_avx2, mat, vec, nrows, len, len, mod);
#if HAVE_AVX512   // TODO handle AVX flags
            n32_ptr res_avx512 = _nmod32_vec_init(nrows);
            _nmod32_vec_mdot_split_avx512(res_avx512, mat, vec, nrows, len, len, mod);

            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k] || res_avx2[k] != correct[k] || res_avx512[k] != correct[k])
#else
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k] || res_avx2[k] != correct[k])
#endif
                {
                    flint_printf("%ld\n", i);
                    TEST_FUNCTION_FAIL("mdot_split, m = %wu, len = %wd\n", m, len);
                }

            _nmod32_vec_clear(res);
            _nmod32_vec_clear(res_avx2);
#if HAVE_AVX512   // TODO handle AVX flags
            _nmod32_vec_clear(res_avx512);
#endif
        }

        {  // dot2_split
            n32_ptr res = _nmod32_vec_init(nrows);
            n32_ptr res_avx2 = _nmod32_vec_init(nrows);
            n32_ptr res3_avx2 = _nmod32_vec_init(nrows);
            _nmod32_vec_mdot2_split(res, mat, vec, nrows, len, len, mod);
            _nmod32_vec_mdot2_split_avx2(res_avx2, mat, vec, nrows, len, len, mod);
            _nmod32_vec_mdot3_split_avx2(res3_avx2, mat, vec, nrows, len, len, mod);
#if HAVE_AVX512   // TODO handle AVX flags
            n32_ptr res_avx512 = _nmod32_vec_init(nrows);
            n32_ptr res4_avx512 = _nmod32_vec_init(nrows);
            _nmod32_vec_mdot2_split_avx512(res_avx512, mat, vec, nrows, len, len, mod);
            _nmod32_vec_mdot4_split_avx512(res4_avx512, mat, vec, nrows, len, len, mod);

            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k] ||
                    res_avx2[k] != correct[k] ||
                    res_avx512[k] != correct[k] ||
                    res3_avx2[k] != correct[k] ||
                    res4_avx512[k] != correct[k])
#else
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k] ||
                    res_avx2[k] != correct[k] ||
                    res3_avx2[k] != correct[k])
#endif
                {
                    flint_printf("%ld\n", i);
                    TEST_FUNCTION_FAIL("mdot_split, m = %wu, len = %wd\n", m, len);
                }

            _nmod32_vec_clear(res);
            _nmod32_vec_clear(res_avx2);
            _nmod32_vec_clear(res3_avx2);
#if HAVE_AVX512   // TODO handle AVX flags
            _nmod32_vec_clear(res_avx512);
            _nmod32_vec_clear(res4_avx512);
#endif
        }

        if (acc8)
        {  // dot_msolve
            n32_ptr res = _nmod32_vec_init(nrows);
            n32_ptr res_native = _nmod32_vec_init(nrows);
            _nmod32_vec_mdot_msolve_via_dot_avx2(res, mat, vec, nrows, len, len, mod);
            _nmod32_vec_mdot_msolve_native_avx2(res_native, mat, vec, nrows, len, len, mod);

            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k] || res_native[k] != correct[k])
                {
                    flint_printf("%ld\n", i);
                    TEST_FUNCTION_FAIL("mdot_msolve, m = %wu, len = %wd\n", m, len);
                }
            _nmod32_vec_clear(res);
            _nmod32_vec_clear(res_native);
        }

        _nmod32_vec_clear(vec);
        _nmod32_vec_clear(mat);
        _nmod_vec_clear(vec64);
        _nmod_vec_clear(correct);
        nmod_mat_clear(mat64);
    }

    TEST_FUNCTION_END(state);
}

int main()
{
    test_nmod32_vec_mdot();
    return 0;
}

