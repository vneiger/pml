/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/ulong_extras.h>
#include <flint/test_helpers.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

#include "nmod32_vec.h"

TEST_FUNCTION_START(nmod32_vec_mdot, state)
{
    flint_rand_set_seed(state, time(NULL), time(NULL)+12984125L);

    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong nrows = n_randint(state, 100) + 1;
        slong len = n_randint(state, 1000) + 1;
        const slong stride = len;

        const long nmax8 = UWORD(1515531528);  // slightly less than floor(2**30.5)
        ulong m = 0;
        if (i < 100)
            m = nmax8;
        else
            while (m == 0)
                m = n_randtest_not_zero(state) % nmax8;

        nmod_t mod;
        nmod_init(&mod, m);

        n32_ptr vec = _nmod32_vec_init(len);
        n32_ptr mat = _nmod32_vec_init(nrows*len);
        nn_ptr vec64 = _nmod_vec_init(len);
        nmod_mat_t mat64;
        nmod_mat_init(mat64, nrows, len, mod.n);

        if (i == 0)  /* maximul value for all entries */
        {
            for (slong k = 0; k < len; k++)
                vec[k] = mod.n - 1;
            for (slong k = 0; k < nrows*len; k++)
                mat[k] = mod.n - 1;
        }
        else
        {
            _nmod32_vec_rand(vec, state, len, mod);
            _nmod32_vec_rand(mat, state, nrows*len, mod);
        }

        for (slong k = 0; k < len; k++)
            vec64[k] = vec[k];
        for (slong r = 0; r < nrows; r++)
            for (slong k = 0; k < len; k++)
                nmod_mat_entry(mat64, r, k) = mat[r*stride + k];

        nn_ptr correct = _nmod_vec_init(nrows);
        nmod_mat_mul_nmod_vec(correct, mat64, vec64, len);

        n32_ptr res = _nmod32_vec_init(nrows);

        {  // dot_split
            _nmod32_vec_mdot_split(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot_split, m = %wu, len = %wd\n", m, len);

#if PML_HAVE_AVX2
            _nmod32_vec_mdot_split_avx2(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot_split_avx2, m = %wu, len = %wd\n", m, len);
#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
            _nmod32_vec_mdot_split_avx512(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot_split_avx512, m = %wu, len = %wd\n", m, len);
#endif  /* PML_HAVE_AVX512 */
        }

        {  // dot2_split
            _nmod32_vec_mdot2_split(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot2_split, m = %wu, len = %wd\n", m, len);

#if PML_HAVE_AVX2
            _nmod32_vec_mdot2_split_avx2(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot2_split_avx2, m = %wu, len = %wd\n", m, len);

            _nmod32_vec_mdot3_split_avx2(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot3_split_avx2, m = %wu, len = %wd\n", m, len);
#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
            _nmod32_vec_mdot2_split_avx512(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot2_split_avx512, m = %wu, len = %wd\n", m, len);

            _nmod32_vec_mdot3_split_avx512(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot3_split_avx512, m = %wu, len = %wd\n", m, len);
#endif  /* PML_HAVE_AVX512 */
        }

        {  // dot_msolve
#if PML_HAVE_AVX2
            _nmod32_vec_mdot_msolve_via_dot_avx2(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot_msolve_avx2, m = %wu, len = %wd\n", m, len);
            _nmod32_vec_mdot_msolve_native_avx2(res, mat, vec, nrows, len, len, mod);
            for (slong k = 0; k < nrows; k++)
                if (res[k] != correct[k])
                    TEST_FUNCTION_FAIL("mdot_msolve_native_avx2, m = %wu, len = %wd\n", m, len);
#endif  /* PML_HAVE_AVX2 */
        }

        _nmod32_vec_clear(res);
        _nmod32_vec_clear(vec);
        _nmod32_vec_clear(mat);
        _nmod_vec_clear(vec64);
        _nmod_vec_clear(correct);
        nmod_mat_clear(mat64);
    }

    TEST_FUNCTION_END(state);
}
