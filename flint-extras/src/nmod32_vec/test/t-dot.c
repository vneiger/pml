/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod32_vec.h"

TEST_FUNCTION_START(nmod32_vec_dot, state)
{
    flint_rand_set_seed(state, time(NULL), time(NULL)+12984125L);

    int i;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        slong len = n_randint(state, 2000) + 1;

        // temporarily disabling ifma (not finished, currently only for len multiple of 32)
        const int ifma = 0;
        if (ifma)
        {
            len += 32;
            len = (len / 32) * 32;
        }

        const long nmax8 = UWORD(1515531528);  // slightly less than floor(2**30.5)
        ulong m = 0;
        if (i < 1000)
            m = nmax8;
        else
            while (m == 0)
                m = n_randtest_not_zero(state) % nmax8;

        nmod_t mod;
        nmod_init(&mod, m);

        const dot_params_t params = _nmod_vec_dot_params(len, mod);
        // force computation of pow2_precomp, as flint's dot might not need it (and then will not compute it)
        ulong pow2_precomp;
        NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);
        ulong pow2_precomp_ifma;
        NMOD_RED(pow2_precomp_ifma, (UWORD(1) << 52), mod);

        n32_ptr x = _nmod32_vec_init(len);
        n32_ptr y = _nmod32_vec_init(len);
        nn_ptr xx = _nmod_vec_init(len);
        nn_ptr yy = _nmod_vec_init(len);

        if (i == 0)  /* maximum value for all entries */
        {
            for (slong k = 0; k < len; k++)
            {
                x[k] = mod.n - 1;
                y[k] = mod.n - 1;
            }
        }
        else
        {
            _nmod32_vec_rand(x, state, len, mod);
            _nmod32_vec_rand(y, state, len, mod);
        }

        for (slong k = 0; k < len; k++)
        {
            xx[k] = x[k];
            yy[k] = y[k];
        }
        ulong correct = _nmod_vec_dot(xx, yy, len, mod, params);

        {  // dot_split
            ulong res = _nmod32_vec_dot_split(x, y, len, mod, pow2_precomp);

            if (res != correct)
                TEST_FUNCTION_FAIL("dot_split, m = %wu, len = %wd\n", m, len);
        }

#if PML_HAVE_AVX2
        {  // dot_split_avx2
            ulong res = _nmod32_vec_dot_split_avx2(x, y, len, mod, pow2_precomp);

            if (res != correct)
                TEST_FUNCTION_FAIL("dot_split_avx2, m = %wu, len = %wd\n", m, len);
        }

        {  // dot_msolve_avx2
            ulong res = _nmod32_vec_dot_msolve_avx2(x, y, len, mod.n);

            if (res != correct)
                TEST_FUNCTION_FAIL("dot_msolve_avx2, m = %wu, len = %wd\n", m, len);
        }
#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
        {  // dot_split_avx512
            ulong res = _nmod32_vec_dot_split_avx512(x, y, len, mod, pow2_precomp);

            if (res != correct)
                TEST_FUNCTION_FAIL("dot_split_avx512, m = %wu, len = %wd\n", m, len);
        }

        if (ifma)
        {  // dot_ifma_avx512
            ulong res = _nmod32_vec_dot_ifma_avx2(x, y, len, mod, pow2_precomp_ifma);

            if (res != correct)
                TEST_FUNCTION_FAIL("dot_ifma_avx2, m = %wu, len = %wd\n", m, len);
        }

        if (ifma)
        {  // dot_ifma_avx512
            ulong res = _nmod32_vec_dot_ifma_avx512(x, y, len, mod, pow2_precomp_ifma);

            if (res != correct)
                TEST_FUNCTION_FAIL("dot_ifma_avx512, m = %wu, len = %wd\n", m, len);
        }
#endif  /* PML_HAVE_AVX512 */

        _nmod32_vec_clear(x);
        _nmod32_vec_clear(y);
        _nmod_vec_clear(xx);
        _nmod_vec_clear(yy);
    }

    TEST_FUNCTION_END(state);
}
