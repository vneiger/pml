/*
   Copyright 2024 (C) Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
   */

#include <immintrin.h>
#include <flint/flint.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

#define NB_ITER 8192

void _nmod_vec_rand_not_zero(nn_ptr vec, flint_rand_t state, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
        vec[k] = 1 + n_randint(state, mod.n - 1);
}

void _nmod_vec_inv_naive(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
        res[k] = nmod_inv(vec[k], mod);
}


void sample_inv_naive(void * FLINT_UNUSED(arg), ulong count)
{
    const ulong len = NB_ITER;
    nn_ptr vec = (nn_ptr) flint_malloc(len*sizeof(ulong));
    nn_ptr res = _nmod_vec_init(len);
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS-2) + 3;  // 3...64
        ulong n = n_randprime(state, bits, 1);
        nmod_t mod;
        nmod_init(&mod, n);

        _nmod_vec_rand_not_zero(vec, state, len, mod);

        prof_start();
        _nmod_vec_inv_naive(res, vec, len, mod);
        prof_stop();
    }

    flint_free(vec);
    flint_free(res);
    FLINT_TEST_CLEAR(state);
}

void sample_inv(void * FLINT_UNUSED(arg), ulong count)
{
    const ulong len = NB_ITER;
    nn_ptr vec = (nn_ptr) flint_malloc(len*sizeof(ulong));
    nn_ptr res = _nmod_vec_init(len);
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS-2) + 3;  // 3...64
        ulong n = n_randprime(state, bits, 1);
        nmod_t mod;
        nmod_init(&mod, n);

        _nmod_vec_rand_not_zero(vec, state, len, mod);

        prof_start();
        _nmod_vec_inv(res, vec, len, mod);
        prof_stop();
    }

    flint_free(vec);
    flint_free(res);
    FLINT_TEST_CLEAR(state);
}

void sample_inv2(void * FLINT_UNUSED(arg), ulong count)
{
    const ulong len = NB_ITER;
    nn_ptr vec = (nn_ptr) flint_malloc(len*sizeof(ulong));
    nn_ptr res = _nmod_vec_init(len);
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS-2) + 3;  // 3...64
        ulong n = n_randprime(state, bits, 1);
        nmod_t mod;
        nmod_init(&mod, n);

        _nmod_vec_rand_not_zero(vec, state, len, mod);

        prof_start();
        _nmod_vec_inv2(res, vec, len, mod);
        prof_stop();
    }

    flint_free(vec);
    flint_free(res);
    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;

    flint_printf("nmod_vec_inv time:\n");

    prof_repeat(&min, &max, sample_inv_naive, NULL);
    flint_printf("   - naive: %.3f us / %.3f us\n",
            min, max);

    prof_repeat(&min, &max, sample_inv, NULL);
    flint_printf("   - ver1: %.3f us / %.3f us\n",
            min, max);

    prof_repeat(&min, &max, sample_inv2, NULL);
    flint_printf("   - ver2: %.3f us / %.3f us\n",
            min, max);

    return 0;
}

