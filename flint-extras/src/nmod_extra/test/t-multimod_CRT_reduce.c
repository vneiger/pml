/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>
#include <flint/fmpz.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"


/*----------------------------------------------------------------*/
/* reduces N elements of n_bits size modulo num_primes primes     */
/*----------------------------------------------------------------*/
int check_nmod_multimod_CRT_reduce(ulong N, ulong num_primes, ulong n_bits, flint_rand_t state)
{
    nmod_multimod_CRT_t CRT;
    nn_ptr *vec_residues;
    nn_ptr input;
    ulong i, j;
    ulong n;
    nmod_t mod;

    n = n_urandint(state, 1L << n_bits);
    nmod_init(&mod, n);

    input = _nmod_vec_init(N);
    _nmod_vec_rand(input, state, N, mod);

    nmod_multimod_CRT_init(CRT, n, num_primes);

    vec_residues = (nn_ptr *) malloc(num_primes * sizeof(nn_ptr *));
    for (i = 0; i < num_primes; i++)
        vec_residues[i] = _nmod_vec_init(N);

    int res = 1;
    nmod_multimod_CRT_reduce(vec_residues, input, N, CRT);
    i = 0;
    while (res && i < N)
    {
        j = 0;
        while (res && j < num_primes)
        {
            res = res && (vec_residues[j][i] == input[i] % CRT->mod_primes[j].n);
            j++;
        }
        i++;
    }

    _nmod_vec_clear(input);
    for (i = 0; i < num_primes; i++)
        _nmod_vec_clear(vec_residues[i]);
    free(vec_residues);
    nmod_multimod_CRT_clear(CRT);

    return res;
}

TEST_FUNCTION_START(nmod_multimod_CRT_reduce, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong N = 1 + n_randint(state, 500);

        result = check_nmod_multimod_CRT_reduce(N, 1, 30, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 1, 30\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 1, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 1, 60\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 2, 30, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 2, 30\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 2, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 2, 60\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 3, 30, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 3, 30\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 3, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 3, 60\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 4, 30, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 4, 30\n", N);

        result = check_nmod_multimod_CRT_reduce(N, 4, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 4, 60\n", N);
    }

    TEST_FUNCTION_END(state);
}
