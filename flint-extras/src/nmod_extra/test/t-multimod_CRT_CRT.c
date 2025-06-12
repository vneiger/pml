/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/fmpz.h>
#include <flint/test_helpers.h>

#include "nmod_extra.h"


/*----------------------------------------------------------------*/
/* CRT in size num_bits, reduced mod n of size p_bits             */
/*----------------------------------------------------------------*/
int check_nmod_multimod_CRT_CRT(ulong N, ulong num_bits, ulong p_bits, flint_rand_t state)
{
    nmod_multimod_CRT_t CRT;
    nn_ptr *vec_residues;
    nn_ptr output;
    ulong i, j, num_primes;
    ulong n;
    fmpz * input;

    n = n_urandint(state, UWORD(1) << p_bits);
    input = (fmpz *) malloc(N * sizeof(fmpz));

    for (i = 0; i < N; i++)
    {
        fmpz_init(input + i);
        fmpz_randbits(input + i, state, num_bits);
        fmpz_abs(input + i, input + i);
    }

    if (num_bits > 4*49)
    {
        printf("not enough primes available\n");
        exit(-1);
    }
    num_primes = 1;
    if (num_bits > 49)
        num_primes = 2;
    if (num_bits > 2*49)
        num_primes = 3;
    if (num_bits > 3*49)
        num_primes = 4;

    nmod_multimod_CRT_init(CRT, n, num_primes);

    vec_residues = (nn_ptr *) malloc(num_primes * sizeof(nn_ptr));
    for (i = 0; i < num_primes; i++)
    {
        vec_residues[i] = _nmod_vec_init(N);
        for (j = 0; j < N; j++)
            vec_residues[i][j] = fmpz_fdiv_ui(input + j, CRT->mod_primes[i].n);
    }

    output = _nmod_vec_init(N);
    nmod_multimod_CRT_CRT(output, vec_residues, N, CRT);

    int res = 1;
    while (res && i < N)
    {
        res = res && (output[i] == fmpz_fdiv_ui(input + i, n));
        i++;
    }

    _nmod_vec_clear(output);

    for (i = 0; i < num_primes; i++)
        _nmod_vec_clear(vec_residues[i]);
    free(vec_residues);

    nmod_multimod_CRT_clear(CRT);

    for (i = 0; i < N; i++)
        fmpz_clear(input + i);
    free(input);

    return res;
}

TEST_FUNCTION_START(nmod_multimod_CRT_CRT, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong N = 1 + n_randint(state, 500);

        result = check_nmod_multimod_CRT_CRT(N, 30, 25, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 30, 25\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 50, 25, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 50, 25\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 100, 25, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 100, 25\n", N);
        /* TODO failing */
        /* result = check_nmod_multimod_CRT_CRT(N, 150, 25, state); */
        /* if (!result) TEST_FUNCTION_FAIL("N = %wu, 150, 25\n", N); */

        result = check_nmod_multimod_CRT_CRT(N, 30, 49, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 30, 49\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 50, 49, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 50, 49\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 100, 49, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 100, 49\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 150, 49, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 150, 49\n", N);

        result = check_nmod_multimod_CRT_CRT(N, 30, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 30, 60\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 50, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 50, 60\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 100, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 100, 60\n", N);
        result = check_nmod_multimod_CRT_CRT(N, 150, 60, state);
        if (!result) TEST_FUNCTION_FAIL("N = %wu, 150, 60\n", N);
    }

    TEST_FUNCTION_END(state);
}
