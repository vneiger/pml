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

#include "nmod_vec_extra.h"
#include "fmpz_extra.h"

/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* uses num_primes of size n_bits                               */
/* reduces 100 random integer modulo all primes                 */
/* input initially reduced modulo the product of primes         */
/*--------------------------------------------------------------*/
int check_fmpz_multimod_naive_reduce(ulong num_primes, ulong n_bits, flint_rand_t state)
{
    fmpz_multimod_naive_t mmod;
    nn_ptr primes, residues;

    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, n_bits);
    fmpz_multimod_naive_init(mmod, primes, num_primes);

    fmpz_t A;
    fmpz_init(A);

    fmpz_randtest(A, state, fmpz_bits(mmod->prod) - 1);
    fmpz_multimod_naive_reduce(residues, A, mmod);

    ulong res = 1;
    ulong i = 0;
    while (res && i < num_primes)
    {
        res = res && (residues[i] == fmpz_fdiv_ui(A, primes[i]));
        if (!res)
        {
            printf("A=");
            fmpz_print(A);
            printf("\n");
            printf("error with i=%lu, num_primes=%lu and n_bits=%lu\n", i, num_primes, n_bits);
            printf("%lu %lu\n", residues[i], fmpz_fdiv_ui(A, primes[i]));
            printf("p=%lu\n", primes[i]);
            exit(-1);
        }
        i++;
    }

    fmpz_clear(A);

    fmpz_multimod_naive_clear(mmod);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);

    return res;
}

TEST_FUNCTION_START(fmpz_multimod_naive_reduce, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong num_primes = 1;  /* corner case for i == 0 */
        if (i == 1)
            num_primes = 1000;  /* corner case for i == 1 */
        if (i > 1)
            num_primes = 1 + n_randint(state, 200);

        result = check_fmpz_multimod_naive_reduce(num_primes, 60, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 60\n", num_primes);
        result = check_fmpz_multimod_naive_reduce(num_primes, 50, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 50\n", num_primes);
        result = check_fmpz_multimod_naive_reduce(num_primes, 29, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 29\n", num_primes);
    }

    TEST_FUNCTION_END(state);
}
