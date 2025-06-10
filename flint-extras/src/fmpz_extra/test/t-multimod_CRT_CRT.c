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
/* computes a CRT with num_primes of n_bits size                */
/*--------------------------------------------------------------*/
int check_fmpz_multimod_CRT_CRT(ulong num_primes, ulong n_bits, flint_rand_t state)
{
    fmpz_t comb;
    fmpz_multimod_CRT_t mCRT; 
    nn_ptr primes, residues;

    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    fmpz_init(comb);

    nmod_vec_primes(primes, num_primes, n_bits);
    for (ulong i = 0; i < num_primes; i++)
    	residues[i] = n_randint(state, primes[i]);

    fmpz_multimod_CRT_init(mCRT, primes, num_primes);
    fmpz_multimod_CRT_CRT(comb, residues, mCRT);

    int res = 1;
    ulong i = 0;
    while (res && i < num_primes)
    {
        res = res && (fmpz_fdiv_ui(comb, primes[i]) == residues[i]);
        i++;
    }
    
    fmpz_multimod_CRT_clear(mCRT);

    fmpz_clear(comb);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);

    return res;
}

TEST_FUNCTION_START(fmpz_multimod_CRT_CRT, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong num_primes = 1;  /* corner case for i == 0 */
        if (i == 1)
            num_primes = 1000;  /* corner case for i == 1 */
        if (i > 1)
            num_primes = 1 + n_randint(state, 200);

        result = check_fmpz_multimod_CRT_CRT(num_primes, 60, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 60\n", num_primes);
        result = check_fmpz_multimod_CRT_CRT(num_primes, 50, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 50\n", num_primes);
        result = check_fmpz_multimod_CRT_CRT(num_primes, 29, state);
        if (!result) TEST_FUNCTION_FAIL("#primes = %wu, p_bits = 29\n", num_primes);
    }

    TEST_FUNCTION_END(state);
}
