/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_multiply.h"

/* a modulus of one of several kinds, chosen so as to exercise the
   different regimes of the plan: direct transform modulo an FFT prime,
   one to four CRT primes, composite moduli, tiny moduli */
static ulong _random_modulus_matmul(flint_rand_t state)
{
    switch (n_randint(state, 8))
    {
        case 0:  /* FFT prime of 30..50 bits with 2-adicity >= 12 */
        {
            ulong bits = 30 + n_randint(state, 21);   /* 30..50 */
            ulong k = 12 + n_randint(state, 8);       /* 2-adicity 12..19 */
            ulong c;
            slong tries = 0;
            do
            {
                /* c odd in [2^(bits-k-2), 2^(bits-k-1)) so that c*2^k+1 has bits bits */
                c = (UWORD(1) << (bits - k - 2)) + n_randint(state, UWORD(1) << (bits - k - 2));
                c |= 1;
                if (++tries > 100000)
                    return UWORD(469762049);   /* 7 * 2^26 + 1 */
            } while (!n_is_prime((c << k) + 1));
            return (c << k) + 1;
        }
        case 1:  /* one of FLINT's own 50-bit FFT primes */
            return UWORD(0x0003f00000000001);
        case 2:  /* tiny modulus */
            return 2 + n_randint(state, 30);
        case 3:  /* composite */
            return 2 + n_randint(state, UWORD(1) << (1 + n_randint(state, 63)));
        case 4:  /* full 64-bit prime */
            return n_randprime(state, 64, 1);
        case 5:  /* 63-bit prime */
            return n_randprime(state, 63, 1);
        default: /* random prime of 2..62 bits */
            return n_randprime(state, 2 + n_randint(state, 61), 1);
    }
}

TEST_FUNCTION_START(nmod_poly_mat_mul_sd_fft_matmul, state)
{
    slong i;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        const ulong modn = _random_modulus_matmul(state);
        slong m, k, n, lenA, lenB;
        nmod_poly_mat_t A, B, C1, C2;

        if (n_randint(state, 10) == 0)  /* larger sizes, small dimensions */
        {
            m = 1 + n_randint(state, 4);
            k = 1 + n_randint(state, 4);
            n = 1 + n_randint(state, 4);
            lenA = 1 + n_randint(state, 3000);
            lenB = 1 + n_randint(state, 3000);
        }
        else
        {
            m = n_randint(state, 9);
            k = n_randint(state, 9);
            n = n_randint(state, 9);
            lenA = n_randint(state, 100);
            lenB = n_randint(state, 100);
        }

        nmod_poly_mat_init(A, m, k, modn);
        nmod_poly_mat_init(B, k, n, modn);
        nmod_poly_mat_init(C1, m, n, modn);
        nmod_poly_mat_init(C2, m, n, modn);

        nmod_poly_mat_randtest(A, state, lenA);
        nmod_poly_mat_randtest(B, state, lenB);
        nmod_poly_mat_randtest(C2, state, 5);

        /* sprinkle zero entries, rows, columns */
        if (m > 0 && k > 0 && n_randint(state, 3) == 0)
        {
            slong r = n_randint(state, m), j;
            for (j = 0; j < k; j++)
                nmod_poly_zero(nmod_poly_mat_entry(A, r, j));
        }
        if (k > 0 && n > 0 && n_randint(state, 3) == 0)
        {
            slong c = n_randint(state, n), j;
            for (j = 0; j < k; j++)
                nmod_poly_zero(nmod_poly_mat_entry(B, j, c));
        }
        if (k > 0 && n_randint(state, 3) == 0)
        {
            slong l = n_randint(state, k), j;
            for (j = 0; j < m; j++)
                nmod_poly_zero(nmod_poly_mat_entry(A, j, l));
        }

        nmod_poly_mat_mul(C1, A, B);
        nmod_poly_mat_mul_sd_fft_matmul(C2, A, B);

        if (!nmod_poly_mat_equal(C1, C2))
            TEST_FUNCTION_FAIL("modn = %wu, m = %wd, k = %wd, n = %wd, lenA = %wd, lenB = %wd\n",
                               modn, m, k, n, lenA, lenB);

        /* aliasing */
        if (m == k && k == n)
        {
            nmod_poly_mat_set(C2, A);
            nmod_poly_mat_mul_sd_fft_matmul(C2, C2, B);
            if (!nmod_poly_mat_equal(C1, C2))
                TEST_FUNCTION_FAIL("aliasing C == A: modn = %wu, m = %wd, k = %wd, n = %wd\n",
                                   modn, m, k, n);
            nmod_poly_mat_set(C2, B);
            nmod_poly_mat_mul_sd_fft_matmul(C2, A, C2);
            if (!nmod_poly_mat_equal(C1, C2))
                TEST_FUNCTION_FAIL("aliasing C == B: modn = %wu, m = %wd, k = %wd, n = %wd\n",
                                   modn, m, k, n);
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C1);
        nmod_poly_mat_clear(C2);
    }

    TEST_FUNCTION_END(state);
}
