/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/ulong_extras.h>  // for n_randtest_prime
#include <flint/test_helpers.h>  // for n_randtest_prime

#include "nmod_mat_extra.h"

// return values:
// 1 --> wrong dimensions
// 2 --> wrong shape, bottom left block is not zero
// 3 --> P*L*U*Q is not equal to A
int check_pluq(nmod_mat_t LU, slong * P, slong * Q, const nmod_mat_t A, slong rank)
{
    // check dimensions
    if (LU->r != A->r || LU->c != A->c)
        return 1;

    // check LU[rank:,rank:] is zero
    for (slong i = rank; i < LU->r; i++)
        for (slong j = rank; j < LU->c; j++)
            if (nmod_mat_entry(LU, i, j) != 0)
                return 2;

    // check P*L*U*Q == A
    nmod_mat_t L;
    nmod_mat_init(L, LU->r, LU->r, LU->mod.n);
    // compute L and U (the latter stored in LU)
    _nmod_mat_l_u_from_compactlu_zeroin(L, LU, LU, rank);
    // compute L*U, stored in LU
    nmod_mat_mul(LU, L, LU);
    // permute rows and columns
    _perm_inv(P, P, LU->r);
    nmod_mat_permute_rows(LU, P, NULL);
    _perm_inv(Q, Q, LU->c);
    nmod_mat_permute_columns(LU, Q, NULL);
    if (!nmod_mat_equal(LU, A))
        return 3;

    return 0;
}

TEST_FUNCTION_START(nmod_mat_pluq, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m = n_randint(state, 20);
        slong n = n_randint(state, 20);
        ulong mod = n_randtest_prime(state, 0);

        for (slong r = 0; r <= FLINT_MIN(m, n); r++)
        {
            nmod_mat_t A;
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);

            slong d;
            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                nmod_mat_randops(A, state, d);
            }

            // "First" variant
            {
                nmod_mat_t LU;
                nmod_mat_init_set(LU, A);
                slong * P = _perm_init(sizeof(slong) * m);
                slong * Q = _perm_init(sizeof(slong) * n);
                slong rank = nmod_mat_pluq(LU, P, Q);

                if (r != rank)
                {
                    TEST_FUNCTION_FAIL("PLUQ: wrong rank\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);
                }

                result = check_pluq(LU, P, Q, A, rank);
                if (result == 1)  // 1 --> wrong dimensions
                    TEST_FUNCTION_FAIL("PLUQ: wrong dimensions for LU\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);
                else if (result == 2)  // 2 --> wrong shape, bottom left block is not zero
                    TEST_FUNCTION_FAIL("PLUQ: wrong shape for LU\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);
                else if (result == 3)  // 3 --> P*L*U*Q is not equal to A
                    TEST_FUNCTION_FAIL("PLUQ: A  !=  P*L*U*Q\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);

                nmod_mat_clear(LU);
                flint_free(P);
                flint_free(Q);
            }

            // Crout variant
            if (0)
            {
                nmod_mat_t LU;
                nmod_mat_init_set(LU, A);
                slong * P = _perm_init(sizeof(slong) * m);
                slong * Q = _perm_init(sizeof(slong) * n);
                slong rank = nmod_mat_pluq_crout(LU, P, Q);

                if (r != rank)
                {
                    TEST_FUNCTION_FAIL("PLUQ(Crout): Wrong rank\n"
                                "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                m, n, r, rank);
                }

                int result = check_pluq(LU, P, Q, A, rank);
                if (result == 1)  // 1 --> wrong dimensions
                    TEST_FUNCTION_FAIL("PLUQ(Crout): wrong dimensions for LU\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);
                else if (result == 2)  // 2 --> wrong shape, bottom left block is not zero
                    TEST_FUNCTION_FAIL("PLUQ(Crout): wrong shape for LU\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);
                else if (result == 3)  // 3 --> P*L*U*Q is not equal to A
                    TEST_FUNCTION_FAIL("PLUQ(Crout): A  !=  P*L*U*Q\n"
                                       "m = %wu, n = %wu, r = %wu, rank = %wu\n",
                                       m, n, r, rank);

                nmod_mat_clear(LU);
                flint_free(P);
                flint_free(Q);
            }

            nmod_mat_clear(A);
        }
    }

    TEST_FUNCTION_END(state);
}
