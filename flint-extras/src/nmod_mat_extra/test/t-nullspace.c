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

#include "nmod_mat_extra.h"

/** checks X is in the nullspace of A, i.e X*A == 0 */
int is_in_nullspace(nmod_mat_t X, nmod_mat_t A)
{
    nmod_mat_t Z;
    nmod_mat_init(Z, X->r, A->c, A->mod.n);
    nmod_mat_mul(Z, X, A);
    if (nmod_mat_is_zero(Z) != 1)
    {
        printf("not in nullspace\n");
        return 0;
    }
    nmod_mat_clear(Z);
    return 1;
}

/** checks X generates the nullspace of A,
 *  and has full row rank, assuming X*A == 0 */
int basis_of_nullspace(nmod_mat_t X, nmod_mat_t A)
{
    slong rank = nmod_mat_rank(A);
    slong nullity = nmod_mat_rank(X);
    if (nullity != A->r - rank)
    {
        printf("nullspace rank not equal to nullity \n");
        return 0;
    }
    if (nullity != X->r)
    {
        printf("nullspace not full row rank\n");
        return 0;
    }
    return 1;
}

/** checks X is in reduced row echelon form, defining pivot as rightmost
 * nonzero entry in a row; assumes X has full row rank */
int is_rref(nmod_mat_t X)
{
    slong * pivot = malloc(X->r * sizeof(slong));
    for (slong i = 0; i < X->r; ++i)
    {
        pivot[i] = X->c - 1;
        while (pivot[i] >= 0 && nmod_mat_entry(X, i, pivot[i]) == 0)
            --pivot[i];
        if (pivot[i] < 0)
        {
            printf("found zero row in nullspace\n");
            return 0;
        }
        if (i>0 && pivot[i] <= pivot[i-1])
        {
            printf("pivots not increasing in nullspace\n");
            return 0;
        }
        for (slong j = i+1; j < X->r; ++j)
        {
            if (nmod_mat_entry(X,j,pivot[i]))
            {
                printf("entries not zero below pivot in nullspace\n");
                return 0;
            }
        }
    }
    return 1;
}

int check(slong field_prime, slong iterations, slong nrows, slong ncols, flint_rand_t state)
{
    int res = 1;

    nmod_mat_t A;
    nmod_mat_init(A, nrows, ncols, field_prime);

    slong k = 0;
    while (res && k < iterations)
    {
        nmod_mat_randfull(A, state);
        nmod_mat_t X;
        nmod_mat_left_nullspace(X, A);
        // DEBUG -->
        //printf("-- %ld: kernel OK --\n",k);
        //printf("%ldOK\n",k);
        //is_in_nullspace(X, A);
        //printf("%ldOK\n",k);
        //basis_of_nullspace(X, A);
        //printf("%ldOK\n",k);
        //is_rref(X);
        //printf("%ldOK\n",k);
        // <-- DEBUG
        res = res && (is_in_nullspace(X, A) && basis_of_nullspace(X, A) && is_rref(X));

        nmod_mat_clear(X);
        k++;
    }
    nmod_mat_clear(A);
    return res;
}

TEST_FUNCTION_START(nmod_mat_left_nullspace, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 63);
        slong prime = n_randprime(state, bits, 1);

        ulong m = 1 + n_randint(state, 100);
        ulong n = 1 + n_randint(state, 100);

        slong iterations = 10;

        result = check(prime, iterations, m, n, state);
        if (!result) TEST_FUNCTION_FAIL("m = %wu, n = %wu, p = %wu\n", m, n, prime);
    }

    TEST_FUNCTION_END(state);
}
