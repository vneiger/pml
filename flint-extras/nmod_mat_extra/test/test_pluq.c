#include <flint/ulong_extras.h>  // for n_randtest_prime

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

int main()
{
    // test from flint: nmod_mat/test/lu_classical
    flint_printf("Testing 5000 matrices; should take 5-15s\n");
    flint_printf("--> no message below means success\n");
    for (slong i = 0; i < 5000; i++)
    {
        FLINT_TEST_INIT(state);

        slong m = n_randint(state, 20);
        slong n = n_randint(state, 20);
        mp_limb_t mod = n_randtest_prime(state, 0);
        //mp_limb_t mod = 2;

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
            //nmod_mat_rand(A, state);

            // "First" variant
            {
                nmod_mat_t LU;
                nmod_mat_init_set(LU, A);
                slong * P = _perm_init(sizeof(slong) * m);
                slong * Q = _perm_init(sizeof(slong) * n);
                slong rank = nmod_mat_pluq(LU, P, Q);

                if (r != rank)
                {
                    printf("PLUQ: wrong rank\n");
                    return 1;
                }

                int result = check_pluq(LU, P, Q, A, rank);
                if (result != 0)
                {
                    if (result == 1)  // 1 --> wrong dimensions
                        flint_printf("PLUQ: wrong dimensions for LU\n");
                    else if (result == 2)  // 2 --> wrong shape, bottom left block is not zero
                        flint_printf("PLUQ: wrong shape for LU\n");
                    else if (result == 3)  // 3 --> P*L*U*Q is not equal to A
                        flint_printf("PLUQ: A  !=  P*L*U*Q\n");
                    return 1;
                }

                nmod_mat_clear(LU);
                flint_free(P);
                flint_free(Q);
            }

            // Crout variant
            {
                nmod_mat_t LU;
                nmod_mat_init_set(LU, A);
                slong * P = _perm_init(sizeof(slong) * m);
                slong * Q = _perm_init(sizeof(slong) * n);
                slong rank = nmod_mat_pluq_crout(LU, P, Q);

                if (r != rank)
                {
                    printf("PLUQ(Crout): Wrong rank\n");
                    return 1;
                }

                int result = check_pluq(LU, P, Q, A, rank);
                if (result != 0)
                {
                    if (result == 1)  // 1 --> wrong dimensions
                        flint_printf("PLUQ(Crout): wrong dimensions for LU\n");
                    else if (result == 2)  // 2 --> wrong shape, bottom left block is not zero
                        flint_printf("PLUQ(Crout): wrong shape for LU\n");
                    else if (result == 3)  // 3 --> P*L*U*Q is not equal to A
                        flint_printf("PLUQ(Crout): A  !=  P*L*U*Q\n");
                    return 1;
                }

                nmod_mat_clear(LU);
                flint_free(P);
                flint_free(Q);
            }

            nmod_mat_clear(A);
        }
        FLINT_TEST_CLEANUP(state);
    }
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
