#include <stdlib.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/ulong_extras.h>  // for n_randtest_prime

#include "nmod_mat_extra.h"

int main(int argc, char ** argv)
{
    if (argc != 5)
    {
        flint_printf("Usage: %s nbits nrows ncols rank\n", argv[0]);
        return 0;
    }

    const slong nbits = atoi(argv[1]);
    const slong m = atoi(argv[2]);
    const slong n = atoi(argv[3]);

    slong try_rk = atoi(argv[4]);
    if (try_rk > FLINT_MIN(m, n))
    {
        flint_printf("Warning provided rank %d is greater than min(nrows,ncols) = %d, setting rank = %d\n", try_rk, FLINT_MIN(m,n), FLINT_MIN(m,n));
        try_rk = FLINT_MIN(m, n);
    }
    const slong rk = try_rk;

    // list of primes, primes[k] has bitlength k+2
    // sage: for nbits in range(4,61):
    // ....:     n = 2**(nbits-1) + 2**(nbits-2)
    // ....:     print(f"{n.next_prime()}      \t// {nbits} bits")
    slong primes[59] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741, 3221225473, 6442450967, 12884901893, 25769803799, 51539607599, 103079215111, 206158430209, 412316860441, 824633720837, 1649267441681, 3298534883417, 6597069766657, 13194139533349, 26388279066671, 52776558133303, 105553116266509, 211106232533047, 422212465066001, 844424930132057, 1688849860263953, 3377699720527897, 6755399441055827, 13510798882111519, 27021597764223071, 54043195528445957, 108086391056891941, 216172782113783843, 432345564227567621, 864691128455135281};

    flint_rand_t state;
    flint_randinit(state);

    nmod_mat_t mat, LU;
    nmod_mat_init(mat, m, n, primes[nbits-2]);
    if (rk == FLINT_MIN(m, n))
        // not 100% fine: not benchmarking non-generic rank profiles...
        // but makes random filling much faster
        nmod_mat_rand(mat, state);
    else
        nmod_mat_randrank_dense(mat, state, rk);
    nmod_mat_init(LU, m, n, primes[nbits-2]);

    double t;
    clock_t tt;
    long nb_iter;

    // warmup
    printf("warmup...\n");
    t = 0.0;
    nb_iter = 0;
    while (t < 2)
    {
        nmod_mat_set(LU, mat);
        slong * P = _perm_init(LU->r);
        slong * Q = _perm_init(LU->c);
        tt = clock();
        nmod_mat_pluq(LU, P, Q);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 1;
    }
    t /= nb_iter;

    printf("new pluq\tnew pluq crout\tlu_classical\tlu_delayed\tlu_recursive\n");

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_set(LU, mat);
        slong * P = _perm_init(LU->r);
        slong * Q = _perm_init(LU->c);
        tt = clock();
        nmod_mat_pluq(LU, P, Q);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _perm_clear(P);
        _perm_clear(Q);
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_set(LU, mat);
        slong * P = _perm_init(LU->r);
        slong * Q = _perm_init(LU->c);
        tt = clock();
        nmod_mat_pluq_crout(LU, P, Q);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _perm_clear(P);
        _perm_clear(Q);
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4e\t", t);


    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_set(LU, mat);
        slong * P = flint_malloc(LU->r * sizeof(slong));
        tt = clock();
        nmod_mat_lu_classical(P, LU, 0);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        flint_free(P);
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_set(LU, mat);
        slong * P = flint_malloc(LU->r * sizeof(slong));
        tt = clock();
        nmod_mat_lu_classical_delayed(P, LU, 0);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        flint_free(P);
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_set(LU, mat);
        slong * P = flint_malloc(LU->r * sizeof(slong));
        tt = clock();
        nmod_mat_lu_recursive(P, LU, 0);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        flint_free(P);
        nb_iter += 1;
    }
    t /= nb_iter;
    printf("%4e\n", t);

    nmod_mat_clear(mat);
    nmod_mat_clear(LU);
    flint_randclear(state);

}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
