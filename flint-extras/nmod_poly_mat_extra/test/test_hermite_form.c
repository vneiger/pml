#include <stdlib.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>

#include "nmod_poly_mat_forms.h"

#include "testing_collection.h"

#define MAT(i,j) (mat->rows[i] + j)

//#define NOTRANS

// verify Hermite form
int verify_hermite_form(const nmod_poly_mat_t hnf, const slong * pivind, const nmod_poly_mat_t tsf, slong rk, const nmod_poly_mat_t mat, flint_rand_t state)
{
    if (! tsf)
    {
        printf("~~~ upper, rowwise ~~~ transformation not provided, skipping verification\n");
        return 1;
    }

    // testing correctness
    // 0. dimensions
    if (mat->r != hnf->r || mat->c != hnf->c || mat->r != tsf->r || mat->r != tsf->c)
    {
            printf("~~~ upper, rowwise ~~~ INCORRECT: dimension mismatch\n");
            return 0;
    }

    // empty matrices are fine (checked by is_hermite, but done here otherwise
    // tsf does not pass test on tsf_deg >= 0)
    if (mat->r == 0 || mat->c == 0)
        return 1;

    // 1. unimodular transformation
    // upper bound on degree of determinant
    slong tsf_deg = nmod_poly_mat_degree(tsf);
    if (tsf_deg < 0)
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is zero\n");
        return 0;
    }
    ulong degdet = tsf->r * (ulong)tsf_deg;
    if (mat->modulus > 100*degdet)
    {
        if (!nmod_poly_mat_is_unimodular_randomized(tsf, state))
        {
            printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is not unimodular (randomized test)\n");
            return 0;
        }
    }
    else
    {
        if (!nmod_poly_mat_is_unimodular(tsf))
        {
            printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is not unimodular\n");
            return 0;
        }
    }

    // 2. tsf * mat == hnf
    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_mul(tmp, tsf, mat);
    if (! nmod_poly_mat_equal(tmp, hnf))
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: tsf * mat != hnf\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    // 3. hnf is in Hermite form
    // first prune zero rows since this is a requirement of is_hermite
    // (note that zero rows are necessarily the ones at rk ... rdim-1)
    // TODO later simplify, once pruning zero rows  in algorithms
    // + verify [bottom part of tsf] * mat == 0
    nmod_poly_mat_clear(tmp);
    nmod_poly_mat_init(tmp, rk, hnf->c, hnf->modulus);
    for (slong i = 0; i < rk; i++)
        for (slong j = 0; j < hnf->c; j++)
            nmod_poly_set(nmod_poly_mat_entry(tmp, i, j), nmod_poly_mat_entry(hnf, i, j));
    if (! nmod_poly_mat_is_hermite(tmp, ROW_UPPER))
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: hnf not in Hermite form\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    // 4. check pivind
    slong * check_pivind = flint_malloc(rk * sizeof(slong));
    nmod_poly_mat_echelon_pivot_index(check_pivind, tmp, ROW_UPPER);
    for (slong i = 0; i < rk; i++)
    {
        if (pivind[i] != check_pivind[i])
        {
            printf("~~~ verify_hermite_form ~~~ INCORRECT: pivot index is wrong\n");
            nmod_poly_mat_clear(tmp);
            return 0;
        }
    }

    flint_free(check_pivind);
    nmod_poly_mat_clear(tmp);
    return 1;
}

int verify_nongeneric_column_rank_profile(const nmod_poly_mat_t mat)
{
    // take copy of mat
    nmod_poly_mat_t copy;
    nmod_poly_mat_init_set(copy, mat);
    // retrieve the rank of mat
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    const long rk = nmod_poly_mat_weak_popov_iter(copy, NULL, NULL, pivind, NULL, ROW_LOWER);
    // reset copy
    nmod_poly_mat_set(copy, mat);

    // retrieve the rank of first rk columns of mat
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, copy, 0, 0, mat->r, rk);
    const long rk_sub = nmod_poly_mat_weak_popov_iter(view, NULL, NULL, pivind, NULL, ROW_LOWER);

    // clear
    nmod_poly_mat_clear(view);
    nmod_poly_mat_clear(copy);
    flint_free(pivind);

    if (rk == rk_sub)
    {
            printf("~~~ upper, rowwise ~~~ INCORRECT: CRP is actually generic\n");
            return 0;
    }
    return 1;
}

// verify matrix rank profile
// (warning, tests all leading submatrices, can be slow!)
int verify_matrix_rank_profile(const nmod_poly_mat_t mat, const slong * mrp)
{
    if (mat->r == 0 || mat->c == 0)
        for (slong i = 0; i < mat->r; i++)
            if (mrp[i] != -1)
            {
                printf("~~~ verify_matrix_rank_profile ~~~ INCORRECT for empty matrix\n");
                return 0;
            }
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_t sub;
    for (slong i = 0; i < mat->r; i++)
    {
        for (slong j = 0; j < mat->c; j++)
        {
            // create submatrix
            nmod_poly_mat_init(sub, i, j, mat->modulus);
            for (slong ii = 0; ii < i; ii++)
                for (slong jj = 0; jj < j; jj++)
                    nmod_poly_set(nmod_poly_mat_entry(sub,ii,jj), MAT(ii, jj));
            // compute rank of submatrix
            slong rk_sub = nmod_poly_mat_weak_popov_iter(sub, NULL, NULL, pivind, NULL, ROW_LOWER);
            // compute rank of matrix rank profile
            slong rk_mrp = 0;
            for (slong k = 0; k < i; k++)
                if (mrp[k] >= 0 && mrp[k] < j)
                    rk_mrp++;
            // compare
            if (rk_sub != rk_mrp)
            {
                printf("~~~ verify_matrix_rank_profile ~~~ INCORRECT submatrix [:%ld,:%ld]\n", i, j);
                return 0;
            }
            nmod_poly_mat_clear(sub);
        }
    }
    flint_free(pivind);
    return 1;
}

// test one given input for hermite form
int core_test_hermite_form(const nmod_poly_mat_t mat, int time, flint_rand_t state)
{
    // initialize pivind list and rank profiles
    slong max_rank = FLINT_MIN(mat->r, mat->c);
    slong * pivind = flint_malloc(max_rank * sizeof(slong));
    slong * mrp = flint_malloc(mat->r * sizeof(slong));
    slong * mrp2 = flint_malloc(mat->r * sizeof(slong));

    // init copy of mat
    nmod_poly_mat_t hnf;
    nmod_poly_mat_init(hnf, mat->r, mat->c, mat->modulus);

    // init unimodular transformation tsf such that hnf = tsf * mat
    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, mat->r, mat->r, mat->modulus);

    int verif_hnf, verif_mrp;

    if (0)
    { // maxdeg, atomic
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t t0, t1;
        timeit_start(t0);
#ifdef NOTRANS
        slong rk = nmod_poly_mat_uref_maxdeg_atomic(hnf, NULL, pivind);
#else
        slong rk = nmod_poly_mat_uref_maxdeg_atomic(hnf, tsf, pivind);
#endif /* ifdef NOTRANS */
        timeit_stop(t0);
        timeit_start(t1);
#ifdef NOTRANS
        _normalize_uref(hnf, NULL, pivind, rk);
#else
        _normalize_uref(hnf, tsf, pivind, rk);
#endif /* ifdef NOTRANS */
        timeit_stop(t1);
        if (time)
            flint_printf("-- time (maxdeg - atomic - ur): %wd ms  (%wd + %wd)\n", t0->wall+t1->wall, t0->wall, t1->wall);
        timeit_start(t0);
#ifdef NOTRANS
        verif_hnf = verify_hermite_form(hnf, pivind, NULL, rk, mat, state);
#else
        verif_hnf = verify_hermite_form(hnf, pivind, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        if (!verif_hnf)
            printf("HNF -- maxdeg-atomic -- failure.\n");
        timeit_stop(t0);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", t0->wall);
    }

    if (0)
    { // revlex, xgcd
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t t0, t1;
        timeit_start(t0);
#ifdef NOTRANS
        slong rk = nmod_poly_mat_uref_revlex_xgcd(hnf, NULL, pivind, mrp);
#else
        slong rk = nmod_poly_mat_uref_revlex_xgcd(hnf, tsf, pivind, mrp);
#endif /* ifdef NOTRANS */
        timeit_stop(t0);
        timeit_start(t1);
#ifdef NOTRANS
        _normalize_uref(hnf, NULL, pivind, rk);
#else
        _normalize_uref(hnf, tsf, pivind, rk);
#endif /* ifdef NOTRANS */
        timeit_stop(t1);
        if (time)
            flint_printf("-- time (revlex - xgcd - ur): %wd ms  (%wd + %wd)\n", t0->wall+t1->wall, t0->wall, t1->wall);
        timeit_start(t0);
#ifdef NOTRANS
        verif_hnf = verify_hermite_form(hnf, pivind, NULL, rk, mat, state);
#else
        verif_hnf = verify_hermite_form(hnf, pivind, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        if (!verif_hnf)
            printf("HNF -- revlex-xgcd -- HNF failure.\n");
        verif_mrp = verify_matrix_rank_profile(mat, mrp);
        if (!verif_hnf)
            printf("HNF -- revlex-xgcd -- matrix rank profile failure.\n");
        timeit_stop(t0);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", t0->wall);
    }

    if (0)
    { // lex, xgcd
      // Warning! the mrp here is compared against the former
      // -- expect errors if former not computed
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t t0, t1;
        timeit_start(t0);
#ifdef NOTRANS
        slong rk = nmod_poly_mat_uref_lex_xgcd(hnf, NULL, pivind, mrp2);
#else
        slong rk = nmod_poly_mat_uref_lex_xgcd(hnf, tsf, pivind, mrp2);
#endif /* ifdef NOTRANS */
        timeit_stop(t0);
        timeit_start(t1);
#ifdef NOTRANS
        _normalize_uref(hnf, NULL, pivind, rk);
#else
        _normalize_uref(hnf, tsf, pivind, rk);
#endif /* ifdef NOTRANS */
        timeit_stop(t1);
        if (time)
            flint_printf("-- time (lex - xgcd - ur): %wd ms  (%wd + %wd)\n", t0->wall+t1->wall, t0->wall, t1->wall);
        timeit_start(t0);
#ifdef NOTRANS
        verif_hnf = verify_hermite_form(hnf, pivind, NULL, rk, mat, state);
#else
        verif_hnf = verify_hermite_form(hnf, pivind, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        if (!verif_hnf)
            printf("HNF -- lex-xgcd -- HNF failure.\n");
        verif_mrp = 1;
        for (slong i = 0; i < mat->r; i++)
            if (mrp[i] != mrp2[i])
                verif_mrp = 0;
        if (!verif_mrp)
            printf("HNF -- lex-xgcd -- matrix rank profile failure.\n");
        timeit_stop(t0);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", t0->wall);
    }

    if (0)
    { // Kannan-Bachem's algorithm
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t timer;
        timeit_start(timer);
#ifdef NOTRANS
        slong rk = nmod_poly_mat_hnf_ur_revlex_xgcd_delayed_zero(hnf, NULL, pivind, mrp2);
#else
        slong rk = nmod_poly_mat_hnf_ur_revlex_xgcd_delayed_zero(hnf, tsf, pivind, mrp2);
#endif /* ifdef NOTRANS */
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (revlex - xgcd - delayedzero - ur): %wd ms\n", timer->wall);
        timeit_start(timer);
#ifdef NOTRANS
        verif_hnf = verify_hermite_form(hnf, pivind, NULL, rk, mat, state);
#else
        verif_hnf = verify_hermite_form(hnf, pivind, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        if (!verif_hnf)
            printf("Hermite form -- revlex - xgcd - delayedzero - ur -- failure.\n");
        verif_mrp = 1;
        for (slong i = 0; i < mat->r; i++)
            if (mrp[i] != mrp2[i])
                verif_mrp = 0;
        if (!verif_mrp)
            printf("HNF -- revlex-xgcd-delayedzero-ur -- matrix rank profile failure.\n");
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", timer->wall);
    }

    { // Mulder-Storjohann's algorithm
        //slong * rrp = flint_malloc(mat->r * sizeof(slong));
        slong * rrp = NULL;
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t t0, t1;
        timeit_start(t0);
#ifdef NOTRANS
        slong rk = nmod_poly_mat_uref_matrixgcd_iter(hnf, NULL, pivind, rrp, NULL);
#else
        slong rk = nmod_poly_mat_uref_matrixgcd_iter(hnf, tsf, pivind, rrp, NULL);
#endif /* ifdef NOTRANS */
        timeit_stop(t0);
        if (rk <= 0 && !nmod_poly_mat_is_zero(mat))
        {
            verif_hnf = verify_nongeneric_column_rank_profile(mat);
            if (time)
                flint_printf("-- time (matrixgcd - ur): %wd ms  (!early exit!)\n", t0->wall);
        }
        else
        {
            timeit_start(t1);
#ifdef NOTRANS
            _normalize_uref(hnf, NULL, pivind, rk);
#else
            _normalize_uref(hnf, tsf, pivind, rk);
#endif /* ifdef NOTRANS */
            timeit_stop(t1);
            if (time)
                flint_printf("-- time (matrixgcd - ur): %wd ms  (%wd + %wd)\n", t0->wall+t1->wall, t0->wall, t1->wall);
            timeit_start(t0);
#ifdef NOTRANS
            verif_hnf = verify_hermite_form(hnf, pivind, NULL, rk, mat, state);
#else
            verif_hnf = verify_hermite_form(hnf, pivind, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        }
        if (!verif_hnf)
            printf("HNF -- matrix-gcd -- failure.\n");
        timeit_stop(t0);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", t0->wall);
        flint_free(rrp);
    }

    nmod_poly_mat_clear(hnf);
    nmod_poly_mat_clear(tsf);
    flint_free(pivind);
    flint_free(mrp);
    flint_free(mrp2);

    return verif_hnf;
}

//** Test against the whole testing collection */
int collection_test_hermite_form(slong iter, flint_rand_t state)
{
    // input matrix
    nmod_poly_mat_t mat;

    long total_nb_tests =
            iter // number of iterations
            * 5 // number of mats (currently 5)
            * _test_collection_nb_primes
            * _test_collection_nb_minidims
            * _test_collection_nb_minidims
            * _test_collection_nb_minidegs;

    printf("Launching testing collection (%ld cases)\n", total_nb_tests);

    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
        for (slong i_rdims = 0; i_rdims < _test_collection_nb_minidims; i_rdims++)
            for (slong i_cdims = 0; i_cdims < _test_collection_nb_minidims; i_cdims++)
                for (slong i_degs = 0; i_degs < _test_collection_nb_minidegs; i_degs++)
                    for (slong it = 0; it < iter; it++)
                    {
                        const long prime = _test_collection_primes[i_primes];
                        const long rdim = _test_collection_minirdims[i_rdims];
                        const long cdim = _test_collection_minicdims[i_cdims];
                        const long len = _test_collection_minidegs[i_degs];
                        printf("prime %ld, rdim %ld, cdim %ld, length %ld.\n", prime, rdim, cdim, len);

                        nmod_poly_mat_init(mat, rdim, cdim, prime);

                        _test_collection_mat_zero(mat);
                        if (! core_test_hermite_form(mat, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                        _test_collection_mat_uniform(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                        _test_collection_mat_test(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                        _test_collection_mat_sparse(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                        _test_collection_mat_rkdef(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "rkdef"); return 0; }

                        nmod_poly_mat_clear(mat);
                    }

    printf("--> Successful\n");
    return 1;
}

int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [order]\n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    if (argc != 1 && argc != 5)
        return 1;

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    int res = 0;

    if (argc == 1)
    {
        printf("launching test collection...\n");
        res = collection_test_hermite_form(10, state);
    }
    else if (argc == 5)
    {
        slong nbits = atoi(argv[1]);
        slong rdim = atoi(argv[2]);
        slong cdim = atoi(argv[3]);
        slong order = atoi(argv[4]);

        slong prime = n_randprime(state, nbits, 1);
        printf("Launching test with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
               \n\tlen = %ld...\n",prime,rdim,cdim,order);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, rdim, cdim, prime);
        nmod_poly_mat_rand(mat, state, order);

        res = core_test_hermite_form(mat, 1, state);

        nmod_poly_mat_t ref;
        nmod_poly_mat_init(ref, rdim, cdim, prime);
        nmod_poly_t den;
        nmod_poly_init(den, prime);
        timeit_t timer;
        timeit_start(timer);
        nmod_poly_mat_rref(ref, den, mat);
        timeit_stop(timer);
        nmod_poly_mat_clear(ref);
        nmod_poly_clear(den);
        flint_printf("-- time FFLU: %wd ms (for comparison)\n", timer->wall);

        nmod_poly_mat_clear(mat);
    }

    flint_rand_clear(state);

    if (res == 0)
    {
        printf("FAILURE\n");
        return 1;
    }
    else
    {
        printf("SUCCESS\n");
        return 0;
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
