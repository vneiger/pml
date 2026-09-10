/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/


/** Targets the pairwise distinct points case.
 *
 * Targets nmod_poly_mat_mintbasis and nmod_poly_mat_pmintbasis directly.
 * Uses nmod_poly_mat_is_interpolant_basis (verification.c, this same
 * staging area) for the defining-property check.
 * Also checks that nmod_poly_mat_pmintbasis, whose base case IS a direct
 * call to nmod_poly_mat_mintbasis (see interpolant_basis.c), agrees with
 * it bit-for-bit when d does not exceed the recursion threshold.
 *
 * Shift and E diversity: reuses testing_collection.h directly (this file's
 * final home, src/nmod_poly_mat_extra/test/, is the same directory t-mbasis.c
 * / t-pmbasis.c already pull it from), picking one of its 8 shift forms and
 * one of its 7 matrix forms (zero, uniform, unbalanced row/column degree,
 * randtest, sparse, rank-deficient) at random per trial via gen_shift/gen_E
 * below.
 * 
 * Also includes an explicit very-small-prime stress block (prime in
 * {2,3,5,7,11}), mirroring t-mintbasis.c's own -- see that file's header
 * comment for why tiny fields are worth testing directly rather than relying
 * on the main loop's nbits >= 20 floor.
 */

#include <flint/test_helpers.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_interpolant.h"
#include "testing_collection.h"
#include "interpolant_test_utils.h"

/* d is passed explicitly rather than read off E: matches pts's own
   "array may be longer than d, extra ignored" convention -- see this
   header's "Conventions" section in nmod_poly_mat_interpolant.h. */
static int core_test_pmintbasis(const nmod_mat_struct * E, slong m, slong n, ulong prime,
                     const ulong * pts, slong d, const slong * shift0)
{
    slong * sh_pm = flint_malloc(m * sizeof(slong));
    for (slong i = 0; i < m; i++)
        sh_pm[i] = shift0[i];

    nmod_poly_mat_t out_pm;
    nmod_poly_mat_init(out_pm, m, m, prime);
    nmod_poly_mat_pmintbasis(out_pm, sh_pm, pts, E, d);

    int res = 1;

    /* base-case agreement: pmintbasis's own base case is a direct call to
       mintbasis, so when d is within the threshold they must agree
       bit-for-bit */
    if (d <= PMINTBASIS_THRES)
    {
        slong * sh_m = flint_malloc(m * sizeof(slong));
        for (slong i = 0; i < m; i++)
            sh_m[i] = shift0[i];
        nmod_poly_mat_t out_m;
        nmod_poly_mat_init(out_m, m, m, prime);
        nmod_poly_mat_mintbasis(out_m, sh_m, pts, E, d);

        if (!nmod_poly_mat_equal(out_pm, out_m))
            res = 0;
        for (slong i = 0; i < m; i++)
            if (sh_pm[i] != sh_m[i])
                res = 0;

        nmod_poly_mat_clear(out_m);
        flint_free(sh_m);
    }

    /* genuine minimal interpolant basis, w.r.t. the original shift0. */
    if (res && !nmod_poly_mat_is_interpolant_basis(out_pm, pts, E, d, shift0, ROW_LOWER))
        res = 0;

    nmod_poly_mat_clear(out_pm);
    flint_free(sh_pm);

    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_pmintbasis, state)
{
    int i, result;

    /* explicit d=0 regression case -- E is never dereferenced when d=0
       (both nmod_poly_mat_pmintbasis and the M-level dispatcher it falls
       back to return before touching E), so NULL is passed directly. */
    for (i = 0; i < 20; i++)
    {
        slong n = 1 + n_randint(state, 12);
        slong m = 1 + n_randint(state, 12);
        ulong prime = n_randprime(state, 2 + n_randint(state, 60), 1);

        nmod_poly_mat_t out;
        nmod_poly_mat_init(out, m, m, prime);
        slong * shift = flint_malloc(m * sizeof(slong));
        for (slong j = 0; j < m; j++)
            shift[j] = n_randint(state, 10);

        nmod_poly_mat_pmintbasis(out, shift, NULL, NULL, 0);
        if (!nmod_poly_mat_is_one(out))
            TEST_FUNCTION_FAIL("d=0 case: m = %wd, n = %wd, prime = %wu\n", m, n, prime);

        nmod_poly_mat_clear(out);
        flint_free(shift);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong n = 1 + n_randint(state, 16);
        slong m = 1 + n_randint(state, 16); /* any m, n */
        slong d = n_randint(state, 250); /* well beyond a small overridden threshold */

        /* nbits' floor must guarantee some prime of that bit length exceeds
           bound, or the retry loop below never terminates (nbits is fixed
           before the retry starts, so it can't grow to escape a bad draw).
           Tying the floor to bound -- instead of a fixed constant such as
           the 20 this replaces -- lets nbits range down as low as 3-4 bits
           when d is small, closing a real gap this file used to leave
           untested (nothing between the 20-59 bit main range and the
           explicit {2,3,5,7,11} block below). */
        ulong bound = (ulong) (2 * FLINT_MAX(d, 1) + 2);
        ulong nbits = FLINT_BIT_COUNT(bound) + 1 + n_randint(state, 50);
        ulong prime;
        do { prime = n_randprime(state, nbits, 1); } while (prime <= bound);

        nmod_poly_mat_t E_poly;
        nmod_poly_mat_init(E_poly, m, n, prime);
        gen_E(E_poly, d, state);
        nmod_mat_struct * E = poly_mat_to_array(E_poly, d);
        nmod_poly_mat_clear(E_poly);

        ulong * pts = flint_malloc(FLINT_MAX(d, 1) * sizeof(ulong));
        for (slong k = 0; k < d; k++)
        {
            int distinct;
            do
            {
                pts[k] = n_randint(state, prime);
                distinct = 1;
                for (slong j = 0; j < k; j++)
                    if (pts[j] == pts[k])
                        distinct = 0;
            } while (!distinct);
        }

        slong * shift0 = flint_malloc(m * sizeof(slong));
        gen_shift(shift0, m, d, state);

        result = core_test_pmintbasis(E, m, n, prime, pts, d, shift0);

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wd, m = %wd, n = %wd, d = %wd\n", prime, m, n, d);

        free_array(E, d);
        flint_free(pts);
        flint_free(shift0);
    }

    /* explicit very-small-prime stress coverage, mirroring t-mintbasis.c's
       own block (see its header comment for why: the main loop's `nbits`
       floor  never picks these, and the verifier's generation check is
       an exact identity regardless of genericity, so tiny fields -- where
       coincidental rank degeneracies at a point are the rule rather than
       the exception -- are worth exercising directly). `d` is capped at
       `prime` (the maximum number of distinct points available). */
    {
        ulong small_primes[] = {2, 3, 5, 7, 11};
        for (slong sp = 0; sp < (slong)(sizeof(small_primes) / sizeof(small_primes[0])); sp++)
        {
            ulong prime = small_primes[sp];

            for (i = 0; i < 10 * flint_test_multiplier(); i++)
            {
                slong n = 1 + n_randint(state, 4);
                slong m = 1 + n_randint(state, 4); /* any m, n, including n == m */
                slong d = n_randint(state, prime + 1); /* 0 <= d <= prime */

                nmod_poly_mat_t E_poly;
                nmod_poly_mat_init(E_poly, m, n, prime);
                gen_E(E_poly, d, state);
                nmod_mat_struct * E = poly_mat_to_array(E_poly, d);
                nmod_poly_mat_clear(E_poly);

                ulong * pts = flint_malloc(FLINT_MAX(d, 1) * sizeof(ulong));
                for (slong k = 0; k < d; k++)
                {
                    int distinct;
                    do
                    {
                        pts[k] = n_randint(state, prime);
                        distinct = 1;
                        for (slong j = 0; j < k; j++)
                            if (pts[j] == pts[k])
                                distinct = 0;
                    } while (!distinct);
                }

                slong * shift0 = flint_malloc(m * sizeof(slong));
                gen_shift(shift0, m, d, state);

                result = core_test_pmintbasis(E, m, n, prime, pts, d, shift0);

                if (!result)
                    TEST_FUNCTION_FAIL("small-prime case: prime = %wu, m = %wd, n = %wd, d = %wd\n",
                                       prime, m, n, d);

                free_array(E, d);
                flint_free(pts);
                flint_free(shift0);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
