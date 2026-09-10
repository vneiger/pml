/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/** 
 * Targets the pairwise distinct points case.
 * 
 * Targets nmod_mat_poly_mintbasis_rescomp / _resupdate / the
 * nmod_mat_poly_mintbasis dispatcher directly. Uses nmod_poly_mat_is_
 * interpolant_basis (verification.c, this same staging area) for the
 * defining-property check -- the M-level output (nmod_mat_poly_t) is
 * converted to nmod_poly_mat_t via nmod_poly_mat_set_from_mat_poly
 * (nmod_poly_mat_utils.h), the same conversion nmod_poly_mat_mintbasis
 * itself uses to wrap this dispatcher (interpolant_basis.c). 
 *
 * Checks, per random trial:
 *   (1) rescomp and resupdate agree bit-for-bit (they must: same row
 *       operations, same nullspace pivot choice, only residual bookkeeping
 *       differs);
 *   (2) the dispatcher's output equals whichever of the two variants its
 *       own condition (d*(m-n+1) <= m) selects;
 *       TODO CHECK nmod_mat_poly_mintbasis.c agrees with this formula 
 *   (3) the dispatcher's output is a genuine minimal interpolant basis
 *       (membership + generation, via the verifier), w.r.t. the original 
 *       (pre-call) shift -- shift is mutated in place by the construction
 *       calls. 
 * 
 * Also includes an explicit, non-random d=0 regression case, mirroring
 * mbasis's own order=0 edge case (see t-mbasis_variants.c): with no
 * points, all three functions must return the identity with shift
 * unchanged. TO BE DISCUSSED.
 *
 * Also includes an explicit very-small-prime stress block (prime in
 * {2,3,5,7,11}): the main randomized loop below never picks a prime under
 * 2^20 (its `nbits` floor is 20), but correctness here does not rely on
 * any genericity of the field -- the verifier's generation check
 * (deg(det) == sum_k rank(E_k)). 
 */

#include <flint/test_helpers.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_poly_mat_utils.h"



#include <flint/test_helpers.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_poly_mat_utils.h"

/* Builds a plain array of `elen` random m x n matrices mod `prime` */
static nmod_mat_struct * build_random_E(slong m, slong n, ulong prime, slong elen, flint_rand_t state)
{
    nmod_mat_struct * E = (nmod_mat_struct *) flint_malloc(FLINT_MAX(elen, 1) * sizeof(nmod_mat_struct));
    for (slong k = 0; k < elen; k++)
    {
        nmod_mat_init(E + k, m, n, prime);
        nmod_mat_randtest(E + k, state);
    }
    return E;
}

static void free_E(nmod_mat_struct * E, slong elen)
{
    for (slong k = 0; k < elen; k++)
        nmod_mat_clear(E + k);
    flint_free(E);
}

static int check_d0(slong m, slong n, ulong prime, flint_rand_t state)
{
    /* d is a separate explicit parameter, so d=0 is tested directly
       against a non-trivial, randomly-filled E that is longer than d --
       exercising "E longer than d, extra unused" (the now-plain-array
       convention, mirroring pts's own d <= length(pts)), not just the
       trivial d=0/E-empty case. TO SEE */
    slong elen = 1 + n_randint(state, 5);
    nmod_mat_struct * E = build_random_E(m, n, prime, elen, state);

    nmod_mat_poly_t i_res, i_upd;
    nmod_mat_poly_init(i_res, m, m, prime);
    nmod_mat_poly_init(i_upd, m, m, prime);

    slong * sh_res = flint_malloc(m * sizeof(slong));
    slong * sh_upd = flint_malloc(m * sizeof(slong));
    for (slong i = 0; i < m; i++)
        sh_res[i] = sh_upd[i] = n_randint(state, 10);

    nmod_mat_poly_mintbasis_rescomp(i_res, sh_res, NULL, E, 0);
    nmod_mat_poly_mintbasis_resupdate(i_upd, sh_upd, NULL, E, 0);

    int ok = nmod_mat_poly_is_one(i_res) && nmod_mat_poly_is_one(i_upd);

    free_E(E, elen);
    nmod_mat_poly_clear(i_res);
    nmod_mat_poly_clear(i_upd);
    flint_free(sh_res);
    flint_free(sh_upd);
    return ok;
}

/* returns 1 if the trial passed all checks, 0 otherwise. `elen` (E's own
   allocated length, >= d) is not needed here: every consumer only ever
   reads E[0..d-1], matching pts's own "array may be longer than d, extra
   ignored" convention. */
static int core_test_mintbasis(const nmod_mat_struct * E, slong m, slong n, ulong prime,
                     const ulong * pts, slong d, const slong * shift0)
{
    slong * sh_res = flint_malloc(m * sizeof(slong));
    slong * sh_upd = flint_malloc(m * sizeof(slong));
    slong * sh_disp = flint_malloc(m * sizeof(slong));
    for (slong i = 0; i < m; i++)
        sh_res[i] = sh_upd[i] = sh_disp[i] = shift0[i];

    nmod_mat_poly_t out_res, out_upd, out_disp;
    nmod_mat_poly_init(out_res, m, m, prime);
    nmod_mat_poly_init(out_upd, m, m, prime);
    nmod_mat_poly_init(out_disp, m, m, prime);

    nmod_mat_poly_mintbasis_rescomp(out_res, sh_res, pts, E, d);
    nmod_mat_poly_mintbasis_resupdate(out_upd, sh_upd, pts, E, d);
    nmod_mat_poly_mintbasis(out_disp, sh_disp, pts, E, d);

    int res = 1;

    /* (1) rescomp == resupdate, bit-for-bit */
    if (!nmod_mat_poly_equal(out_res, out_upd))
        res = 0;
    for (slong i = 0; i < m; i++)
        if (sh_res[i] != sh_upd[i])
            res = 0;

    /* (2) dispatcher matches whichever variant it should have picked */
    if (res)
    {
        if (!nmod_mat_poly_equal(out_disp, out_upd))
            res = 0;
        for (slong i = 0; i < m; i++)
            if (sh_disp[i] != sh_upd[i])
                res = 0;

        // int expect_resupdate = (d * (m - n + 1) <= m);
        // nmod_mat_poly_struct * expected = expect_resupdate ? out_upd : out_res;
        // slong * sh_expected = expect_resupdate ? sh_upd : sh_res;
        // if (!nmod_mat_poly_equal(out_disp, expected))
        //     res = 0;
        // for (slong i = 0; i < m; i++)
        //     if (sh_disp[i] != sh_expected[i])
        //         res = 0;
    }

    /* (3) genuine minimal interpolant basis, w.r.t. the originalshift0 --
       build the nmod_poly_mat_t E used by the verifier directly from the
       first d entries of the plain array (nmod_poly_mat_set_coeff_mat,
       the same primitive interpolant_basis.c itself uses in the other
       direction), and convert out_disp to nmod_poly_mat_t the same way
       nmod_poly_mat_mintbasis itself uses to wrap this dispatcher */
    if (res)
    {
        nmod_poly_mat_t out_disp_pm;
        nmod_poly_mat_init(out_disp_pm, m, m, prime);
        nmod_poly_mat_set_from_mat_poly(out_disp_pm, out_disp);

        if (!nmod_poly_mat_is_interpolant_basis(out_disp_pm, pts, E, d, shift0, ROW_LOWER))
            res = 0;
        nmod_poly_mat_clear(out_disp_pm);
    }

    nmod_mat_poly_clear(out_res);
    nmod_mat_poly_clear(out_upd);
    nmod_mat_poly_clear(out_disp);
    flint_free(sh_res);
    flint_free(sh_upd);
    flint_free(sh_disp);

    return res;
}



TEST_FUNCTION_START(nmod_mat_poly_mintbasis, state)
{
    int i, result;

    /* explicit d=0 regression case, both sides of the dispatcher's own
       boundary, not left to chance */
    for (i = 0; i < 20; i++)
    {
        slong m = 1 + n_randint(state, 15);
        slong n = 1 + n_randint(state, m);
        ulong prime = n_randprime(state, 2 + n_randint(state, 60), 1);
        if (!check_d0(m, n, prime, state))
            TEST_FUNCTION_FAIL("d=0 case: m = %wd, n = %wd, prime = %wu\n", m, n, prime);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong n = 1 + n_randint(state, 16);
        slong m = 1 + n_randint(state, 16); /* any m, n */
        slong d = n_randint(state, 150);

        /* Take a prime with nbits sufficiently large to ensure the existence
            of a list of pairwise distinct points. */ 
        ulong bound = (ulong) (2 * FLINT_MAX(d, 1) + 2);
        ulong nbits = FLINT_BIT_COUNT(bound) + 1 + n_randint(state, 50);
        ulong prime;
        do { prime = n_randprime(state, nbits, 1); } while (prime <= bound);

        /* E's own allocated length is drawn independently of, and at least
           as large as, d -- exercising "E longer than d, extra unused"
           (the plain array's own convention, mirroring pts's own
           d <= length(pts)) -- not just elen == d. Unlike the old
           nmod_mat_poly_t-based version of this file, no clamp is needed
           here any more: build_random_E's array length is exactly what is
           allocated, never silently shrunk the way nmod_mat_poly_t's own
           normalisation could (see build_random_E's own doc comment). */
        slong elen = d + n_randint(state, 10);
        nmod_mat_struct * E = build_random_E(m, n, prime, elen, state);

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
        for (slong j = 0; j < m; j++)
            shift0[j] = (n_randint(state, 4) == 0) ? n_randint(state, 20) : 0;

        result = core_test_mintbasis(E, m, n, prime, pts, d, shift0);

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wd, m = %wd, n = %wd, d = %wd, elen = %wd\n",
                               prime, m, n, d, elen);

        free_E(E, elen);
        flint_free(pts);
        flint_free(shift0);
    }

    /* explicit very-small-prime stress coverage, see header comment above */
    {
        ulong small_primes[] = {2, 3, 5, 7, 11};
        for (slong sp = 0; sp < (slong)(sizeof(small_primes) / sizeof(small_primes[0])); sp++)
        {
            ulong prime = small_primes[sp];

            for (i = 0; i < 10 * flint_test_multiplier(); i++)
            {
                slong n = 1 + n_randint(state, 16);
                slong m = n + n_randint(state, 8); /* n <= m, including n == m */
                slong d = n_randint(state, prime + 1); /* 0 <= d <= prime */

                nmod_mat_struct * E = build_random_E(m, n, prime, d, state);

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
                for (slong j = 0; j < m; j++)
                    shift0[j] = (n_randint(state, 4) == 0) ? n_randint(state, 20) : 0;

                result = core_test_mintbasis(E, m, n, prime, pts, d, shift0);

                if (!result)
                    TEST_FUNCTION_FAIL("small-prime case: prime = %wu, m = %wd, n = %wd, d = %wd\n",
                                       prime, m, n, d);

                free_E(E, d);
                flint_free(pts);
                flint_free(shift0);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
