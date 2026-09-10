/*
    Copyright (C) 2026 Gilles Villard 
    
    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Targets nmod_mat_poly_mbasis_rescomp / _resupdate / the nmod_mat_poly_mbasis
 * dispatcher directly (the mat_poly level), complementing the existing
 * nmod_poly_mat_extra/test/t-mbasis.c (which exercises the nmod_poly_mat_mbasis
 * wrapper via the defining-property check nmod_poly_mat_is_approximant_basis).
 *
 * Checks, per random trial:
 *   (1) rescomp and resupdate agree bit-for-bit (they must: same row
 *       operations, same nullspace pivot choice, only residual bookkeeping
 *       differs) -- this is the property specific to introducing a second
 *       variant, not just "is it a valid basis";
 *   (2) the dispatcher's output equals whichever of the two variants its
 *       own condition (2*cdim > rdim) selects;
 *   (3) converted to nmod_poly_mat_t, the result is a genuine shift-ordered
 *       weak Popov approximant basis (nmod_poly_mat_is_approximant_basis),
 *       independent of any implementation comparison.
 *
 * FIXME next paragraph and corresponding code to be handled within issue 52
 * 
 * Also includes an explicit, non-random order=0 regression case: a
 * pre-existing bug  left `nsbas`
 * uninitialized when order=0 (the main loop never runs), making the
 * unconditional nmod_mat_clear(nsbas) at the end of the function undefined
 * behaviour -- intermittent, not a deterministic crash, since it depends on
 * whatever happens to be on the stack/heap at that point. Fixed in both
 * rescomp and resupdate by giving nsbas a trivial valid init before the
 * loop. This case is run first, unconditionally, rather than left to chance
 * inside the random loop below. */

#include <flint/test_helpers.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_utils.h"   /* nmod_poly_mat_set_from_mat_poly */

static int check_order0(slong rdim, slong cdim, ulong prime, flint_rand_t state)
{
    nmod_mat_poly_t matp, app_res, app_upd;
    nmod_mat_poly_init(matp, rdim, cdim, prime);
    nmod_mat_poly_rand(matp, state, 1 + n_randint(state, 10));
    nmod_mat_poly_init(app_res, rdim, rdim, prime);
    nmod_mat_poly_init(app_upd, rdim, rdim, prime);

    slong * sh_res = flint_malloc(rdim * sizeof(slong));
    slong * sh_upd = flint_malloc(rdim * sizeof(slong));
    for (slong i = 0; i < rdim; i++)
        sh_res[i] = sh_upd[i] = n_randint(state, 10);

    nmod_mat_poly_mbasis_rescomp(app_res, sh_res, matp, 0);
    nmod_mat_poly_mbasis_resupdate(app_upd, sh_upd, matp, 0);

    int ok = nmod_mat_poly_is_one(app_res) && nmod_mat_poly_is_one(app_upd);

    nmod_mat_poly_clear(matp);
    nmod_mat_poly_clear(app_res);
    nmod_mat_poly_clear(app_upd);
    flint_free(sh_res);
    flint_free(sh_upd);
    return ok;
}

static int core_test_mbasis_variants(nmod_mat_poly_t matp, slong order, const slong * shift)
{
    const slong rdim = matp->r;
    const slong cdim = matp->c;

    slong * sh_res = flint_malloc(rdim * sizeof(slong));
    slong * sh_upd = flint_malloc(rdim * sizeof(slong));
    slong * sh_disp = flint_malloc(rdim * sizeof(slong));
    for (slong i = 0; i < rdim; i++)
        sh_res[i] = sh_upd[i] = sh_disp[i] = shift[i];

    nmod_mat_poly_t app_res, app_upd, app_disp;
    nmod_mat_poly_init(app_res, rdim, rdim, matp->mod.n);
    nmod_mat_poly_init(app_upd, rdim, rdim, matp->mod.n);
    nmod_mat_poly_init(app_disp, rdim, rdim, matp->mod.n);

    nmod_mat_poly_mbasis_rescomp(app_res, sh_res, matp, order);
    nmod_mat_poly_mbasis_resupdate(app_upd, sh_upd, matp, order);
    nmod_mat_poly_mbasis(app_disp, sh_disp, matp, order);

    int res = 1;

    /* (1) rescomp == resupdate, bit-for-bit */
    if (!nmod_mat_poly_equal(app_res, app_upd))
        res = 0;
    for (slong i = 0; i < rdim; i++)
        if (sh_res[i] != sh_upd[i])
            res = 0;

    /* (2) dispatcher matches whichever variant it should have picked */
    if (res)
    {
        int expect_resupdate = (2 * cdim > rdim);
        nmod_mat_poly_struct * expected = expect_resupdate ? app_upd : app_res;
        slong * sh_expected = expect_resupdate ? sh_upd : sh_res;
        if (!nmod_mat_poly_equal(app_disp, expected))
            res = 0;
        for (slong i = 0; i < rdim; i++)
            if (sh_disp[i] != sh_expected[i])
                res = 0;
    }

    /* (3) genuine approximant basis property, converted to poly_mat */
    if (res)
    {
        nmod_poly_mat_t pmat, appbas;
        nmod_poly_mat_init(pmat, rdim, cdim, matp->mod.n);
        nmod_poly_mat_set_from_mat_poly(pmat, matp);
        nmod_poly_mat_init(appbas, rdim, rdim, matp->mod.n);
        nmod_poly_mat_set_from_mat_poly(appbas, app_disp);

        slong * cshift = flint_malloc(rdim * sizeof(slong));
        for (slong i = 0; i < rdim; i++)
            cshift[i] = shift[i];

        if (!nmod_poly_mat_is_approximant_basis(appbas, pmat, order, cshift, ROW_LOWER))
            res = 0;

        nmod_poly_mat_clear(pmat);
        nmod_poly_mat_clear(appbas);
        flint_free(cshift);
    }

    nmod_mat_poly_clear(app_res);
    nmod_mat_poly_clear(app_upd);
    nmod_mat_poly_clear(app_disp);
    flint_free(sh_res);
    flint_free(sh_upd);
    flint_free(sh_disp);

    return res;
}

TEST_FUNCTION_START(nmod_mat_poly_mbasis_variants, state)
{
    int i, result;

    /* explicit order=0 regression case, both sides of the cdim/rdim
       dispatch boundary, not left to chance */
    for (i = 0; i < 20; i++)
    {
        slong rdim = 1 + n_randint(state, 15);
        slong cdim = 1 + n_randint(state, 15);
        ulong prime = n_randprime(state, 2 + n_randint(state, 60), 1);
        if (!check_order0(rdim, cdim, prime, state))
            TEST_FUNCTION_FAIL("order=0 case: rdim = %wd, cdim = %wd, prime = %wu\n",
                               rdim, cdim, prime);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 63);
        ulong rdim = 1 + n_randint(state, 15);
        ulong cdim = 1 + n_randint(state, 15);
        ulong order = n_randint(state, 30);
        ulong len = n_randint(state, 30);

        slong prime = n_randprime(state, nbits, 1);

        nmod_mat_poly_t matp;
        nmod_mat_poly_init(matp, rdim, cdim, prime);
        nmod_mat_poly_rand(matp, state, len);

        slong * shift = flint_malloc(rdim * sizeof(slong));
        for (slong j = 0; j < rdim; j++)
            shift[j] = (n_randint(state, 4) == 0) ? n_randint(state, 20) : 0;

        result = core_test_mbasis_variants(matp, order, shift);

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wd, rdim = %wu, cdim = %wu, order = %wu, len = %wu\n",
                               prime, rdim, cdim, order, len);

        nmod_mat_poly_clear(matp);
        flint_free(shift);
    }

    TEST_FUNCTION_END(state);
}
