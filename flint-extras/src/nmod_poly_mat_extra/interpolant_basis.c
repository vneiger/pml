/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/


#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

/** nmod_mat_poly_mintbasis takes E as a plain, flat array of at 
 *     least `d` already-initialized nmod_mat_t's (nmod_mat_struct *)*/ 
void nmod_poly_mat_mintbasis(nmod_poly_mat_t intbas,
                             slong * shift,
                             const ulong * pts,
                             const nmod_mat_struct * E,
                             slong d)
{
    const slong m = intbas->r;
    nmod_mat_poly_t intbasmp;
    nmod_mat_poly_init(intbasmp, m, m, intbas->modulus);
    nmod_mat_poly_mintbasis(intbasmp, shift, pts, E, d);
    nmod_poly_mat_set_from_mat_poly(intbas, intbasmp);
    nmod_mat_poly_clear(intbasmp);
}



void nmod_poly_mat_pmintbasis(nmod_poly_mat_t intbas,
                              slong * shift,
                              const ulong * pts,
                              const nmod_mat_struct * E,
                              slong d)
{
    if (d <= PMINTBASIS_THRES)
    {
        nmod_poly_mat_mintbasis(intbas, shift, pts, E, d);
        return;
    }

    const slong m = intbas->r;

    const slong n = E[0].c;
    const slong d1 = (d + 1) / 2;
    const slong d2 = d - d1;

    /* No truncated copy of E needed for the first half -- E may be longer
       than d1, extra entries simply unused. */
    nmod_poly_mat_t P1;
    nmod_poly_mat_init(P1, m, m, intbas->modulus);
    nmod_poly_mat_pmintbasis(P1, shift, pts, E, d1);

   /* R_i = P1(pts[d1+i]) * E_{d1+i}, i = 0..d2-1 -- unlike pmbasis's
       residual (a middle product), this is d2 independent evaluate-then-
       multiply operations; no polynomial-arithmetic trick applies here
       (see this file's header comment / nmod_poly_mat_interpolant.h). */
    nmod_mat_struct * R = (nmod_mat_struct *) flint_malloc(FLINT_MAX(d2, 1) * sizeof(nmod_mat_struct));
    nmod_mat_t evalP1;
    nmod_mat_init(evalP1, m, m, intbas->modulus);
    for (slong idx = 0; idx < d2; idx++)
    {
        nmod_poly_mat_evaluate_nmod(evalP1, P1, pts[d1 + idx]);
        nmod_mat_init(R + idx, m, n, intbas->modulus);
        nmod_mat_mul(R + idx, evalP1, E + (d1 + idx));
    }
    nmod_mat_clear(evalP1);

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, m, m, intbas->modulus);
    nmod_poly_mat_pmintbasis(P2, shift, pts + d1, R, d2);
    for (slong idx = 0; idx < d2; idx++)
        nmod_mat_clear(R + idx);
    flint_free(R);

    nmod_poly_mat_multiply(intbas, P2, P1);

    nmod_poly_mat_clear(P1);
    nmod_poly_mat_clear(P2);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

