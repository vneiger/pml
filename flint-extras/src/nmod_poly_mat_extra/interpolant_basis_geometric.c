/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_extra.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

/* Same D&C shape as nmod_poly_mat_pmintbasis (interpolant_basis.c), but
   points are rho_start * rho^k (k=0..d-1), rho = r^2, using one shared,
   already-built G (see nmod_poly_mat_pmintbasis_geometric below).
   rho_start is passed down as a single field element rather than an
   index, so no point array needs to be materialized except at the base
   case. */
static void _pmintbasis_geometric_rec(nmod_poly_mat_t intbas,
                                      slong * shift,
                                      const nmod_mat_struct * E,
                                      slong n,
                                      ulong r,
                                      ulong rho_start,
                                      slong d,
                                      nmod_t mod,
                                      const nmod_geometric_progression_t G)
{
    ulong rho = nmod_mul(r, r, mod);
    const slong m = intbas->r;

    if (d <= PMINTBASIS_THRES)
    {
        ulong * pts = (ulong *) flint_malloc(d * sizeof(ulong));
        ulong cur = rho_start;
        for (slong k = 0; k < d; k++)
        {
            pts[k] = cur;
            cur = nmod_mul(cur, rho, mod);
        }
        nmod_poly_mat_mintbasis(intbas, shift, pts, E, d);
        flint_free(pts);
        return;
    }

    const slong d1 = (d + 1) / 2;
    const slong d2 = d - d1;

    /* No truncated copy of E needed for the first half -- E may be longer
       than d1, extra entries simply unused. */
    nmod_poly_mat_t P1;
    nmod_poly_mat_init(P1, m, m, intbas->modulus);
    _pmintbasis_geometric_rec(P1, shift, E, n, r, rho_start, d1, mod, G);

     /* R_i = P1(zeta*rho^i) * E_{d1+i}, i=0..d2-1, zeta = rho_start*rho^d1
       (the absolute start of the second half's points). Evaluating P1's
       entries at zeta*rho^i is done by rescaling each entry's coefficient
       k by zeta^k (Q(x) := entry(zeta*x)) and then evaluating the rescaled
       polynomial at the *first* d2 points of the shared G -- no offset
       handling needed inside G itself. */
    ulong zeta = nmod_mul(rho_start, nmod_pow_ui(rho, d1, mod), mod);

    /* a row of an order-d1 interpolant basis can have degree up to d1
       (length up to d1+1): the trivial witness M(x)*e_i, where M is the
       degree-d1 polynomial vanishing at all d1 points, is always a valid
       (if non-minimal) solution, so minimality never needs to exceed this
       -- same role X^d1*e_i plays for approximants, just a different
       degree-d1 witness since X^d1 alone isn't in the interpolant module
       (evaluation doesn't vanish automatically the way truncation does). */
    ulong * zeta_pow = (ulong *) flint_malloc((d1 + 1) * sizeof(ulong));
    zeta_pow[0] = 1;
    for (slong k = 1; k <= d1; k++)
        zeta_pow[k] = nmod_mul(zeta_pow[k - 1], zeta, mod);

    ulong * scaled = (ulong *) flint_malloc((d1 + 1) * sizeof(ulong));
    ulong * vs = (ulong *) flint_malloc(d2 * sizeof(ulong));
    ulong * P1vals = (ulong *) flint_malloc(m * m * d2 * sizeof(ulong));

    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
        {
            const nmod_poly_struct * entry = nmod_poly_mat_entry(P1, i, j);
            slong len = nmod_poly_length(entry);
            if (len == 0)
                _nmod_vec_zero(vs, d2);
            else
            {
                for (slong k = 0; k < len; k++)
                    scaled[k] = nmod_mul(nmod_poly_get_coeff_ui(entry, k), zeta_pow[k], mod);
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(vs, scaled, len, G, d2, mod);
            }
            for (slong idx = 0; idx < d2; idx++)
                P1vals[(i * m + j) * d2 + idx] = vs[idx];
        }

    flint_free(scaled);
    flint_free(vs);
    flint_free(zeta_pow);

    /* R, the "new E" for the second half's recursive call, is a plain
       array too -- nmod_mat_mul writes each R_i directly. */
    nmod_mat_struct * R = (nmod_mat_struct *) flint_malloc(FLINT_MAX(d2, 1) * sizeof(nmod_mat_struct));
    nmod_mat_t evalP1;
    nmod_mat_init(evalP1, m, m, intbas->modulus);
    for (slong idx = 0; idx < d2; idx++)
    {
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                nmod_mat_entry(evalP1, i, j) = P1vals[(i * m + j) * d2 + idx];
        nmod_mat_init(R + idx, m, n, intbas->modulus);
        nmod_mat_mul(R + idx, evalP1, E + (d1 + idx));
    }
    nmod_mat_clear(evalP1);
    flint_free(P1vals);

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, m, m, intbas->modulus);
    _pmintbasis_geometric_rec(P2, shift, R, n, r, zeta, d2, mod, G);
    for (slong idx = 0; idx < d2; idx++)
        nmod_mat_clear(R + idx);
    flint_free(R);

    nmod_poly_mat_multiply(intbas, P2, P1);

    nmod_poly_mat_clear(P1);
    nmod_poly_mat_clear(P2);
}

void nmod_poly_mat_pmintbasis_geometric(nmod_poly_mat_t intbas,
                                        slong * shift,
                                        ulong * pts,
                                        const nmod_mat_struct * E,
                                        ulong r,
                                        slong d)
{
    /* d=0: no points, output is the identity, shift unchanged -- matching
       the same convention as every other function in this project TO SEE. */
    if (d == 0)
    {
        nmod_poly_mat_one(intbas);
        return;
    }

    if (r == 0)
    flint_throw(FLINT_ERROR,
                "Exception (nmod_poly_mat_pmintbasis_geometric). "
                "r must be nonzero (r=0 degenerates FLINT's geometric-"
                "progression setup, which computes 1/r unconditionally).\n");

    const slong n = E[0].c;

    nmod_t mod;
    nmod_init(&mod, intbas->modulus);

    nmod_geometric_progression_t G;
    _nmod_geometric_progression_init_function(G, r, d, mod, 1); /* 1 = evaluation only */

    _pmintbasis_geometric_rec(intbas, shift, E, n, r, 1, d, mod, G);

    nmod_geometric_progression_clear(G);

    if (pts)
    {
        ulong rho = nmod_mul(r, r, mod);
        ulong cur = 1;
        for (slong k = 0; k < d; k++)
        {
            pts[k] = cur;
            cur = nmod_mul(cur, rho, mod);
        }
    }
}

/** Targets the pairwise distinct points case.    
 * 
 * Tries `nmod_find_root` (`nmod_extra.h`) first, rather than going
 * straight to `n_primitive_root_prime`, the algorithm only needs an 
 * element of multiplicative order order at least 2d−1, so that 
 *  rho = r^2 has order at least d and the d points 
 * are distinct; 
 * TO SEE: nmod_find_root(2*d) is documented to return order at least 2d, 
 * which suffices. 
 * 
 * Passing `2*d-1` would already suffice; `2*d` is used both for one unit
 * of margin and because `nmod_find_root` then returns 0 exactly when
 * `p <= 2*d+1`, that is, exactly when the precondition `p > 2*d+1` fails 
 *  (one unit of margin is kept too). 
 * `2*d`also matches the 2*len convention used by mul_geometric.c and mulmid.c.
 * */
void nmod_poly_mat_pmintbasis_geometric_auto(nmod_poly_mat_t intbas,
                                             slong * shift,
                                             ulong * pts,
                                             const nmod_mat_struct * E,
                                             slong d)
{
    /** With the best possible r (a primitive root), rho = r^2 has 
     * order (p-1)/2, so d distinct points exist iff (p-1)/2 >= d, 
     * i.e. p >= 2*d+1; we require the documented p > 2*d+1.  */
    if (intbas->modulus <= (ulong) (2 * d + 1))
        flint_throw(FLINT_ERROR,
                    "Exception (nmod_poly_mat_pmintbasis_geometric_auto). "
                    "Modulus %wu too small for d = %wd points "
                    "(requires p > 2*d+1 = %wd).\n",
                    intbas->modulus, d, 2 * d + 1);

    nmod_t mod;
    nmod_init(&mod, intbas->modulus);

    ulong r = nmod_find_root(2 * d, mod);
    
    // if (r == 0)                     /* not reachable under the check above:
    //                                    find_root(2*d) returns 0 iff p <= 2*d+1.
    //                                    Kept as a potential guard in case the check or the
    //                                    argument is ever changed independently. */
    //     r = n_primitive_root_prime(mod.n);

    nmod_poly_mat_pmintbasis_geometric(intbas, shift, pts, E, r, d);
}

