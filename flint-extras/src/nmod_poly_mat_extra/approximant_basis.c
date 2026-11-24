/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly_mat_multiply.h"  // for middle product
#include "nmod_poly_mat_approximant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "impl.h"

/* nmod_poly_mat_mbasis: see nmod_mat_poly_extra/nmod_mat_poly_mbasis.c */

void _nmod_poly_mat_approximant_basis_krylov_leaf(nmod_poly_mat_t P,
                                                  slong * shift,
                                                  const nmod_poly_mat_t F,
                                                  const slong * orders)
{
    nmod_poly_mat_approximant_basis_krylov(P, F, orders, shift);

    /* P is shift-Popov: the pivot of row i is P[i,i] */
    for (slong i = 0; i < P->r; i++)
        shift[i] += nmod_poly_degree(nmod_poly_mat_entry(P, i, i));
}

void _nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                            slong * shift,
                            const nmod_poly_mat_t ipmat,
                            slong order,
                            double krylov_cutoff,
                            slong thres)
{
    FLINT_ASSERT(order >= 0);

    /* no condition: the identity, shift unchanged (the recursion below never
       creates such a node, see the assertion on thres) */
    if (order <= 0 || ipmat->c == 0)
    {
        nmod_poly_mat_one(appbas);
        return;
    }

    FLINT_ASSERT(thres >= 1);

    /* Krylov leaf: replaces the subtree below this node */
    if (_nmod_poly_mat_krylov_wanted(ipmat->r, ipmat->c, ipmat->c * order, krylov_cutoff))
    {
        slong * orders = flint_malloc(ipmat->c * sizeof(slong));
        for (slong j = 0; j < ipmat->c; j++)
            orders[j] = order;
        _nmod_poly_mat_approximant_basis_krylov_leaf(appbas, shift, ipmat, orders);
        flint_free(orders);
        return;
    }

    nmod_poly_mat_t pmat;

    nmod_poly_mat_init(pmat, ipmat->r, ipmat->c, ipmat->modulus);
    nmod_poly_mat_set_trunc(pmat,ipmat,order);

    if (order <= thres)
    {
        nmod_poly_mat_mbasis(appbas, shift, pmat, order);
        nmod_poly_mat_clear(pmat);
        return;
    }

    const long order1 = order>>1;
    const long order2 = order - order1;
    nmod_poly_mat_t appbas2, residual;

    nmod_poly_mat_init(appbas2, pmat->r, pmat->r, pmat->modulus);
    nmod_poly_mat_init(residual, pmat->r, pmat->c, pmat->modulus);

    _nmod_poly_mat_pmbasis(appbas, shift, pmat, order1, krylov_cutoff, thres);

    nmod_poly_mat_mulmid(residual, appbas, pmat, order1, order);

    _nmod_poly_mat_pmbasis(appbas2, shift, residual, order2, krylov_cutoff, thres);

    nmod_poly_mat_multiply(appbas, appbas2, appbas);

    nmod_poly_mat_clear(appbas2);
    nmod_poly_mat_clear(residual);
    nmod_poly_mat_clear(pmat);
}

void nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           slong order)
{
    _nmod_poly_mat_pmbasis(appbas, shift, pmat, order,
                           _nmod_poly_mat_pmbasis_krylov_cutoff(pmat->r),
                           _nmod_poly_mat_pmbasis_threshold(pmat->r, pmat->c, pmat->modulus));
}
