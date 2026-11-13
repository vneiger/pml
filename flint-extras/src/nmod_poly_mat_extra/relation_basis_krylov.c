/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_poly_mat_relation_basis.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_interpolant.h"

/*
    The three functions below form the m x sigma matrix E whose block of
    columns j holds the first orders[j] coefficients of column j of F (or
    of its expansion at the point pts[j]), and compute the kernel of the
    Krylov matrix of E with respect to the Jordan matrix with blocks
    (pts[j], orders[j]) (JNSV 2017, Section 1.2): for such a block, the
    product p * e of a polynomial p by the coefficient vector e of a
    polynomial f is the coefficient vector of p(x + pts[j]) f mod x^orders[j].
*/

void
nmod_poly_mat_approximant_basis_krylov(nmod_poly_mat_t P,
        const nmod_poly_mat_t F, const slong * orders, const slong * shift)
{
    slong m = F->r, n = F->c, i, j, sigma, nb, off;
    slong * sizes;
    nmod_mat_t E;

    if (P->r != m || P->c != m)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_approximant_basis_krylov: P must be F->r x F->r\n");

    sigma = 0;
    for (j = 0; j < n; j++)
    {
        if (orders[j] < 0)
            flint_throw(FLINT_ERROR, "nmod_poly_mat_approximant_basis_krylov: negative order\n");
        sigma += orders[j];
    }

    sizes = flint_malloc((n + 1) * sizeof(slong));
    nmod_mat_init(E, m, sigma, nmod_poly_mat_modulus(F));

    nb = 0;
    off = 0;
    for (j = 0; j < n; j++)
    {
        if (orders[j] == 0)
            continue;

        for (i = 0; i < m; i++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(F, i, j);
            _nmod_vec_set(nmod_mat_entry_ptr(E, i, off), f->coeffs,
                          FLINT_MIN(f->length, orders[j]));
        }

        sizes[nb++] = orders[j];
        off += orders[j];
    }

    if (!nmod_poly_mat_krylov_kernel_basis_jordan(P, E, NULL, sizes, nb,
                                                  shift, NULL))
        flint_throw(FLINT_ERROR, "nmod_poly_mat_approximant_basis_krylov: unexpected failure\n");

    nmod_mat_clear(E);
    flint_free(sizes);
}

void
nmod_poly_mat_interpolant_basis_krylov(nmod_poly_mat_t P, const nmod_mat_t V,
        nn_srcptr pts, const slong * shift)
{
    if (P->r != V->r || P->c != V->r)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_interpolant_basis_krylov: P must be V->r x V->r\n");

    if (!nmod_poly_mat_krylov_kernel_basis_jordan(P, V, pts, NULL, V->c,
                                                  shift, NULL))
        flint_throw(FLINT_ERROR, "nmod_poly_mat_interpolant_basis_krylov: unexpected failure\n");
}

void
nmod_poly_mat_mpade_basis_krylov(nmod_poly_mat_t P, const nmod_poly_mat_t F,
        nn_srcptr pts, const slong * orders, const slong * shift)
{
    slong m = F->r, n = F->c, i, j, sigma, nb, off;
    slong * sizes;
    nn_ptr eig, tmp;
    nmod_mat_t E;
    nmod_t mod;

    if (P->r != m || P->c != m)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_mpade_basis_krylov: P must be F->r x F->r\n");

    nmod_init(&mod, nmod_poly_mat_modulus(F));

    sigma = 0;
    for (j = 0; j < n; j++)
    {
        if (orders[j] < 0)
            flint_throw(FLINT_ERROR, "nmod_poly_mat_mpade_basis_krylov: negative order\n");
        sigma += orders[j];
    }

    sizes = flint_malloc((n + 1) * sizeof(slong));
    eig = _nmod_vec_init(n + 1);
    nmod_mat_init(E, m, sigma, mod.n);

    nb = 0;
    off = 0;
    for (j = 0; j < n; j++)
    {
        if (orders[j] == 0)
            continue;

        for (i = 0; i < m; i++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(F, i, j);
            slong len = f->length;

            if (len == 0)
                continue;

            /* coefficients of f(x + pts[j]) mod x^orders[j] */
            tmp = _nmod_vec_init(len);
            _nmod_vec_set(tmp, f->coeffs, len);
            _nmod_poly_taylor_shift(tmp, n_mod2_preinv(pts[j], mod.n, mod.ninv), len, mod);
            _nmod_vec_set(nmod_mat_entry_ptr(E, i, off), tmp,
                          FLINT_MIN(len, orders[j]));
            _nmod_vec_clear(tmp);
        }

        eig[nb] = pts[j];
        sizes[nb++] = orders[j];
        off += orders[j];
    }

    if (!nmod_poly_mat_krylov_kernel_basis_jordan(P, E, eig, sizes, nb,
                                                  shift, NULL))
        flint_throw(FLINT_ERROR, "nmod_poly_mat_mpade_basis_krylov: unexpected failure\n");

    nmod_mat_clear(E);
    _nmod_vec_clear(eig);
    flint_free(sizes);
}
