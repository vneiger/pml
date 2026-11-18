/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

/*
    With A = X[:, 0:lim] and its PLE decomposition (nmod_mat_ple, which
    reveals the rank profile matrix of A): the rows of A in its row rank
    profile are the rows P[0], ..., P[r-1], and the pivot of row P[t] is
    column c_t = Q[t], with c_0 < ... < c_{r-1}; that is, c_t is the first
    nonzero entry of row P[t] once reduced by the rows above it, as in a
    row by row elimination. With L1 (r x r, unit lower triangular) and E
    (r x lim, echelon) from A[P[0..r), :] = L1 E, and U = E[:, Q[0..r)]
    (upper triangular, invertible), the reduced echelon form of these rows
    is Y = U^-1 L1^-1 X[P[0..r), :]: on the first lim columns, the identity
    on the pivot columns and U^-1 E on the others; on the carried columns,
    U^-1 L1^-1 applied to the carried part of these rows (one triangular
    solve with each factor). The rows of Y are in the order of their pivot
    columns.
*/
slong
_nmod_mat_rref_rows(nmod_mat_t X, slong lim, slong * Q, char * indep)
{
    slong k = X->r, n = X->c, ncar = n - lim, nz, r, t, j;
    slong * P, * Qc;
    nmod_mat_t A, Z, U, L1, Zc;

    if (k == 0)
        return 0;

    if (lim == 0)
    {
        memset(indep, 0, k);
        return 0;
    }

    P = flint_malloc((k + lim) * sizeof(slong));
    Qc = P + k;

    /* the PLE decomposition of the first lim columns, in place; the rows
       of the carried columns are not moved */
    nmod_mat_window_init(A, X, 0, 0, k, lim);
    r = nmod_mat_ple(P, Qc, A, 0);
    nmod_mat_window_clear(A);

    memset(indep, 0, k);
    for (t = 0; t < r; t++)
        indep[P[t]] = 1;

    if (r == 0)
    {
        flint_free(P);
        return 0;
    }

    /* Z = [E2 | carried columns of the rows P[0..r)], where E2 holds the
       entries of E in the columns without pivot (zero before the pivot of
       each row, where A holds L) */
    nz = lim - r;
    nmod_mat_init(Z, r, nz + ncar, X->mod.n);

    for (t = 0; t < r; t++)
    {
        nn_srcptr x = nmod_mat_entry_ptr(X, t, 0);
        nn_ptr z = nmod_mat_entry_ptr(Z, t, 0);

        for (j = 0; j < nz; j++)
            z[j] = (Qc[r + j] > Qc[t]) ? x[Qc[r + j]] : 0;
        if (ncar > 0)
            memcpy(z + nz, nmod_mat_entry_ptr(X, P[t], lim), ncar * sizeof(ulong));
    }

    /* carried columns: L1^-1 */
    if (ncar > 0 && r > 1)
    {
        nmod_mat_window_init(L1, X, 0, 0, r, r);
        nmod_mat_window_init(Zc, Z, 0, nz, r, nz + ncar);
        nmod_mat_solve_tril(Zc, L1, Zc, 1);
        nmod_mat_window_clear(Zc);
        nmod_mat_window_clear(L1);
    }

    /* U^-1, with U read on and above the diagonal */
    for (j = 0; j < r && Qc[j] == j; j++)
        ;

    if (j == r)
    {
        nmod_mat_window_init(U, X, 0, 0, r, r);
        nmod_mat_solve_triu(Z, U, Z, 0);
        nmod_mat_window_clear(U);
    }
    else
    {
        nmod_mat_init(U, r, r, X->mod.n);
        for (t = 0; t < r; t++)
            for (j = t; j < r; j++)
                nmod_mat_entry(U, t, j) = nmod_mat_entry(X, t, Qc[j]);
        nmod_mat_solve_triu(Z, U, Z, 0);
        nmod_mat_clear(U);
    }

    /* the reduced rows */
    for (t = 0; t < r; t++)
    {
        nn_ptr x = nmod_mat_entry_ptr(X, t, 0);
        nn_srcptr z = nmod_mat_entry_ptr(Z, t, 0);

        for (j = 0; j < r; j++)
            x[Qc[j]] = (j == t);
        for (j = 0; j < nz; j++)
            x[Qc[r + j]] = z[j];
        if (ncar > 0)
            memcpy(x + lim, z + nz, ncar * sizeof(ulong));

        Q[t] = Qc[t];
    }

    nmod_mat_clear(Z);
    flint_free(P);

    return r;
}

slong
_nmod_mat_row_rank_profile(char * indep, nmod_mat_t X)
{
    slong k = X->r, r, t;
    slong * P;

    if (k == 0)
        return 0;

    P = flint_malloc(k * sizeof(slong));
    r = nmod_mat_ple(P, NULL, X, 0);

    memset(indep, 0, k);
    for (t = 0; t < r; t++)
        indep[P[t]] = 1;

    flint_free(P);
    return r;
}
