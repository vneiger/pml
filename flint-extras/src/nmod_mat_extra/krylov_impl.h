/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_KRYLOV_IMPL_H
#define NMOD_MAT_KRYLOV_IMPL_H

#include <flint/nmod_mat.h>

#include "nmod_mat_extra.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Internal interface of the Krylov functions (krylov_*.c).

    Notation: E is m x n with rows E_0, ..., E_{m-1}, M is n x n, s is a
    shift (m integers), d is a list of degree bounds (m nonnegative
    integers). The Krylov matrix K_{s,d}(E, M) is formed by the rows
    E_i M^j for 0 <= j <= d_i, sorted by increasing (s_i + j, i). The pair
    (i, j) designates the row E_i M^j.

    The multiplication matrix M is described by an "operator": either a
    dense matrix, or a Jordan matrix given by the eigenvalues and sizes of
    its blocks. A Jordan block of size t with eigenvalue x is x I + Z with
    Z the upper shift matrix (ones on the superdiagonal), so that on the
    row vector v of length t, (v (x I + Z))_c = x v_c + v_{c-1} (with
    v_{-1} = 0).
*/

#define _NMOD_MAT_KRYLOV_DENSE 0     /* dense matrix */
#define _NMOD_MAT_KRYLOV_SHIFT 1     /* Jordan, all eigenvalues zero */
#define _NMOD_MAT_KRYLOV_DIAGONAL 2  /* Jordan, all blocks of size one */
#define _NMOD_MAT_KRYLOV_JORDAN 3    /* Jordan, general */

typedef struct
{
    int kind;
    slong n;                    /* dimension */
    nmod_t mod;
    const nmod_mat_struct * M;  /* dense kind: the matrix */
    slong logpow;               /* dense kind: powers M^(2^t), t <= logpow,
                                   may be used to generate rows (see
                                   _nmod_mat_krylov_basis_iter); 0: M only,
                                   negative: automatic choice (default) */
    slong powruns;              /* chunks with more than powruns indices
                                   only use M (0: no limit); ignored by
                                   the automatic choice */
    slong nblocks;              /* Jordan kinds: number of blocks */
    nn_ptr eig;                 /* eigenvalues (reduced), NULL for SHIFT */
    nn_ptr eig_pre;             /* Shoup precomputations, or NULL */
    slong * off;                /* block offsets, nblocks + 1 entries */
}
_nmod_mat_krylov_op_struct;

typedef _nmod_mat_krylov_op_struct _nmod_mat_krylov_op_t[1];

/* the operator only references M (dense case); eig and sizes are copied;
   eig == NULL means all eigenvalues zero, sizes == NULL all sizes one */
void _nmod_mat_krylov_op_init_dense(_nmod_mat_krylov_op_t op,
        const nmod_mat_t M);
void _nmod_mat_krylov_op_init_jordan(_nmod_mat_krylov_op_t op,
        nn_srcptr eig, const slong * sizes, slong nblocks, nmod_t mod);
void _nmod_mat_krylov_op_clear(_nmod_mat_krylov_op_t op);

/* out = in * M, single row, out and in must not overlap (Jordan kinds) */
void _nmod_mat_krylov_op_apply_vec(nn_ptr out, nn_srcptr in,
        const _nmod_mat_krylov_op_t op);

/*
    Binary min-heap of row indices i, ordered by (key[i], i); h has room
    for all indices, *len is the current number of elements.
*/
FLINT_FORCE_INLINE int
_nmod_mat_krylov_less(const slong * key, slong a, slong b)
{
    return key[a] < key[b] || (key[a] == key[b] && a < b);
}

FLINT_FORCE_INLINE void
_nmod_mat_krylov_heap_push(slong * h, slong * len, slong i, const slong * key)
{
    slong c = (*len)++, p;

    while (c > 0)
    {
        p = (c - 1) / 2;
        if (!_nmod_mat_krylov_less(key, i, h[p]))
            break;
        h[c] = h[p];
        c = p;
    }

    h[c] = i;
}

FLINT_FORCE_INLINE slong
_nmod_mat_krylov_heap_pop(slong * h, slong * len, const slong * key)
{
    slong top = h[0], last, c = 0, ch, n;

    n = --(*len);
    last = h[n];

    while ((ch = 2 * c + 1) < n)
    {
        if (ch + 1 < n && _nmod_mat_krylov_less(key, h[ch + 1], h[ch]))
            ch++;
        if (!_nmod_mat_krylov_less(key, h[ch], last))
            break;
        h[c] = h[ch];
        c = ch;
    }

    if (n > 0)
        h[c] = last;

    return top;
}

/*
    Shifts equivalent to shifts (NULL: zero shifts) for degree bounds at
    most D: writes in t (m entries) shifts in [0, (m - 1)(D + 1)] such that
    the order (t_i + j, i) for 0 <= j <= D coincides with the order
    (s_i + j, i), by replacing the gaps larger than D + 1 between the
    sorted shift values by D + 1. This avoids overflows with arbitrary
    slong shifts.
*/
void _nmod_mat_krylov_normalize_shifts(slong * t, const slong * shifts,
        slong m, slong D);

/*
    Fills d (m entries) with the degree bounds: degrees[i], or n if
    degrees is NULL; if clamp is nonzero, bounds larger than n are
    replaced by n. Throws if a bound is negative.
*/
void _nmod_mat_krylov_degrees(slong * d, slong m, const slong * degrees,
        slong n, int clamp);

/*
    Priority order: writes in (rows[k], exps[k]), k < N = sum (d_i + 1),
    the pairs (i, j) with 0 <= j <= d_i sorted by increasing (s_i + j, i);
    returns N. shifts == NULL means zero shifts. The arrays must have
    room for N entries.
*/
slong _nmod_mat_krylov_order(slong * rows, slong * exps, slong m,
        const slong * shifts, const slong * degrees);

/* number of rows sum (d_i + 1) */
slong _nmod_mat_krylov_nrows(slong m, const slong * degrees);

/*
    Builds K_{s,d}(E, M) in K (N x n, N = sum (d_i + 1)); if rows and exps
    are not NULL, they receive the pair (i, j) of each row of K.
*/
void _nmod_mat_krylov_matrix(nmod_mat_t K, slong * rows, slong * exps,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees);

/*
    Reduced row echelon form of a k x n matrix X, processing its rows in
    order (prime modulus), with pivots searched in the first lim columns
    only (the other columns are only carried along, for instance to record
    the row operations). Row a of the input is independent if its first lim
    entries are not in the span of those of rows 0, ..., a-1, which is
    recorded in indep[a]. On output, the first k' rows of X form a reduced
    row echelon basis: row t has a one in column Q[t] < lim, and the other
    k' - 1 rows are zero in that column, with Q[0] < Q[1] < ... (the pivot
    Q[t] is that of one of the independent rows, the first nonzero entry
    of this row once reduced by the rows above it); each of these rows is a
    linear combination of the independent input rows only. The remaining
    rows of X are destroyed. Returns k', the rank of the first lim columns.
    X can be a window. Based on nmod_mat_ple (rank profile matrix of the
    first lim columns) and two triangular solves.
*/
slong _nmod_mat_rref_rows(nmod_mat_t X, slong lim, slong * Q, char * indep);

/*
    Row rank profile of X: indep[a] = 1 if row a is not in the span of
    rows 0, ..., a-1, and 0 otherwise. Returns the rank; X is destroyed.
*/
slong _nmod_mat_row_rank_profile(char * indep, nmod_mat_t X);

/*
    Computes the row rank profile of K_{s,d}(E, M) assuming the degree
    bounds are valid (that is, for each i, E_i M^{d_i} is in the span of
    the rows which precede it in K_{s,d}): returns its rank r, writes the
    corresponding rows (the Krylov basis, in increasing priority) in the
    first r rows of B (which must have at least min(n, N) rows and n
    columns, where N = sum (d_i + 1)) and their pairs (i, j) in rows and
    exps (room for r entries).

    The following outputs are optional (NULL to ignore):
      - deltas[i] (length m): the minimal degree delta_i, that is,
        1 + the largest j such that (i, j) is in the profile, or 0;
      - T (m x n): row i is E_i M^{delta_i} (the target rows);
      - G (n x n): its leading r x r submatrix is set to the inverse of
        B restricted to the columns piv (requires piv), so that any row
        vector v in the row space of B is v[piv] G B;
      - piv (length n): r columns such that B restricted to these columns
        is invertible;
      - ok: set to 1 if the degree bounds were verified to be valid, to 0
        otherwise (some row (i, d_i) is in the profile; then the profile
        computed may differ from that of K_{s,d}, and deltas[i] = d_i + 1
        and row i of T is undefined for such i).

    Degree bounds larger than n are replaced by n; degrees == NULL means
    n for all rows, shifts == NULL means zero shifts. The parameter chunk
    is the maximal number of rows processed at once (0: default).

    For a dense operator, a run of c rows (i, j0), ..., (i, j0 + c - 1) of
    one index in a chunk is generated from the rows at distance
    e = min(2^logpow, largest power of two <= c), as (i, j) = (i, j - e) M^e,
    with the powers M^e computed by squaring when first needed (as in
    Keller-Gehrig's algorithm); only M and the highest power computed so
    far are kept, and the runs which would need a lower power use M.
    logpow = 0 means that each row is the previous one times M. This pays
    off when the products of a few rows by M are limited by the memory
    bandwidth (large n) and few indices are active (Hermite-type shifts),
    which the automatic choice targets.
*/
slong _nmod_mat_krylov_basis_iter(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t G, slong * piv, int * ok,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees, slong chunk);

/*
    Same output (B, rows, exps, deltas, ok) by the Keller-Gehrig-type
    doubling algorithm of Jeannerod, Neiger, Schost, Villard (JSC 2017,
    Algorithm KrylovRankProfile), relying on the same assumption on the
    degree bounds.
*/
slong _nmod_mat_krylov_basis_doubling(nmod_mat_t B, slong * rows,
        slong * exps, slong * deltas, int * ok,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees);

/*
    Row rank profile of the full matrix K_{s,d}(E, M), built explicitly
    (no assumption on the degree bounds; degrees == NULL means n).
*/
slong _nmod_mat_krylov_basis_naive(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees);

/*
    Kernel of the Krylov matrix in compact form: on success (return value
    r >= 0), writes the profile (rows, exps, r entries), the minimal
    degrees deltas (m entries) and the m x r relation matrix in the first
    r columns of R (m x n): row i of R contains the coefficients of
    E_i M^{delta_i} on the rows of the Krylov basis. Returns -1 if the
    degree bounds are found insufficient (some row (i, d_i) in the
    profile), in which case the outputs are undefined.
*/
slong _nmod_mat_krylov_kernel(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees, slong chunk);

#ifdef __cplusplus
}
#endif

#endif
