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

    Notation: V is m x D with rows V_0, ..., V_{m-1}, M is D x D, s is a
    shift (m integers), beta is a list of degree bounds (m nonnegative
    integers). The Krylov matrix K_{s,beta}(V, M) is formed by the rows
    V_i M^e for 0 <= e <= beta_i, sorted by increasing (s_i + e, i). The pair
    (i, e) designates the row V_i M^e.

    The multiplication matrix M is described by an "operator": either a
    dense matrix, or a Jordan matrix given by the eigenvalues and sizes of
    its blocks. The Jordan block b, of size D_b and eigenvalue x_b, is
    x_b I + Z with Z the upper shift matrix (ones on the superdiagonal), so
    that on the row vector v of length D_b, (v (x_b I + Z))_c = x_b v_c +
    v_{c-1} (with v_{-1} = 0).
*/

#define _NMOD_MAT_KRYLOV_DENSE 0     /* dense matrix */
#define _NMOD_MAT_KRYLOV_SHIFT 1     /* Jordan, all eigenvalues zero */
#define _NMOD_MAT_KRYLOV_DIAGONAL 2  /* Jordan, all blocks of size one */
#define _NMOD_MAT_KRYLOV_JORDAN 3    /* Jordan, general */

typedef struct
{
    int kind;
    slong D;                    /* dimension */
    nmod_t mod;
    const nmod_mat_struct * M;  /* dense kind: the matrix */
    slong logpow;               /* dense kind: powers M^(2^t), t <= logpow,
                                   may be used to generate rows (see
                                   _nmod_mat_krylov_basis_iter); 0: M only,
                                   negative: automatic choice (default) */
    slong powruns;              /* chunks with more than powruns indices
                                   only use M (0: no limit); ignored by
                                   the automatic choice */
    slong mpdeg;                /* degree of the minimal polynomial of M
                                   (Jordan kinds; D for the dense kind):
                                   degree bounds are capped at mpdeg */
    int predict;                /* budget of the indices in a chunk (see
                                   _nmod_mat_krylov_basis_iter): 1 (default)
                                   from the generic profile until an index
                                   stops before its predicted rows, 0 never
                                   (max(1, f_i) only), 2 always (tests) */
    slong nblocks;              /* Jordan kinds: number of blocks */
    nn_ptr eig;                 /* eigenvalues (reduced), NULL for SHIFT */
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
    slong top = h[0], last, c = 0, ch, hl;

    hl = --(*len);
    last = h[hl];

    while ((ch = 2 * c + 1) < hl)
    {
        if (ch + 1 < hl && _nmod_mat_krylov_less(key, h[ch + 1], h[ch]))
            ch++;
        if (!_nmod_mat_krylov_less(key, h[ch], last))
            break;
        h[c] = h[ch];
        c = ch;
    }

    if (hl > 0)
        h[c] = last;

    return top;
}

/*
    A shift equivalent to s = shift (NULL: zero shift) for degree bounds at
    most bmax: writes in t (m entries) a shift with entries in
    [0, (m - 1)(bmax + 1)] such that the order (t_i + e, i) for
    0 <= e <= bmax coincides with the order (s_i + e, i), by replacing the
    gaps larger than bmax + 1 between the sorted shift entries by bmax + 1.
    This avoids overflows with arbitrary slong shift entries.
*/
void _nmod_mat_krylov_normalize_shift(slong * t, const slong * shift,
        slong m, slong bmax);

/*
    Fills beta (m entries) with the degree bounds: degrees[i], or D if
    degrees is NULL; if clamp is nonzero, bounds larger than D are
    replaced by D. Throws if a bound is negative.
*/
void _nmod_mat_krylov_degrees(slong * beta, slong m, const slong * degrees,
        slong D, int clamp);

/*
    Priority order: writes in (rows[k], exps[k]), k < N = sum (beta_i + 1),
    the pairs (i, e) with 0 <= e <= beta_i sorted by increasing (s_i + e, i);
    returns N. shift == NULL means zero shift. The arrays must have
    room for N entries.
*/
slong _nmod_mat_krylov_order(slong * rows, slong * exps, slong m,
        const slong * shift, const slong * degrees);

/* number of rows sum (beta_i + 1) */
slong _nmod_mat_krylov_nrows(slong m, const slong * degrees);

/*
    Builds K_{s,beta}(V, M) in K (N x D, N = sum (beta_i + 1)); if rows and exps
    are not NULL, they receive the pair (i, e) of each row of K.
*/
void _nmod_mat_krylov_matrix(nmod_mat_t K, slong * rows, slong * exps,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees);

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
    Computes the row rank profile of K_{s,beta}(V, M) assuming the degree
    bounds are valid (that is, for each i, V_i M^{beta_i} is in the span of
    the rows which precede it in K_{s,beta}): returns its rank r, writes the
    corresponding rows (the Krylov basis, in increasing priority) in the
    first r rows of B (which must have at least min(D, N) rows and D
    columns, where N = sum (beta_i + 1)) and their pairs (i, e) in rows and
    exps (room for r entries).

    The following outputs are optional (NULL to ignore):
      - deltas[i] (length m): the minimal degree delta_i, that is,
        1 + the largest e such that (i, e) is in the profile, or 0;
      - T (m x D): row i is V_i M^{delta_i} (the target rows);
      - G (D x D): its leading r x r submatrix is set to the inverse of
        B restricted to the columns piv (requires piv), so that any row
        vector v in the row space of B is v[piv] G B;
      - piv (length D): r columns such that B restricted to these columns
        is invertible;
      - ok: set to 1 if the degree bounds were verified to be valid, to 0
        otherwise (some row (i, beta_i) is in the profile; then the profile
        computed may differ from that of K_{s,beta}, and deltas[i] = beta_i + 1
        and row i of T is undefined for such i).

    Degree bounds larger than D are replaced by D, and by op->mpdeg if
    larger (this does not change the profile); degrees == NULL means D for
    all rows, shift == NULL means zero shift. The parameter chunk is the
    maximal number of rows processed at once (0: default). Index i takes at
    most max(1, f_i) rows in a chunk, f_i being its number of rows already
    in the profile, or max(1, f_i, p_i - f_i) as long as the prediction
    holds (op->predict) and f_i >= 2^logpow, p_i being its number of rows
    (i, e), e < beta_i, among the first D rows in priority order: these rows
    form the profile for generic inputs.

    For a dense operator, a run of c rows (i, e0), ..., (i, e0 + c - 1) of
    one index in a chunk is generated from the rows at distance
    q = min(2^logpow, largest power of two <= min(c, e0)) (or 1 if e0 = 0),
    as (i, e) = (i, e - q) M^q,
    with the powers M^q computed by squaring when first needed (as in
    Keller-Gehrig's algorithm); only M and the highest power computed so
    far are kept, and the runs which would need a lower power use M.
    logpow = 0 means that each row is the previous one times M. This pays
    off when the products of a few rows by M are limited by the memory
    bandwidth (large D) and few indices are active (Hermite-type shift),
    which the automatic choice targets.
*/
slong _nmod_mat_krylov_basis_iter(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t G, slong * piv, int * ok,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk);

/*
    Same, with the basis kept in row echelon form instead of reduced row
    echelon form (each chunk is reduced by a triangular solve and one
    product, then decomposed by nmod_mat_ple): no update of the previous
    rows, and a factorization of the pivot block instead of G. On output
    (r = the returned rank), with C = B[:, piv] (r x r): the leading r x r
    block of LU holds U (on and above the diagonal, invertible) and L
    (strictly below, unit lower triangular), such that row t of L U is the
    row bpos[t] of C (bpos: room for D entries). Then the relations of a
    row vector v in the row space of B are v = x B with
    x[bpos[t]] = (v[piv] U^-1 L^-1)[t]. LU can be NULL (then bpos is
    ignored), for the basis only.
*/
slong _nmod_mat_krylov_basis_iter_lu(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t LU, slong * bpos, slong * piv,
        int * ok, const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk);

/*
    Same output (B, rows, exps, deltas, ok) by the Keller-Gehrig-type
    doubling algorithm of Jeannerod, Neiger, Schost, Villard (JSC 2017,
    Algorithm KrylovRankProfile), relying on the same assumption on the
    degree bounds.
*/
slong _nmod_mat_krylov_basis_doubling(nmod_mat_t B, slong * rows,
        slong * exps, slong * deltas, int * ok,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees);

/*
    Row rank profile of the full matrix K_{s,beta}(V, M), built explicitly
    (no assumption on the degree bounds; degrees == NULL means D).
*/
slong _nmod_mat_krylov_basis_naive(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees);

/*
    Kernel of the Krylov matrix in compact form: on success (return value
    r >= 0), writes the profile (rows, exps, r entries), the minimal
    degrees deltas (m entries) and the m x r relation matrix in the first
    r columns of R (m x D): row i of R contains the coefficients of
    V_i M^{delta_i} on the rows of the Krylov basis. Returns -1 if the
    degree bounds are found insufficient (some row (i, beta_i) in the
    profile), in which case the outputs are undefined.

    The relations R = T[:, piv] C^{-1}, with C = B[:, piv] the r x r
    invertible submatrix of the basis B and T the target rows, are
    obtained by one of two methods (the output is the same):
      - _NMOD_MAT_KRYLOV_KERNEL_TRACK: C^{-1} is maintained along the
        elimination (_nmod_mat_krylov_basis_iter with G), then one product
        T[:, piv] C^{-1}; extra cost about r^3 whatever m is;
      - _NMOD_MAT_KRYLOV_KERNEL_LU: the basis is kept in row echelon form
        (_nmod_mat_krylov_basis_iter_lu), which gives C = Pi^T L U along
        the elimination, then two triangular solves with m right-hand
        sides (nmod_mat_solve_triu_right, nmod_mat_solve_tril_right);
      - _NMOD_MAT_KRYLOV_KERNEL_AUTO: LU unless m is large compared with D
        (see krylov_kernel.c).
*/
#define _NMOD_MAT_KRYLOV_KERNEL_AUTO 0
#define _NMOD_MAT_KRYLOV_KERNEL_TRACK 1
#define _NMOD_MAT_KRYLOV_KERNEL_LU 2

slong _nmod_mat_krylov_kernel(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk, int method);

#ifdef __cplusplus
}
#endif

#endif
