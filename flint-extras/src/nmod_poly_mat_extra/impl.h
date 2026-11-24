/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_EXTRA_IMPL_H
#define NMOD_POLY_MAT_EXTRA_IMPL_H

#include <flint/nmod_types.h>
#include <flint/nmod_poly.h>  /* for geometric_progression_t */

#include "pml.h"               /* for PML_HAVE_MACHINE_VECTORS */
#include "nmod_mat_poly.h"     /* for nmod_mat_poly_t */

/* Two families of thresholds for nmod_poly_mat_multiply and
   nmod_poly_mat_mulmid (mul.c, mulmid.c): whether nmod_mat_mul has FLINT's
   AVX512-IFMA kernel, which mirrors NMOD_MAT_HAVE_MUL_U52 in FLINT's
   src/nmod_mat/impl.h (PML is assumed to be compiled for the same ISA as
   FLINT). They were measured on machines with that kernel (Zen 4, Ice Lake)
   and without it (Apple M4): the cost of nmod_mat_mul -- what the
   Vandermonde routines and sd_fft_matmul spend their time in -- relative to
   the fft_small transforms differs by up to 2x between them. */
#if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512IFMA__)
# define PML_NMOD_MAT_MUL_IFMA 1
#else
# define PML_NMOD_MAT_MUL_IFMA 0
#endif

/* With several transforms, the moduli up to this size form a class of their
   own, for which nmod_mat_mul is fastest: up to 40 bits with the IFMA
   kernel, up to 32 bits otherwise (on the M4, a 40-bit prime is closer to a
   50-bit one) */
#if PML_NMOD_MAT_MUL_IFMA
# define PML_NMOD_MAT_MUL_SMALL_BITS 40
#else
# define PML_NMOD_MAT_MUL_SMALL_BITS 32
#endif

/* Middle products with a short left operand and few columns go to products
   of coefficient matrices (_nmod_poly_mat_mulmid_coeffs, mulmid.c), for
   pmat1 of size at least k x k, k >= 2, and pmat2 with n columns, when
   - (len1 n)^2 <= NMOD_POLY_MAT_MULMID_COEFFS_K k, for a modulus of at most
     NMOD_POLY_MAT_MULMID_COEFFS_BITS bits;
   - (len1 n)^2 <= NMOD_POLY_MAT_MULMID_COEFFS_K_LARGE k and
     len1 n <= NMOD_POLY_MAT_MULMID_COEFFS_MAX_LARGE, for larger moduli, where
     nmod_mat_mul has no SIMD kernel in FLINT and the products of constant
     matrices lose sooner against the FFT routines. */
#ifndef NMOD_POLY_MAT_MULMID_COEFFS_K
# define NMOD_POLY_MAT_MULMID_COEFFS_K 64
#endif

#ifndef NMOD_POLY_MAT_MULMID_COEFFS_BITS
# define NMOD_POLY_MAT_MULMID_COEFFS_BITS 52
#endif

#ifndef NMOD_POLY_MAT_MULMID_COEFFS_K_LARGE
# define NMOD_POLY_MAT_MULMID_COEFFS_K_LARGE 16
#endif

#ifndef NMOD_POLY_MAT_MULMID_COEFFS_MAX_LARGE
# define NMOD_POLY_MAT_MULMID_COEFFS_MAX_LARGE 32
#endif

/* Hermite form helpers */
void _atomic_solve_pivot_collision_uechelon_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                                    slong pi1, slong pi2, slong j);
void _pivot_collision_xgcd_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                                      slong pi, slong ii, slong j,
                                                      nmod_poly_t g, nmod_poly_t u, nmod_poly_t v,
                                                      nmod_poly_t pivg, nmod_poly_t nonzg);
void _reduce_against_pivot_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                            slong i, slong j, slong ii,
                                            nmod_poly_t u, nmod_poly_t v);
ulong _normalize_pivot_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j);
void _normalize_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong * pivind, slong rk);

/* weak Popov form helpers */
void _atomic_solve_pivot_collision_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                           slong pi1, slong pi2, slong j);
void _reduce_against_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                   slong i, slong j, slong ii,
                                   nmod_poly_t u, nmod_poly_t v);
ulong _normalize_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j);

/* comparator for qsort on slong's, increasing order */
static inline int _slong_comparator(const void * a, const void * b)
{
    const slong x = *(const slong *) a;
    const slong y = *(const slong *) b;
    return (x > y) - (x < y);
}

/* Views (nmod_poly_mat_utils.c): a view of a polynomial matrix pmat is a
 * polynomial matrix whose entries are shallow copies of (parts of) entries
 * of pmat. No coefficient is copied: pmat must be left unchanged while the
 * view is in use, and the view must be cleared by _nmod_poly_mat_view_clear
 * (not by nmod_poly_mat_clear). */

/* view = the submatrix of pmat with rows ridx[0..nr) and columns cidx[0..nc)
 * (ridx, resp. cidx, NULL: the first nr rows, resp. nc columns, in order) */
void _nmod_poly_mat_view_init_submatrix(nmod_poly_mat_t view, const nmod_poly_mat_t pmat,
                                        const slong * ridx, slong nr,
                                        const slong * cidx, slong nc);

/* view = (pmat mod x^h) div x^k (requires k < h), with normalized entries */
void _nmod_poly_mat_view_init_coeffs(nmod_poly_mat_t view, const nmod_poly_mat_t pmat,
                                     slong k, slong h);

void _nmod_poly_mat_view_clear(nmod_poly_mat_t view);

/* Evaluation at a geometric progression, by chunks of points
 * (nmod_poly_mat_utils.c), as in nmod_poly_mat_det_geometric and
 * nmod_poly_mat_linsolve_geometric. */

/* vals[(i * pmat->c + j) * npts + k] = pmat[i][j](z q^k), for all entries
 * (i, j) and 0 <= k < npts, where q = r^2 for the ratio r of G. G must be
 * prepared for evaluation, with G->len at least npts and at least the length
 * of the entries of pmat. Each entry f is evaluated as f(z x) at 1, q, ...,
 * q^(npts-1): its coefficients are scaled by the powers of z (no scaling if
 * z == 1). For the points q^k0, ..., q^(k0 + npts - 1), take z = q^k0. */
void _nmod_poly_mat_evaluate_geometric(nn_ptr vals, const nmod_poly_mat_t pmat, ulong z,
                                       const nmod_geometric_progression_t G, slong npts);

/* Number of points per chunk for the evaluation of pmat at npts points by
 * _nmod_poly_mat_evaluate_geometric: `chunk` if it is positive, and
 * otherwise such that a chunk holds at most max(r c len, 2^22) values, where
 * r x c are the dimensions of pmat and len the largest length of its entries
 * (longer chunks make fewer and longer products, which are more efficient);
 * in any case, at least len and at most npts (npts must be at least len). */
slong _nmod_poly_mat_evaluate_geometric_chunk_length(const nmod_poly_mat_t pmat,
                                                     slong npts, slong chunk);

/* For a square n x n matrix A and an n x k matrix B, modulo a prime: if A is
 * invertible, sets c = det(A) and X = A^{-1} B and returns 1; otherwise
 * returns 0 (c and X undefined). A is overwritten by its LU decomposition
 * (nmod_mat_lu). perm (n entries) and PB (n x k) are temporary storage
 * provided by the caller, so that repeated calls do not allocate. X may
 * alias B, but not A or PB. (nmod_poly_mat_det.c) */
int _nmod_mat_det_solve_lu(ulong * c, nmod_mat_t X, nmod_mat_t A, const nmod_mat_t B,
                           slong * perm, nmod_mat_t PB);

/* ------------------------------------------------------------------------ */
/* Krylov leaves of the approximant and interpolant bases                   */
/* ------------------------------------------------------------------------ */

/* The direct Krylov functions (nmod_poly_mat_approximant_basis_krylov,
 * nmod_poly_mat_interpolant_basis_krylov) cost about D^2 (D + m) operations
 * for an m-row problem with D conditions (the sum of the orders), whatever
 * the shift and however D is split into columns, whereas the cost of pmbasis
 * and pmintbasis grows with m and decreases with the number of columns n.
 * Krylov is chosen when
 *     n < m   and   D <= c m / n^(1/3),
 * with two values of the cutoff c:
 * - _nmod_poly_mat_pmbasis_krylov_cutoff(m), inside pmbasis and pmintbasis
 *   (one order d, D = n d): the test is made at each node of the recursion, so
 *   that Krylov replaces the whole subtree below the first node where it
 *   holds;
 * - NMOD_POLY_MAT_APPROXIMANT_BASIS_KRYLOV_CUTOFF, in the divide and conquer
 *   for orders per column (nmod_poly_mat_approximant_basis and
 *   nmod_poly_mat_interpolant_basis), where the alternative to a direct call
 *   is a split whose two halves have dense m x m products, and where nodes
 *   with one order go to pmbasis; there n is the average number of
 *   conditions per level, ceil(D / dmax) for orders d_b <= dmax (resp.
 *   ceil(D / number of points) for interpolants), which is n for one order
 *   and is smaller than the number of columns when few columns have the
 *   largest orders.
 * 2026-11 Values measured with p-approximant_basis (modes leaf, ileaf, orders,
 * interp) on Zen 4, Ice Lake, Cascade Lake and Apple M4, at 30, 50 and 60 bits
 * (m = 16 to 256, n = 1 to m/2, D = m/4 to 32m), with the previous mbasis and
 * mintbasis and the base case at order 32 (see the note below for the current
 * values).
 * Inside pmbasis and pmintbasis, c = 6 was the best of 4, 6, 8, 12, 16 on
 * all twelve machine and modulus pairs (1.02 to 1.04 of the best value on
 * average, 1.04 to 1.09 for 8). For the direct calls, c = 8 is the best of
 * 4, 8, 16 or within 0.03 of it; c = 4 would be better for m = 256 (about
 * 1.2x at D = 4m) and worse for m = 16.
 * TODO a cutoff depending on m for the direct calls; at the smallest orders
 * with n = m/4 to m/2 the leaf can lose (PR table, 2026-11, the four machines,
 * 30 to 60 bits: pmintbasis 64 x 16 at d = 8 down to 0.56x, 32 x 16 at d = 4
 * down to 0.61x, pmbasis 64 x 16 at d = 8 down to 0.72x, all below 2 ms),
 * but n < m/4 would lose more than it saves (at n = m/4 the leaf is 1.07 to
 * 1.52 times faster on average over all orders). */
/* FIXME NOTE With the rewritten mbasis and mintbasis (2026-11, tune, m = 2 to 64,
 * 31 and 60 bits; Zen 4, Comet Lake, Apple M4, Emerald Rapids VM), the leaf
 * of pmbasis and pmintbasis loses for small m: no leaf for m <= 4, c = 3 for
 * 4 < m < 64 and c = 6 from m = 64 on (the value above, measured up to
 * m = 256) are within 1.5% to 2.6% of the better choice on average on each
 * machine and for pmbasis, pmintbasis, pmintbasis_geometric, against 5% to
 * 13% for c = 6 everywhere.
 * When the shift-Popov form is requested (nmod_poly_mat_approximant_basis and
 * nmod_poly_mat_interpolant_basis with POPOV), a Krylov leaf at the top gives
 * it directly, whereas otherwise a second basis is computed: there the cutoff
 * is NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF_POPOV for all m. With the per-m
 * cutoff instead, POPOV was 1.25 to 1.85 times slower in that range (Sapphire
 * Rapids VM, 3 x 1 to 32 x 2, 31 and 60 bits, up to 2.5 ms); 4 or 5 may be
 * slightly better than 6 (TODO tune). */
#ifndef NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF
# define NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF 6.0
#endif
#ifndef NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF_POPOV
# define NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF_POPOV 6.0
#endif
#ifndef NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF_MEDIUM
# define NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF_MEDIUM 3.0
#endif

/* the cutoff of the Krylov leaf of pmbasis and pmintbasis for m rows */
static inline double
_nmod_poly_mat_pmbasis_krylov_cutoff(slong m)
{
    if (m <= 4)
        return 0.0;
    return (m < 64) ? NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF_MEDIUM
                    : NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF;
}

#ifndef NMOD_POLY_MAT_APPROXIMANT_BASIS_KRYLOV_CUTOFF
# define NMOD_POLY_MAT_APPROXIMANT_BASIS_KRYLOV_CUTOFF 8.0
#endif

/* whether a node with m rows, n columns (or conditions per point) and D
 * conditions in total goes to Krylov: n < m and D^3 n <= (cutoff m)^3; never
 * for n = 0, D = 0 or cutoff <= 0 */
static inline int
_nmod_poly_mat_krylov_wanted(slong m, slong n, slong D, double cutoff)
{
    double a;

    if (n < 1 || n >= m || D < 1 || cutoff <= 0.0)
        return 0;

    a = cutoff * (double) m;
    return (double) D * (double) D * (double) D * (double) n <= a * a * a;
}

/* Krylov leaves with the contract of pmbasis and pmintbasis: shift is the
 * input shift on entry and the shift-row degree of the output on exit; the
 * output is the shift-Popov basis, whose shift-row degree is shift[i] +
 * deg(P[i,i]) (approximant_basis.c, interpolant_basis.c). */

/* approximant basis of F at order orders[j] for column j */
void _nmod_poly_mat_approximant_basis_krylov_leaf(nmod_poly_mat_t P, slong * shift,
                                                  const nmod_poly_mat_t F,
                                                  const slong * orders);

/* interpolant basis of (V, pts), one point per column of V */
void _nmod_poly_mat_interpolant_basis_krylov_leaf(nmod_poly_mat_t P, slong * shift,
                                                  const nmod_mat_t V,
                                                  nn_srcptr pts);

/* interpolant basis of (E, pts) as in nmod_poly_mat_pmintbasis: d matrices
 * E[k] of size m x n at the points pts[0], ..., pts[d-1] */
void _nmod_poly_mat_pmintbasis_krylov_leaf(nmod_poly_mat_t P, slong * shift,
                                           const ulong * pts,
                                           const nmod_mat_struct * E, slong d);

/* R[:, b] = P(pts[b]) V[:, b] for each column b of V, from the coefficients
 * of P: R = sum_k P_k V diag(pts)^k by Horner's rule, deg(P) + 1 products of
 * constant matrices; with c columns, R is P->r x c and V is P->c x c, and R
 * must not alias V (interpolant_basis.c) */
void _nmod_poly_mat_evaluate_mul_columns(nmod_mat_t R, const nmod_poly_mat_t P,
                                         const nmod_mat_t V, nn_srcptr pts);

/* Whether the residual of an interpolant basis with m rows and n columns per
 * point is computed from the coefficients of the first basis
 * (_nmod_poly_mat_evaluate_mul_columns) rather than by one evaluation of it
 * per point: n^2 <= m and 16 n <= m for a modulus of at most
 * NMOD_POLY_MAT_MULMID_COEFFS_BITS bits (see interpolant_basis.c); above,
 * where nmod_mat_mul has no SIMD kernel in FLINT, n = 1 with m >= 32, or
 * n^2 <= m and 64 n <= m. Comes from measurements on 2026-11-16 on 4
 * Intel / AMD / Apple machines. On the residual step alone (2026-11-17, Ice
 * Lake, Cascade Lake, Apple M4; 30, 50, 60 bits; m = 16 to 256, n up to m/8,
 * deg(P1) = n d2 / m and 4 times smaller or larger), this rule takes 1.05 to
 * 1.16 times the time of the better of the two on average (at worst 2.4 to
 * 6.9): the better choice also depends on deg(P1), as in the rule below.
 * TODO a rule depending on deg(P1) */
static inline int
_nmod_poly_mat_evaluate_mul_columns_wanted(slong m, slong n, ulong modulus)
{
    if (FLINT_BIT_COUNT(modulus) <= NMOD_POLY_MAT_MULMID_COEFFS_BITS)
        return (n * n <= m && 16 * n <= m);
    return (n == 1) ? (m >= 32) : (n * n <= m && 64 * n <= m);
}

/* The same choice in pmintbasis_geometric, where the alternative is the
 * evaluation of all the entries of the first basis at the geometric points
 * (_nmod_poly_mat_evaluate_geometric), whose cost does not depend on its
 * degree deg: the coefficients when n^(3/2) (deg + 1) <= 2 m for a modulus of
 * at most NMOD_POLY_MAT_MULMID_COEFFS_BITS bits, and <= m above.
 *
 * Measured 2026-11, the residual step alone (Ice Lake, Cascade Lake, Apple
 * M4; 30, 50, 60 bits; m = 16 to 256, n up to m/8, d2 = 16 to 8192, deg = n d2 / m
 * and 4 times smaller or larger): this rule takes 1.00 to 1.06 times the time
 * of the better of the two on average (at worst 2.7); the ratio geometric /
 * coefficients fits m^(0.4 to 0.8) / (n^0.95 (deg + 1)^0.65) up to a constant.
 *
 * On the whole of pmintbasis_geometric (PR table, 2026-11-17, same machines,
 * m = 16 to 512), this takes 0.65 to 1.05 times the time of pmintbasis at
 * 30 and 50 bits (median over the orders of each m x n). At 50 bits on x86,
 * where the geometric evaluation takes twice its 30-bit time and the
 * products of coefficients the same, the geometric version is up to 1.8
 * times slower than pmintbasis at the largest orders with few columns (32 x 1
 * at d = 4096), and the constant 3 would be better there (residual step: 1.01
 * instead of 1.03 to 1.06 on average).
 * TODO a constant depending on the modulus size and the architecture
 **/
static inline int
_nmod_poly_mat_evaluate_mul_columns_wanted_geometric(slong m, slong n,
                                                     slong deg, ulong modulus)
{
    const double c = (FLINT_BIT_COUNT(modulus) <= NMOD_POLY_MAT_MULMID_COEFFS_BITS) ? 2.0 : 1.0;
    const double y = (double) (deg + 1);
    /* n^(3/2) (deg + 1) <= c m, squared */
    return (double) n * n * n * y * y <= c * c * (double) m * (double) m;
}

/* The order up to which pmbasis calls mbasis directly (its base case), for
 * an m x n input modulo `modulus`; the base cases then have orders from
 * about half of it. With the mbasis that stores the rows as vectors of
 * coefficients, the base case pays off up to about 512 conditions n*order:
 * the threshold is 512/n, at least 32 (it was 32), except 192 for n = 2 (and
 * 512 for n = 0).
 *
 * TODO provisional, to be re-tuned more thoroughly
 * Some notes: Measured on a Xeon (Emerald Rapids) VM, with the Krylov
 * leaf off, m = 2 to 64, n = 1 to 3m/4, orders 128 to 2500, moduli of 31
 * and 60 bits: 1.04 times the time with the best threshold among 16 to 512
 * on average, at most 1.45 (32: 1.58 on average, at most 3.7). Provisional,
 * to be measured on other machines; the best value also grows with the
 * modulus (384/n at 31 bits, 768/n at 60 bits). */
static inline slong
_nmod_poly_mat_pmbasis_threshold(slong m, slong n, ulong modulus)
{
    (void) m;
    (void) modulus;
    if (n <= 2)
        return (n <= 1) ? 512 : 192;
    return FLINT_MAX(32, 512 / n);
}

/* The number of points up to which pmintbasis calls mintbasis directly, for
 * m x n matrices E_k modulo `modulus`: 1024/n, at least 32 and at most 512,
 * except 128 for m <= 2 and 256 for m = 3 (same measurements as above: 1.04
 * times the best threshold on average, at most 1.6; 32: 1.59 on average, at
 * most 4.0). Provisional. */
static inline slong
_nmod_poly_mat_pmintbasis_threshold(slong m, slong n, ulong modulus)
{
    (void) modulus;
    if (m <= 3)
        return (m <= 2) ? 128 : 256;
    return FLINT_MAX(32, FLINT_MIN(512, 1024 / FLINT_MAX(n, 1)));
}

/* The same for pmintbasis_geometric, whose recursion costs less (the
 * residual by evaluation at geometric points): the threshold of pmbasis,
 * and at most 128 for m <= 4 (same measurements: 1.04 times the best
 * threshold on average, at most 1.35; 32: 1.46 on average, at most 3.8).
 * Provisional. */
static inline slong
_nmod_poly_mat_pmintbasis_geometric_threshold(slong m, slong n, ulong modulus)
{
    const slong t = _nmod_poly_mat_pmbasis_threshold(m, n, modulus);
    return (m <= 4) ? FLINT_MIN(t, 128) : t;
}

/* pmbasis, pmintbasis and pmintbasis_geometric with the Krylov cutoff
 * (krylov_cutoff = 0: no Krylov leaf) and the threshold of the base case as
 * parameters; nmod_poly_mat_pmbasis, nmod_poly_mat_pmintbasis and
 * nmod_poly_mat_pmintbasis_geometric use _nmod_poly_mat_pmbasis_krylov_cutoff
 * and the thresholds above. For profiling, tuning and tests. */
void _nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas, slong * shift,
                            const nmod_poly_mat_t pmat, slong order,
                            double krylov_cutoff, slong thres);

void _nmod_poly_mat_pmintbasis(nmod_poly_mat_t intbas, slong * shift,
                               const ulong * pts, const nmod_mat_struct * E,
                               slong d, double krylov_cutoff, slong thres);

void _nmod_poly_mat_pmintbasis_geometric(nmod_poly_mat_t intbas, slong * shift,
                                         ulong * pts, const nmod_mat_struct * E,
                                         ulong r, slong d,
                                         double krylov_cutoff, slong thres);

/* Divide and conquer for orders per column (approximant_basis_general.c)
 * and for one point per column (interpolant_basis_general.c), behind
 * nmod_poly_mat_approximant_basis and nmod_poly_mat_interpolant_basis. A
 * node with one order for all its columns (resp. as many columns at each of
 * its points) goes to pmbasis (resp. pmintbasis); otherwise a node with
 * few conditions goes to the direct Krylov function; otherwise, for
 * approximants whose orders are close to each other, the node goes to
 * pmbasis on its input with each column multiplied by x^(dmax - d_b), at
 * order dmax; otherwise the node is split into two halves of about D/2
 * conditions each. */

/* close orders: pmbasis on the input raised to the largest order when
 * n dmax <= (1 + eps) D, eps = NMOD_POLY_MAT_APPROXIMANT_BASIS_ORDER_SLACK
 * (never if eps is negative). 2026-11 measurements (Zen 4, Ice Lake, Cascade
 * Lake, Apple M4; 30, 50, 60 bits): 0.25 is the best of -1, 0.25, 0.5, 1 on
 * average; with two distinct orders (d, d, d/2, d/2: n dmax = 1.33 D), the
 * split gives two pmbasis calls and eps must be at most 0.25, while many
 * distinct orders (m/2 columns in [d/4, d]: 1.6 D) prefer eps = 1.
 * TODO a slack depending on the number of distinct orders */
#ifndef NMOD_POLY_MAT_APPROXIMANT_BASIS_ORDER_SLACK
# define NMOD_POLY_MAT_APPROXIMANT_BASIS_ORDER_SLACK 0.25
#endif

typedef struct
{
    double leaf_cutoff;    /* cutoff of the Krylov leaf of pmbasis, pmintbasis;
                            * negative: _nmod_poly_mat_pmbasis_krylov_cutoff(m) */
    double direct_cutoff;  /* cutoff of the direct Krylov calls */
    double order_slack;    /* close orders, see above (approximants) */
    int split;             /* approximants: 0 to cut at a quasi-uniform order
                            * tau (orders min(d_b, tau) in the first half),
                            * 1 to cut each order in half (interpolants always
                            * cut at D/2) */
}
_nmod_poly_mat_relation_basis_dac_params_struct;

typedef _nmod_poly_mat_relation_basis_dac_params_struct
        _nmod_poly_mat_relation_basis_dac_params_t[1];

/* the parameters used by nmod_poly_mat_approximant_basis and
 * nmod_poly_mat_interpolant_basis */
static inline void
_nmod_poly_mat_relation_basis_dac_params_init(_nmod_poly_mat_relation_basis_dac_params_t prm)
{
    prm->leaf_cutoff = -1.0;
    prm->direct_cutoff = NMOD_POLY_MAT_APPROXIMANT_BASIS_KRYLOV_CUTOFF;
    prm->order_slack = NMOD_POLY_MAT_APPROXIMANT_BASIS_ORDER_SLACK;
    prm->split = 0;
}

/* the cutoff of the Krylov leaf for m rows given by prm */
static inline double
_nmod_poly_mat_relation_basis_leaf_cutoff(const _nmod_poly_mat_relation_basis_dac_params_t prm,
                                          slong m)
{
    return (prm->leaf_cutoff < 0.0) ? _nmod_poly_mat_pmbasis_krylov_cutoff(m) : prm->leaf_cutoff;
}

/* shift-ordered weak Popov approximant basis of F at order orders[j] for
 * column j, with the contract of pmbasis (shift in, shifted row degree out;
 * shift entries plus the sum of the orders at most WORD_MAX) */
void _nmod_poly_mat_approximant_basis(nmod_poly_mat_t P, slong * shift,
                                      const nmod_poly_mat_t F,
                                      const slong * orders,
                                      const _nmod_poly_mat_relation_basis_dac_params_t prm);

/* shift-ordered weak Popov interpolant basis of (V, pts), one point per
 * column of V, with the contract of pmintbasis (shift in, shifted row degree
 * out); the columns of V must be sorted by point (equal points adjacent) */
void _nmod_poly_mat_interpolant_basis(nmod_poly_mat_t P, slong * shift,
                                      const nmod_mat_t V, nn_srcptr pts,
                                      const _nmod_poly_mat_relation_basis_dac_params_t prm);

/* the s-Popov basis from an s-ordered weak Popov basis P of the same module,
 * with shifted row degree rdeg (in place): with delta = rdeg - s, a basis R
 * computed by `basis` with shift -delta has (-delta)-row degree 0, and P =
 * lm(R)^(-1) R, with lm(R) its (-delta)-leading matrix; `basis` computes a
 * shift-ordered weak Popov basis of the module with the contract of pmbasis
 * (shift in, shifted row degree out), `data` is passed to it */
void _nmod_poly_mat_popov_from_second_relation_basis(nmod_poly_mat_t P,
        const slong * shift, const slong * rdeg,
        void (* basis)(nmod_poly_mat_t, slong *, const void *), const void * data);


/* ------------------------------------------------------------------------ */
/* a few functions to help tests/profiles be more complete                  */
/* ------------------------------------------------------------------------ */

/* Same as nmod_poly_mat_set_trunc_from_mat_poly, with explicit control over
 * the block kernel and the loop schedule of the underlying transposition.
 *
 * `kern` is one of NMOD_MAT_POLY_CONV_{SCALAR,VEC4,VEC8} (see
 * nmod_mat_poly_extra/impl.h); any other value selects the default, narrowed
 * if its block does not fit inside the problem.
 *
 * `dmaj` is 1 for the entry-major schedule (the stores into the output
 * polynomials are long sequential streams) and 0 for the coefficient-major
 * one (the loads from the input matrices are long sequential streams); any
 * other value (e.g. -1) selects the default.  Note that the two schedules
 * carry the opposite names in the other direction: `dmaj` always means
 * "sweep the destination rows sequentially". */
void _nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                            const nmod_mat_poly_t matp,
                                            slong order,
                                            int kern,
                                            int dmaj);

/* Same as _nmod_poly_mat_mul_geometric_precomp, with the soft bound (in
 * bytes) on the memory used for the constant matrices given explicitly;
 * 0 selects the default. Below the bound the product is performed in one
 * go, above it the rows of A and then the columns of B are processed by
 * groups. Exposed so that the tests can exercise the grouping at sizes
 * that do not need hundreds of megabytes. */
void _nmod_poly_mat_mul_geometric_precomp_bounded(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G,
                                          ulong membytes);

/* Same for the middle product _nmod_poly_mat_mulmid_geometric_precomp. */
void _nmod_poly_mat_mulmid_geometric_precomp_bounded(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi,
                                             nmod_geometric_progression_t G,
                                             ulong membytes);

/* Same as _nmod_poly_mat_mulmid_sd_fft_direct and
 * _nmod_poly_mat_mulmid_sd_fft_matmul, with the soft bound (in bytes) on the
 * memory used for the transforms given explicitly; 0 selects the default.
 * Exposed so that the tests can exercise the grouping of rows and columns. */
void _nmod_poly_mat_mulmid_sd_fft_direct_bounded(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1, slong len1,
                                        const nmod_poly_mat_t pmat2, slong len2,
                                        slong nlo, slong nhi, ulong membytes);
void _nmod_poly_mat_mulmid_sd_fft_matmul_bounded(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1, slong len1,
                                        const nmod_poly_mat_t pmat2, slong len2,
                                        slong nlo, slong nhi, ulong membytes);

/* Precomputed transforms of a right operand for the middle products of
 * mul_sd_fft_direct.c (see nmod_poly_mat_mulmid_precomp_right_init in
 * nmod_poly_mat_multiply.h): the plan for the window [nlo, nhi) of products
 * A * B with A having entries of length at most lenA and B (fixed) entries
 * of length at most lenB, and the transforms of B. Returns NULL if
 * PML_HAVE_MACHINE_VECTORS is not set, if the plan fails, if the arguments
 * are degenerate (B with no row or column, lenA or lenB not positive, empty
 * window, nhi > lenA + lenB - 1), or if the transforms would take more than
 * maxbytes bytes (0: no bound). B is not used after the call. The product
 * computes the window [nlo, nhi') of A * B, for nlo < nhi' <= nhi, into C
 * (A->r x B->c, not aliasing A); A must have B->r columns and entries of
 * length at most lenA. */
void * _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right_init(const nmod_poly_mat_t B,
                                                              slong lenA, slong lenB,
                                                              slong nlo, slong nhi,
                                                              ulong maxbytes);
void _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right(nmod_poly_mat_t C,
                                                       const nmod_poly_mat_t A,
                                                       const void * pre, slong nhi);
void _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right_clear(void * pre);

/* ------------------------------------------------------------------------ */
/* middle products of the fft_small based multiplications                   */
/* ------------------------------------------------------------------------ */

#if PML_HAVE_MACHINE_VECTORS

#include <flint/fft_small.h>

/*
    Whether the fft_small plans of the routines below use a single transform
    (np == 1) rather than several CRT primes, for transforms of depth about
    `depth` and a bound len_bound on the number of products accumulated onto
    one coefficient (k * min(lenA, lenB) for an inner dimension k). This
    follows the choice made in fft_small_plan_init_nmod:
    - the direct transform modulo p, which these plans force whenever p
      allows it (see _sd_fft_mulmid_plan_init): p of at least 20 bits and
      fft_small_prime_supports_depth (p of at most 50 bits satisfying the
      bounds of its modular multiplication, 2^max(depth, 4) dividing p - 1);
      primality is not checked, it is left to the plan;
    - otherwise, a single CRT prime of the context when it is large enough,
      at least len_bound * 2^(2 * bits(p)): typically for p below 20 bits
      and products that are not too long, where the direct transform has
      nothing to gain (a single transform is used anyway).
    Used by the dispatch of nmod_poly_mat_multiply and _nmod_poly_mat_mulmid.
*/
static inline int _pml_fft_small_single_transform(ulong modn, ulong depth,
                                                  ulong len_bound)
{
    const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);
    const mpn_ctx_struct * R;

    if (modbits >= 20 && fft_small_prime_supports_depth(modn, depth))
        return 1;

    R = get_default_mpn_ctx();
    return flint_mpn_cmp_ui_2exp(crt_data_prod_primes(R->crts),
                                 R->crts[0].coeff_len,
                                 len_bound, 2 * modbits) >= 0;
}

/*
    Middle products by transposition, shared by mul_sd_fft_direct.c and
    mul_sd_fft_matmul.c (only with PML_HAVE_MACHINE_VECTORS).

    The window [nlo, nhi) of the product of s, of length sn, and l is the
    transposed product (see _nmod_poly_mul_tr_mpn_ctx in FLINT's
    fft_small/nmod_poly_mul.c)

        z[j] = sum_{i < sn} rev_sn(s)[i] * g[i + j],     0 <= j < nhi - nlo,

    where g = x^(sn - 1 - nlo) * l, the negative powers dropped, is read on
    its first gn = nhi - nlo + sn - 1 coefficients. Indeed rev_sn(s)[i] =
    s[sn - 1 - i] and g[i + j] = l[nlo + j - (sn - 1 - i)], so that z[j] is
    the coefficient nlo + j of s*l. The transposed product costs a product
    of lengths sn and nhi - nlo: its transforms have length about gn, against
    about lenA + lenB - 1 (or a power of two N fitting the window, see
    mul_sd_fft_direct.c) for the window of the full product.

    For matrices, with s the side of the shorter entries (sn = min(lenA,
    lenB)), only the two ends of the pipeline change:

        - the entries of the s side are reversed to length sn and
          transformed forward, as in a product (fft_small_fft_nmod);
        - the entries of the other side are shifted into g and sent through
          the transpose of the inverse transform (fft_small_ifft_tr_nmod),
          at the full truncation ztrunc;
        - the pointwise matrix products, and the normalization factor folded
          into the transforms of B, are unchanged: by linearity the
          transposed products accumulate over the inner dimension exactly as
          the products do;
        - each entry of C is sent through the transpose of the forward
          transform (fft_small_fft_tr), truncated to the coefficients
          wanted, and exported from [0, nhi - nlo).

    The ranges agree with the ones the kernels expect: fft_small_ifft_tr_nmod
    produces entries in (-3n, 3n) as fft_small_fft_nmod does, and
    fft_small_fft_tr accepts entries in (-2n, 2n) as fft_small_ifft does.
    The number of primes is unchanged too: each coefficient of C is still a
    sum of at most k * sn products of two residues.

    The transposed product is used only when its transform is strictly
    shorter than the one of the window of the full product, as computed by
    the window selection of both plans (this is the rule of FLINT's
    _nmod_poly_mulmid_fft_small); at equal length the window is kept, since
    it spares the reversals and the shifts. It is never shorter for a full
    product (then gn >= lenA + lenB - 1).
*/

typedef struct
{
    int on;         /* whether the transposed product is used */
    int revA;       /* if on: 1 if A is the reversed side s, 0 if B is */
    slong sn;       /* the reversal length, min(lenA, lenB) */
    slong lead;     /* sn - 1 - nlo: g = x^lead * l (lead < 0: drop -lead) */
    slong gn;       /* nhi - nlo + sn - 1: the length of g that is read */
}
_sd_fft_tr_struct;

/*
    Initialize the plan of the window [nlo, nhi) of the product of entries
    of lengths at most lenA and lenB, or of the transposed product when its
    transform is strictly shorter, and fill tr accordingly. Returns 0 if
    fft_small_plan_init_nmod fails, which needs a bound
    len_bound = k * min(lenA, lenB) beyond the capacity of the eight primes.

    The direct single-prime transform modulo p is forced whenever p allows
    it (direct_len = UWORD_MAX): its setup cost, proportional to the
    transform length, is amortized over the m*k + k*n + m*n transforms of a
    matrix product, unlike in a single polynomial product.

    xtrunc_max, the largest truncation an operand is transformed with, only
    matters to the window plan for the wraparound test (operands must fit in
    the cyclic length). Below one block the operands are transformed at the
    full transform length (see _fft_small_op_trunc), so their length itself
    is what must fit; for a plain product this makes no difference and the
    historical value is kept. The transposed plan has no wraparound (its
    window starts at 0), and takes the value of _nmod_poly_mul_tr_mpn_ctx.
*/
static inline int _sd_fft_mulmid_plan_init(fft_small_plan_t P, _sd_fft_tr_struct * tr,
                                           mpn_ctx_struct * R,
                                           slong lenA, slong lenB, slong nlo, slong nhi,
                                           ulong len_bound, nmod_t mod)
{
    const ulong zn = (ulong) lenA + lenB - 1;
    const ulong lmax = n_max(lenA, lenB);
    const ulong xtrunc_max = (nlo == 0 || lmax >= BLK_SZ)
                             ? n_max(n_round_up(lenA, BLK_SZ), n_round_up(lenB, BLK_SZ))
                             : lmax;
    const slong sn = FLINT_MIN(lenA, lenB);
    const ulong mm = (ulong) (nhi - nlo);
    const ulong gn = (ulong) sn + mm - 1;
    const ulong xtrunc_tr = n_max(_fft_small_len_trunc((ulong) sn),
                                  _fft_small_len_trunc(mm));

    tr->on = 0;
    tr->revA = 0;
    tr->sn = sn;
    tr->lead = sn - 1 - nlo;
    tr->gn = (slong) gn;

    if (gn < zn)
    {
        fft_small_plan_t Pw, Pt;
        _fft_small_plan_set_window_nmod(Pw, nlo, nhi, zn, xtrunc_max);
        _fft_small_plan_set_window_nmod(Pt, 0, mm, gn, xtrunc_tr);
        tr->on = (Pt->ztrunc < Pw->ztrunc);
    }

    if (tr->on)
    {
        tr->revA = (lenA <= lenB);
        return fft_small_plan_init_nmod(P, R, 0, mm, gn, xtrunc_tr, len_bound,
                                        2 * NMOD_BITS(mod), mod, UWORD_MAX);
    }

    return fft_small_plan_init_nmod(P, R, nlo, nhi, zn, xtrunc_max, len_bound,
                                    2 * NMOD_BITS(mod), mod, UWORD_MAX);
}

/* number of words of the scratch the two transforms below need */
static inline ulong _sd_fft_tr_scratch_len(const _sd_fft_tr_struct * tr)
{
    return tr->on ? (ulong) tr->gn : 0;
}

/* forward transform of the entry pol (assumed nonzero) of the reversed
   side, reversed to length sn, into X; t is scratch of gn >= sn words */
static inline void _sd_fft_tr_fft_rev(fft_small_op_t X, ulong * t,
                                      const nmod_poly_struct * pol,
                                      const _sd_fft_tr_struct * tr,
                                      nmod_t mod, const fft_small_plan_t P)
{
    _nmod_poly_reverse(t, pol->coeffs, pol->length, tr->sn);
    fft_small_fft_nmod(X, t, tr->sn,
                       n_min(_fft_small_op_trunc((ulong) tr->sn, P), P->ztrunc),
                       mod, P);
}

/* whether the shift g of an entry of length len of the other side is zero:
   then its transform is zero, and nothing is computed */
static inline int _sd_fft_tr_shift_is_zero(slong len, const _sd_fft_tr_struct * tr)
{
    return (len == 0) || (tr->lead < 0 && len <= -tr->lead);
}

/* the transpose of the inverse transform of the shift g of the entry pol of
   the other side, assumed not to be zero, into X; t is scratch of gn
   words, used only when zeros have to be prepended */
static inline void _sd_fft_tr_ifft_tr_shift(fft_small_op_t X, ulong * t,
                                            const nmod_poly_struct * pol,
                                            const _sd_fft_tr_struct * tr,
                                            nmod_t mod, const fft_small_plan_t P)
{
    const slong len = pol->length;

    if (tr->lead <= 0)
    {
        const slong off = -tr->lead;
        fft_small_ifft_tr_nmod(X, pol->coeffs + off,
                               FLINT_MIN(len - off, tr->gn), mod, P);
    }
    else
    {
        /* gn - lead = nhi > 0 */
        const slong glen = FLINT_MIN(len, tr->gn - tr->lead);
        _nmod_vec_zero(t, tr->lead);
        _nmod_vec_set(t + tr->lead, pol->coeffs, glen);
        fft_small_ifft_tr_nmod(X, t, tr->lead + glen, mod, P);
    }
}

/* the coefficients [0, zn) of a transposed product from its pointwise
   product X (domain FFT_SMALL_OP_PRODUCT), 0 < zn <= nhi - nlo; z receives
   zn coefficients */
static inline void _sd_fft_tr_export(ulong * z, slong zn, fft_small_op_t X,
                                     nmod_t mod, const fft_small_plan_t P)
{
    fft_small_fft_tr(X, n_min(_fft_small_op_trunc((ulong) zn, P), P->ztrunc), P);
    fft_small_export_nmod_range(z, X, 0, zn, mod, P);
}

#endif  /* PML_HAVE_MACHINE_VECTORS */

#endif  /* NMOD_POLY_MAT_EXTRA_IMPL_H */
