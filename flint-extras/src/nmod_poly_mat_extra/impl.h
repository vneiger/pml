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

/* ------------------------------------------------------------------------ */
/* middle products of the fft_small based multiplications                   */
/* ------------------------------------------------------------------------ */

#if PML_HAVE_MACHINE_VECTORS

#include <flint/fft_small.h>

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
