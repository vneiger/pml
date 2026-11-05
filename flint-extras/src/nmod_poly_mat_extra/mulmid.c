/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/mpn_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_utils.h"
#include "impl.h"                 /* for PML_NMOD_MAT_MUL_IFMA */

/*------------------------------------------------------------*/
/* naive: multiply, truncate, shift                           */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid_naive(nmod_poly_mat_t res,
                                 const nmod_poly_mat_t pmat1, slong len1,
                                 const nmod_poly_mat_t pmat2, slong len2,
                                 slong nlo, slong nhi)
{
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    nmod_poly_mat_multiply(res, pmat1, pmat2);
    nmod_poly_mat_shift_right_trunc(res, res, nlo, nhi - nlo);
}



/*------------------------------------------------------------*/
/* general interfaces                                         */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                           const nmod_poly_mat_t pmat1, slong len1,
                           const nmod_poly_mat_t pmat2, slong len2,
                           slong nlo, slong nhi)
{
    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    /* TODO handle constant pmat1 or pmat2 properly */

#if PML_HAVE_MACHINE_VECTORS
    /*
        2026-10-01 Thresholds fitted on square balanced middle products
        (nlo = nhi/2, len1 = nlo + 1, len2 = nhi, the shape of the
        approximant basis algorithms; nhi a power of two, and 1.5 times
        one) on Zen 4, one thread, against FLINT-dev with the new
        nmod_mat_mul kernels, over a 21- and a 50-bit FFT prime, a 30-bit
        and a 60-bit prime. Against the previous rule (geometric from
        dimension 8 on, naive below), on the same measurements: average
        slowdown with respect to the best routine 2.29 -> 1.04, 95th
        percentile 8.0 -> 1.3.

        2026-10-15 Refitted on Zen 4 against FLINT main (#2842 merged,
        SIMD transpose #2871) and the 2026-10-15 product thresholds, same
        moduli and shapes. naive is the product followed by a truncation,
        and the product of short operands now goes to the Vandermonde
        routines, which made naive the best choice for short windows:
        modulo a 30-bit prime (np == 2) up to nhi = min(dim, 128), a 1.5-5x
        gain over the previous rule there. Over 1096 points: average
        slowdown 1.13 -> 1.016, 95th percentile 2.44 -> 1.15, worst
        11.1 -> 1.46.

        2026-10-16 thresholds fitted on Zen 4 with the transposed Vandermonde
        routines (same shapes; a 21- and a 50-bit FFT prime, a 30-bit and
        a generic 50-bit prime, a 60-bit prime; 1760 points). They replace
        naive and geometric: on short windows they are faster than both
        everywhere except at dimension 2, where all are within a few
        microseconds. Average slowdown with respect to the best routine,
        against the previous rule (which was fitted on a Xeon): 1.025 ->
        1.014 for a single transform, 1.010 -> 1.004 with several
        transforms and p of at most 52 bits, 1.058 -> 1.002 above 52 bits;
        95th percentile at most 1.13.

        2026-10-17 Refitted on Zen 4 with the transposed fft_small middle
        products and the corrected point-major layout of sd_fft_matmul
        (whose transforms had been slowed down severalfold by cache
        conflicts, at dimensions such as 16, 32 and 64). sd_fft_matmul is
        now the fastest fft_small routine from dimension 16-32, against
        96-128 before, and it takes over a good part of the windows of
        the Vandermonde routines. Same moduli; shapes lo=8 hi=16, lo=6
        hi=12 and lo=4 hi=12 of p-mulmid_tune (the last one is the
        residual of a generic pmbasis after the shift of
        nmod_poly_mat_mulmid); budget 5 s; the 60-bit data stops at
        dimension 512 of the first shape, and was completed by a separate
        run (budget 15 s) from which the rule above 52 bits was adjusted.
        Over 2074 points, slowdown with respect to the best routine,
        previous rule -> new rule: average 1.144 -> 1.006 for a single
        transform, 1.027 -> 1.010 for p <= 40 bits, 1.021 -> 1.006 for
        p <= 52 bits, 1.075 -> 1.005 above 52 bits; 95th percentile
        1.60 -> 1.05; worst 2.55 -> 1.42.

        2026-10-19 Refitted on an Ice Lake server (Xeon Gold 6354, with
        AVX512-IFMA, budget 10 s) together with the Zen 4 data above, and
        separately on an Apple M4 Max (budget 5 s); seven moduli (21- and
        50-bit FFT primes, 30-, 40-, 50-, 52- and 60-bit primes), shapes
        lo/hi = 8/16, 6/12, 4/12 and 12/16 of p-mulmid_tune. The two
        families differ (see impl.h): without the IFMA kernel of
        nmod_mat_mul, the Vandermonde routines and sd_fft_matmul are slower
        relative to sd_fft_direct (on the M4, sd_fft_direct wins up to
        dimension 32-48 except for small single-transform moduli, and
        vandermonde2 covers much longer windows with several transforms),
        and the table of one family loses 7-9% on average on the other.
        geometric is the fastest at 4 points out of about 4400 on the M4,
        by at most 8%, and not elsewhere. Slowdown with respect to the best
        routine, previous rule -> new rule: Ice Lake average 1.020 ->
        1.010 (95th percentile 1.14 -> 1.07), M4 1.066 -> 1.013 (1.51 ->
        1.10), Zen 4 1.007 -> 1.010.

        The rule: the Vandermonde routines up to nhi = 16 at dimension at
        most 4, and above up to nhi = min(max(lo, a2 * dim / 2), cap),
        vandermonde1 up to vdm1, vandermonde2 above (they need a prime
        modulus, and naive takes their place otherwise up to nhi = 16);
        the windowed fft_small routines everything else, sd_fft_matmul from
        dimension mm_dim except for long windows (nhi above mm_nhi) below
        dimension mm_long_dim, where its grouping under the memory bound
        costs more than it saves, sd_fft_direct otherwise; see the tables
        below. They transform at length about nhi (a power of two N
        fitting the window, or the transposed product when it is shorter,
        see impl.h).
    */
    {
        const ulong modn = pmat1->modulus;
        const slong dim = n_cbrt(pmat1->r * pmat1->c * pmat2->c);
        const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);
        const int single_prime = _pml_fft_small_single_transform(modn,
                                        FLINT_BIT_COUNT((ulong) (len1 + len2 + 2)),
                                        (ulong) pmat1->c * FLINT_MIN(len1, len2));
        /* 0, 1: a single transform, p of at most 32 bits or more; several
           transforms and 2: p of at most PML_NMOD_MAT_MUL_SMALL_BITS bits,
           3: at most 52 bits, 4: more */
        const int kind = single_prime ? (modbits <= 32 ? 0 : 1)
                                      : (modbits <= PML_NMOD_MAT_MUL_SMALL_BITS ? 2
                                         : (modbits <= 52 ? 3 : 4));

        /* thin shapes, i.e. a row vector times a matrix (pmat1->r == 1) or
           an inner dimension 1: the thresholds below were fitted on square
           shapes, and the Vandermonde routines have little to amortize
           their evaluations on, except for very few points. The rule:
           - vandermonde1 for short products, with at most vdm_npts points,
             of not too small matrices (dimension at least 6, e.g. 1 x 16
             times 16 x 16), where the plain product is dominated by the
             overhead of its many short polynomial products;
           - otherwise the product (FLINT's nmod_poly_mat_mul for such
             shapes, see nmod_poly_mat_multiply) followed by a truncation,
             when the shorter operand has length at most naive_len: its cost
             is about len1 * len2 per entry, against about the transform
             length of the longer one for sd_fft_direct;
           - sd_fft_direct otherwise.
           A matrix times a column vector (pmat2->c == 1) keeps the rule for
           square shapes: there, the Vandermonde routines are faster than
           naive modulo a 30-bit prime.

           2026-11 Measured with mulmid_bench.sh on Zen 4 and Ice Lake (with
           AVX512-IFMA) and Apple M4, 30-, 50- and 60-bit primes, 1 x m times
           m x m for m = 2..128, operands of length 4..64, low and balanced
           windows, against the rule for square shapes:
           - with IFMA, vandermonde1 is the fastest up to about 7 points,
             and for the 50-bit prime (several transforms, kind 3) up to
             about 31 points; beyond, it is up to 2-3x slower than naive or
             sd_fft_direct; without IFMA (M4), it only helps for very few
             points (3 for p above PML_NMOD_MAT_MUL_SMALL_BITS bits);
           - naive is the fastest for balanced operands up to length 16,
             except with IFMA for p above PML_NMOD_MAT_MUL_SMALL_BITS bits,
             where sd_fft_direct is faster from length 16 (by 1.2-1.7x for
             m >= 16); for unbalanced operands (e.g. lengths 4 and 32 or 64,
             as in the low product of Dixon's lifting), naive was 2-4x faster
             than all the others on a Xeon with IFMA;
           - gains against the rule for square shapes: up to 5-20x for
             m <= 4, typically 1.5-2.7x elsewhere.
           Exception: with IFMA and the 50-bit prime (kind 3), for unbalanced
           operands (the shorter at most a quarter of the longer) from
           dimension 6, the rule for square shapes is kept: there, its
           Vandermonde routines were 1.3-1.6x faster than both naive and
           sd_fft_direct on Ice Lake (low products of Dixon's lifting with the
           default order, lengths deg and 8 deg, deg = 4..32, m = 16..128). */
        if (pmat1->r == 1 || pmat1->c == 1)
        {
            const slong npts = FLINT_MIN(len1, len2)
                               + FLINT_MIN(nhi, len1 + len2 - 1) - nlo - 1;
#if PML_NMOD_MAT_MUL_IFMA
            const slong vdm_npts = (kind == 3) ? 31 : 7;
            const slong naive_len = (modbits <= PML_NMOD_MAT_MUL_SMALL_BITS) ? 16 : 8;
            const int square_rule = (kind == 3 && dim >= 6
                && 4 * FLINT_MIN(len1, len2) <= FLINT_MAX(len1, len2));
#else
            const slong vdm_npts = (modbits <= PML_NMOD_MAT_MUL_SMALL_BITS) ? 7 : 3;
            const slong naive_len = 16;
            const int square_rule = 0;
#endif
            if (! square_rule)
            {
                if (dim >= 6 && npts <= vdm_npts && NMOD_POLY_CAN_USE_VANDERMONDE1(modn, npts)
                    && n_is_prime(modn))
                    _nmod_poly_mat_mulmid_vandermonde1(res, pmat1, len1, pmat2, len2, nlo, nhi);
                else if (FLINT_MIN(len1, len2) <= naive_len)
                    _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi);
                else
                    _nmod_poly_mat_mulmid_sd_fft_direct(res, pmat1, len1, pmat2, len2, nlo, nhi);
                return;
            }
        }
        /* the Vandermonde routines up to nhi = min(max(lo, a2 * dim / 2), cap),
           vandermonde1 up to vdm1; sd_fft_matmul from dimension mm_dim, for
           nhi at most mm_nhi unless the dimension is at least mm_long_dim */
        typedef struct { slong lo, a2, cap, vdm1, mm_dim, mm_nhi, mm_long_dim; } _thr_t;
#if PML_NMOD_MAT_MUL_IFMA
        static const _thr_t thr[5] = {
            { 24,  4,   24, 24, 16, WORD_MAX,  0 },
            {  0,  4,   24, 24, 16, WORD_MAX,  0 },
            {  0, 16,  384, 32, 24, WORD_MAX,  0 },
            {  0, 32, 1024, 32, 16,     3072, 64 },
            {  0,  4,   12,  9, 16,     3072, 64 } };
#else
        static const _thr_t thr[5] = {
            {  0,  4,   24, 12, 12, WORD_MAX,   0 },
            {  0,  1,   12,  9, 48, WORD_MAX,   0 },
            { 48, 16, 1024, 16, 48,     1024, 128 },
            {  0,  8,  384,  9, 96, WORD_MAX,   0 },
            { 48,  4,  384,  9, 48,      768, 128 } };
#endif
        const _thr_t * t = thr + kind;
        const slong vdm_nhi = (dim <= 4) ? 16
                              : FLINT_MIN(FLINT_MAX(t->lo, t->a2 * dim / 2), t->cap);
        const slong vdm1_nhi = t->vdm1;

        /* the Vandermonde points must be distinct with invertible
           differences, hence the primality test, which comes last so as to
           be reached only for the windows where they are wanted */
        if (nhi <= vdm_nhi)
        {
            const slong npts = FLINT_MIN(len1, len2)
                               + FLINT_MIN(nhi, len1 + len2 - 1) - nlo - 1;
            if (nhi <= vdm1_nhi && NMOD_POLY_CAN_USE_VANDERMONDE1(modn, npts)
                && n_is_prime(modn))
            {
                _nmod_poly_mat_mulmid_vandermonde1(res, pmat1, len1, pmat2, len2, nlo, nhi);
                return;
            }
            if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, npts) && n_is_prime(modn))
            {
                _nmod_poly_mat_mulmid_vandermonde2(res, pmat1, len1, pmat2, len2, nlo, nhi);
                return;
            }
            /* composite modulus, or not enough points */
            if (nhi <= 16)
            {
                _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi);
                return;
            }
        }

        if (dim >= t->mm_dim && (nhi <= t->mm_nhi || dim >= t->mm_long_dim))
            _nmod_poly_mat_mulmid_sd_fft_matmul(res, pmat1, len1, pmat2, len2, nlo, nhi);
        else
            _nmod_poly_mat_mulmid_sd_fft_direct(res, pmat1, len1, pmat2, len2, nlo, nhi);
        return;
    }
#endif /* PML_HAVE_MACHINE_VECTORS */

    /* TODO rough thresholds, not finely tuned */
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    if (NMOD_POLY_CAN_USE_GEOMETRIC(pmat1->modulus, nhi)
        && ((pmat1->r >= 8 && pmat2->c >= 2) || (pmat1->r >= 2 && pmat2->c >= 8))
        && (len1 <= nlo+1 || len2 <= nlo+1)
        && n_is_prime(pmat1->modulus))
        _nmod_poly_mat_mulmid_geometric(res, pmat1, len1, pmat2, len2, nlo, nhi);

    else
#endif
        _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi);
}

/* view of the coefficients [k, h) of the entries of pmat (k < h), i.e.
 * view = (pmat mod x^h) div x^k, normalized, without copying coefficients:
 * the entries of view point into those of pmat, which must be left unchanged
 * while view is in use; view must be cleared by _mulmid_view_clear */
static void _mulmid_view_init(nmod_poly_mat_t view, const nmod_poly_mat_t pmat, slong k, slong h)
{
    nmod_poly_mat_init(view, pmat->r, pmat->c, pmat->modulus);
    for (slong i = 0; i < pmat->r; i++)
    {
        for (slong j = 0; j < pmat->c; j++)
        {
            const nmod_poly_struct * p = nmod_poly_mat_entry(pmat, i, j);
            slong len = FLINT_MIN(p->length, h) - k;
            while (len > 0 && p->coeffs[k + len - 1] == 0)
                len--;
            if (len > 0)
            {
                nmod_poly_struct * v = nmod_poly_mat_entry(view, i, j);
                v->coeffs = p->coeffs + k;
                v->alloc = len;
                v->length = len;
            }
        }
    }
}

static void _mulmid_view_clear(nmod_poly_mat_t view)
{
    for (slong i = 0; i < view->r; i++)
    {
        for (slong j = 0; j < view->c; j++)
        {
            nmod_poly_struct * v = nmod_poly_mat_entry(view, i, j);
            v->coeffs = NULL;
            v->alloc = 0;
            v->length = 0;
        }
    }
    nmod_poly_mat_clear(view);
}

void nmod_poly_mat_mulmid(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                          slong nlo, slong nhi)
{
    PML_ASSERT(nlo >= 0);
    PML_ASSERT(nhi >= 0);

    const slong len1 = nmod_poly_mat_max_length(pmat1);
    const slong len2 = nmod_poly_mat_max_length(pmat2);
    nhi = FLINT_MIN(nhi, len1 + len2 - 1);

    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    /* reductions, as views of the operands (no coefficient is copied):
     * - an operand longer than nhi is truncated at nhi;
     * - then, as in the shapes for which the routines were tuned, at most
     *   one operand drops its useless low coefficients: if len1 <= len2 and
     *   len1 <= nlo, the coefficients of pmat2 below k2 = nlo - len1 + 1
     *   do not contribute to [nlo, nhi); otherwise, if len2 <= nlo, those of
     *   pmat1 below k1 = nlo - len2 + 1 do not contribute.
     * Dropping the useless coefficients of both operands would be valid as
     * well, but it shortens the window without changing the number of
     * points of the evaluation-based routines, and the dispatch, which uses
     * nhi, then makes worse choices: e.g. 2x slower on Zen 4 for 64 x 64
     * matrices of length 32 and the window [48, 63), modulo a 30-bit prime.
     * One has k1 < h1 and k2 < h2. */
    const slong h1 = FLINT_MIN(len1, nhi);
    const slong h2 = FLINT_MIN(len2, nhi);
    slong k1 = 0, k2 = 0;
    if (h1 <= h2 && h1 <= nlo)
        k2 = nlo - h1 + 1;
    else if (h2 <= nlo)
        k1 = nlo - h2 + 1;
    const int reduce1 = (k1 > 0 || h1 < len1);
    const int reduce2 = (k2 > 0 || h2 < len2);

    if (!reduce1 && !reduce2)
    {
        _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2, len2, nlo, nhi);
        return;
    }

    /* the operands are replaced by views of their useful coefficients,
     * which the routines called below cannot detect as aliases of res */
    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, pmat1->r, pmat2->c, pmat1->modulus);
        nmod_poly_mat_mulmid(tmp, pmat1, pmat2, nlo, nhi);
        nmod_poly_mat_swap_entrywise(res, tmp);
        nmod_poly_mat_clear(tmp);
        return;
    }

    nmod_poly_mat_t view1, view2;
    if (reduce1)
        _mulmid_view_init(view1, pmat1, k1, h1);
    if (reduce2)
        _mulmid_view_init(view2, pmat2, k2, h2);

    _nmod_poly_mat_mulmid(res, reduce1 ? view1 : pmat1, h1 - k1,
                               reduce2 ? view2 : pmat2, h2 - k2,
                               nlo - k1 - k2, nhi - k1 - k2);

    if (reduce1)
        _mulmid_view_clear(view1);
    if (reduce2)
        _mulmid_view_clear(view2);
}

/*------------------------------------------------------------*/
/* middle products by a fixed right operand                   */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid_precomp_right_init(nmod_poly_mat_mulmid_precomp_right_t pre,
                                              const nmod_poly_mat_t pmat2,
                                              slong len1, slong nlo, slong nhi,
                                              int method)
{
    PML_ASSERT(nlo >= 0);

    const slong len2 = nmod_poly_mat_max_length(pmat2);

    pre->pmat2 = pmat2;
    pre->len1 = len1;
    pre->nlo = nlo;
    pre->nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    pre->shift = 0;
    pre->fft = NULL;

    if (method == 0 || len1 <= 0 || len2 == 0 || nlo >= pre->nhi
        || pmat2->r == 0 || pmat2->c == 0)
        return;

    nhi = pre->nhi;

    /* as in nmod_poly_mat_mulmid: pmat2 is truncated at nhi, and when the
     * entries of pmat1 are not longer than those of pmat2 nor than nlo, its
     * coefficients below nlo - len1 + 1 are dropped; the products then
     * compute the window [nlo - shift, nhi - shift) */
    const slong h1 = FLINT_MIN(len1, nhi);
    const slong h2 = FLINT_MIN(len2, nhi);
    const slong k2 = (h1 <= h2 && h1 <= nlo) ? nlo - h1 + 1 : 0;

    /* automatic choice: stored transforms whenever they are available,
     * within a memory bound of the larger of 1 GiB and 16 times the
     * coefficients of pmat2 (method 1: no bound)
     * TODO for short operands, products of short polynomials or the
     * Vandermonde routines are faster than the stored transforms for small
     * dimensions (see the rule of nmod_poly_mat_solve_series_dixon); for
     * now, the callers decide */
    const ulong maxbytes = (method == 1) ? 0
        : n_max(UWORD(1) << 30, 16 * (ulong) pmat2->r * pmat2->c * (h2 - k2) * sizeof(ulong));

    if (k2 > 0 || h2 < len2)
    {
        nmod_poly_mat_t view2;
        _mulmid_view_init(view2, pmat2, k2, h2);
        pre->fft = _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right_init(view2,
                                h1, h2 - k2, nlo - k2, nhi - k2, maxbytes);
        _mulmid_view_clear(view2);
    }
    else
        pre->fft = _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right_init(pmat2,
                                h1, h2, nlo, nhi, maxbytes);

    pre->len1 = h1;
    pre->shift = k2;
}

void nmod_poly_mat_mulmid_precomp_right_init(nmod_poly_mat_mulmid_precomp_right_t pre,
                                             const nmod_poly_mat_t pmat2,
                                             slong len1, slong nlo, slong nhi)
{
    _nmod_poly_mat_mulmid_precomp_right_init(pre, pmat2, len1, nlo, nhi, -1);
}

void nmod_poly_mat_mulmid_precomp_right(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1,
                                        const nmod_poly_mat_mulmid_precomp_right_t pre,
                                        slong nhi)
{
    nhi = FLINT_MIN(nhi, pre->nhi);

    if (pre->fft == NULL || nhi <= pre->nlo
        || nmod_poly_mat_max_length(pmat1) > pre->len1)
    {
        nmod_poly_mat_mulmid(res, pmat1, pre->pmat2, pre->nlo, nhi);
        return;
    }

    if (pmat1->r == 0 || pre->pmat2->c == 0)
        return;

    if (res == pmat1)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, pmat1->r, pre->pmat2->c, pmat1->modulus);
        _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right(tmp, pmat1, pre->fft, nhi - pre->shift);
        nmod_poly_mat_swap_entrywise(res, tmp);
        nmod_poly_mat_clear(tmp);
        return;
    }

    _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right(res, pmat1, pre->fft, nhi - pre->shift);
}

void nmod_poly_mat_mulmid_precomp_right_clear(nmod_poly_mat_mulmid_precomp_right_t pre)
{
    _nmod_poly_mat_mulmid_sd_fft_direct_precomp_right_clear(pre->fft);
    pre->fft = NULL;
}
