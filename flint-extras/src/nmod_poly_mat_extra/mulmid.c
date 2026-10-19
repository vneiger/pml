/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/longlong.h>   /* for flint_ctz */
#include <flint/mpn_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_multiply.h"
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
    nmod_poly_mat_shift_right(res, res, nlo);
    nmod_poly_mat_truncate(res, nhi - nlo);
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
        const int single_prime = (modbits <= 50)
            && ((flint_bitcnt_t) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) (len1 + len2 + 2)));
        /* 0, 1: a single transform, p of at most 32 bits or more; several
           transforms and 2: p of at most PML_NMOD_MAT_MUL_SMALL_BITS bits,
           3: at most 52 bits, 4: more */
        const int kind = single_prime ? (modbits <= 32 ? 0 : 1)
                                      : (modbits <= PML_NMOD_MAT_MUL_SMALL_BITS ? 2
                                         : (modbits <= 52 ? 3 : 4));
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

void nmod_poly_mat_mulmid(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                          slong nlo, slong nhi)
{
    PML_ASSERT(nlo >= 0);
    PML_ASSERT(nhi >= 0);

    slong len1 = nmod_poly_mat_max_length(pmat1);
    slong len2 = nmod_poly_mat_max_length(pmat2);
    nhi = FLINT_MIN(nhi, len1 + len2 - 1);

    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (len1 > nhi)
    {
        nmod_poly_mat_t pmat1_tmp;
        nmod_poly_mat_init(pmat1_tmp, pmat1->r, pmat1->c, pmat1->modulus);
        nmod_poly_mat_set_trunc(pmat1_tmp, pmat1, nhi);
        nmod_poly_mat_mulmid(res, pmat1_tmp, pmat2, nlo, nhi);
        nmod_poly_mat_clear(pmat1_tmp);
        return;
    }

    if (len2 > nhi)
    {
        nmod_poly_mat_t pmat2_tmp;
        nmod_poly_mat_init(pmat2_tmp, pmat2->r, pmat2->c, pmat2->modulus);
        nmod_poly_mat_set_trunc(pmat2_tmp, pmat2, nhi);
        nmod_poly_mat_mulmid(res, pmat1, pmat2_tmp, nlo, nhi);
        nmod_poly_mat_clear(pmat2_tmp);
        return;
    }

    /* len1 <= len2 and len(pmat1) <= nlo: shift away useless coefficients of pmat2 */
    if (len1 <= len2 && len1 <= nlo)
    {
        nmod_poly_mat_t pmat2_tmp;
        nmod_poly_mat_init(pmat2_tmp, pmat2->r, pmat2->c, pmat2->modulus);
        nmod_poly_mat_shift_right(pmat2_tmp, pmat2, nlo - len1 + 1);
        _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2_tmp, len2 - nlo + len1 - 1, len1 - 1, nhi - nlo + len1 - 1);
        nmod_poly_mat_clear(pmat2_tmp);
        return;
    }

    /* if len(pmat2) <= nlo (implies len1 > len2), shift away useless coefficients of pmat1 */
    if (len2 <= nlo)
    {
        nmod_poly_mat_t pmat1_tmp;
        nmod_poly_mat_init(pmat1_tmp, pmat1->r, pmat1->c, pmat1->modulus);
        nmod_poly_mat_shift_right(pmat1_tmp, pmat1, nlo - len2 + 1);
        _nmod_poly_mat_mulmid(res, pmat1_tmp, len1 - nlo + len2 - 1, pmat2, len2, len2 - 1, nhi - nlo + len2 - 1);
        nmod_poly_mat_clear(pmat1_tmp);
        return;
    }

    _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2, len2, nlo, nhi);
}
