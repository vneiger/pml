/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/longlong.h>   /* for flint_ctz */
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_extra.h"  /* for NMOD_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_multiply.h"
#include "impl.h"                 /* for PML_NMOD_MAT_MUL_IFMA */

void nmod_poly_mat_multiply(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2)
{
    slong len1 = nmod_poly_mat_max_length(pmat1);
    slong len2 = nmod_poly_mat_max_length(pmat2);

    /* TODO once available, call multiply with constant */
    /* TODO once available, call mat-vec | vec-mat multiply */

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (pmat1->r == 1 || pmat2->c == 1  /* vec-mat or mat-vec */
        || pmat1->c == 1                /* independent products */
        || len1 == 1 || len2 == 1)      /* constant lhs or rhs */
    {
        nmod_poly_mat_mul(res, pmat1, pmat2);
        return;
    }

    /* FIXME decide where to handle this: currently done at multiply levels... */
    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, pmat1->r, pmat2->c, pmat1->modulus);
        nmod_poly_mat_multiply(tmp, pmat1, pmat2);
        nmod_poly_mat_swap_entrywise(res, tmp);
        nmod_poly_mat_clear(tmp);
        return;
    }

    /* from here, all dimensions and lengths are >= 2 */
    /* TODO thresholds essentially tuned from square cases
     * -> would be safe to check that they are ok for fat vectors */
    /* TODO see impact of FLINT+BLAS on thresholds */

    /*
        TODO 2026-09-13
        NOTE : measurements were focused on "small" matrices (but products
        already taking up to about 60s), with dimension <= 512, on square
        shapes only.
    */

    const slong dim = n_cbrt(pmat1->r * pmat1->c * pmat2->c);
    const slong len = len1 + len2 - 1;
    const ulong modn = pmat1->modulus;

#if PML_HAVE_MACHINE_VECTORS
    /*
        Six routines compete. Evaluation-interpolation at the roots of
        unity of fft_small, with the pointwise stage either a kernel on the
        transforms (nmod_poly_mat_mul_sd_fft_direct) or nmod_mat_mul at each
        point (nmod_poly_mat_mul_sd_fft_matmul); at small points with the
        evaluations themselves done by nmod_mat_mul against a Vandermonde
        matrix (nmod_poly_mat_mul_vandermonde1, and vandermonde2 which
        uses the points in pairs +-x); Waksman's algorithm; and a geometric
        progression (nmod_poly_mat_mul_geometric).

        The rule, for dimension above 4 (Waksman below, as before):
        vandermonde1 up to result length vdm1_len; vandermonde2, for the
        tiers of lengths up to 63, 127, ..., 2047, from the dimensions
        vdm2_dim; otherwise sd_fft_matmul from dimension mm_dim, except for
        long products (length above mm_len) below dimension mm_long_dim,
        and sd_fft_direct; see the tables below.

        The modulus enters through two properties: whether fft_small uses
        a single transform (single_prime), and the size of p: nmod_mat_mul
        modulo p itself, which the Vandermonde routines (and sd_fft_matmul
        for a single transform) spend their time in, is fastest up to 32
        or 40 bits and within the fast kernels up to 52 bits.
    */

    const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);

    /* the cheap part of what makes fft_small use a single transform
     * rather than several CRT primes, i.e. np == 1 (primality is
     * left to the plan, which checks it) */
    const int single_prime = (modbits <= 50)
        && ((slong) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) len + 3));

    /* 0, 1: a single transform, p of at most 32 bits or more; several
       transforms and 2: p of at most PML_NMOD_MAT_MUL_SMALL_BITS bits,
       3: at most 52 bits, 4: more */
    const int kind = single_prime ? (modbits <= 32 ? 0 : 1)
                                  : (modbits <= PML_NMOD_MAT_MUL_SMALL_BITS ? 2
                                     : (modbits <= 52 ? 3 : 4));
    /* vandermonde1 up to result length vdm1_len; vandermonde2 for result
       lengths up to 63, 127, ..., 2047 from the dimensions vdm2_dim (N:
       never); sd_fft_matmul from dimension mm_dim, for result lengths up to
       mm_len unless the dimension is at least mm_long_dim */
#define N WORD_MAX
    typedef struct { slong vdm1_len, vdm2_dim[6], mm_dim, mm_len, mm_long_dim; } _thr_t;
#if PML_NMOD_MAT_MUL_IFMA
    static const _thr_t thr[5] = {
        { 23, { N, N,  N,   N,  N, N }, 16,    N,  0 },
        { 23, { N, N,  N,   N,  N, N }, 16,    N,  0 },
        { 31, { 8, 8, 16, 128,  N, N }, 16,    N,  0 },
        { 31, { 8, 8,  8,  16, 32, N }, 16, 3071, 64 },
        { 11, { N, N,  N,   N,  N, N }, 16, 3071, 64 } };
#else
    static const _thr_t thr[5] = {
        { 23, { 512,  N,   N,   N, N, N }, 12,    N,   0 },
        {  5, { 512,  N,   N,   N, N, N }, 48,    N,   0 },
        { 23, {   8,  8,  16,  32, 64, N }, 48, 1023, 128 },
        { 23, {  24, 32,  48, 128, N, N }, 96,    N,   0 },
        { 11, {   8, 32, 128, 256, N, N }, 48,  767, 128 } };
#endif
#undef N
    const _thr_t * t = thr + kind;

    /* Waksman divides by 2, and the Vandermonde routines invert differences
       of their points, which the cardinality tests of the CAN_USE macros
       do not ensure for a composite modulus (with modn = 1000, Waksman is
       silently wrong and vandermonde2 raises "Impossible inverse"). The
       fft_small routines work modulo their own primes and are correct for
       any modulus, so they take those cases. The primality test is only
       reached for the shapes where a Vandermonde routine is wanted. */
    if (dim <= 4)
    {
        if (len < 31 && (modn & 1))
        {
            nmod_poly_mat_mul_waksman(res, pmat1, pmat2);
            return;
        }
    }
    else if (len <= t->vdm1_len)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE1(modn, len) && n_is_prime(modn))
        {
            nmod_poly_mat_mul_vandermonde1(res, pmat1, pmat2);
            return;
        }
    }
    else
    {
        slong i, tier_len;
        for (i = 0, tier_len = 63; i < 6; i++, tier_len = 2 * tier_len + 1)
        {
            if (len <= tier_len && dim >= t->vdm2_dim[i])
            {
                if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && n_is_prime(modn))
                {
                    nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
                    return;
                }
                break;
            }
        }
    }

    if (dim >= t->mm_dim && (len <= t->mm_len || dim >= t->mm_long_dim))
        nmod_poly_mat_mul_sd_fft_matmul(res, pmat1, pmat2);
    else
        nmod_poly_mat_mul_sd_fft_direct(res, pmat1, pmat2);
    return;

#endif /* PML_HAVE_MACHINE_VECTORS */

    /* the evaluation points of geometric and vandermonde2 need inverses of
       their differences, which the cardinality tests of the CAN_USE macros
       do not ensure for a composite modulus; Waksman needs an odd one (see
       above). The primality test (< 1us for 64-bit moduli) comes
       last, so that it is only reached when one of these routines is
       wanted, never for the small products. */
    /* TODO this is (as in some other places) a situation where one would like
       to have a suitable point available in some nmod context provided to
       functions (this points would give a good geometric progression) */

    if (dim > 12)
    {
        if (NMOD_POLY_CAN_USE_GEOMETRIC(modn, len) && len > 300 && n_is_prime(modn))
            nmod_poly_mat_mul_geometric(res, pmat1, pmat2);
        else if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && n_is_prime(modn))
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
        /* FIXME small fields (cannot use geom/vdm): should probably call waksman in some cases... */
    }

    else if (dim > 10)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && len < 200 && n_is_prime(modn))
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
            nmod_poly_mat_mul_waksman(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
    }

    else if (dim > 8)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && len < 120 && n_is_prime(modn))
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
            nmod_poly_mat_mul_waksman(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
    }

    else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
        nmod_poly_mat_mul_waksman(res, pmat1, pmat2);

    else
        nmod_poly_mat_mul(res, pmat1, pmat2);
}

