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

/* TODO would benefit from automatic tuning */
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
        - vandermonde1 up to result length 15 (31 for p <= 40 bits, 11
          above 52 bits);
        - vandermonde2, for several transforms and p <= 52 bits, up to
          length 63 from dimension 8, 127 from 12, then 255 from 128 and
          511 from 384 for p <= 40 bits, 255 from 16 and 511 from 96 for
          p <= 52 bits;
        - otherwise sd_fft_matmul from dimension 16 (24 for p <= 52 bits
          with several transforms, 32 above 52 bits), except for long
          products at moderate dimensions (length above 2047 below
          dimension 96 for p <= 40 bits, above 1023 below dimension 64
          for larger p), where its grouping under the memory bound costs
          more than it saves; sd_fft_direct otherwise.

        The modulus enters through two properties: whether fft_small uses
        a single transform (single_prime), and, when it does not, the size
        of p: nmod_mat_mul modulo p itself, which the Vandermonde routines
        spend their time in, is fastest up to 40 bits and within the fast
        kernels up to 52 bits.
    */

    const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);

    /* the cheap part of what makes fft_small use a single transform
     * rather than several CRT primes, i.e. np == 1 (primality is
     * left to the plan, which checks it) */
    const int single_prime = (modbits <= 50)
        && ((slong) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) len + 3));

    /* 0: one transform; several transforms and 1: p of at most 40 bits,
       2: p of at most 52 bits, 3: larger p */
    const int kind = single_prime ? 0 : (modbits <= 40 ? 1 : (modbits <= 52 ? 2 : 3));
    /* vandermonde1 up to this result length */
    static const slong vdm1_len[4] = { 15, 31, 15, 11 };
    /* vandermonde2 for result lengths up to 63, 127, 255, 511 from these
       dimensions (WORD_MAX: never) */
    static const slong vdm2_dim[4][4] = {
        { WORD_MAX, WORD_MAX, WORD_MAX, WORD_MAX },
        { 8, 12, 128, 384 },
        { 8, 12, 16, 96 },
        { WORD_MAX, WORD_MAX, WORD_MAX, WORD_MAX } };
    /* sd_fft_matmul from dimension matmul_dim, for result lengths up to
       matmul_len unless the dimension is at least matmul_long_dim */
    static const slong matmul_dim[4] = { 16, 16, 24, 32 };
    static const slong matmul_len[4] = { WORD_MAX, 2047, 1023, 1023 };
    static const slong matmul_long_dim[4] = { 0, 96, 64, 64 };

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
    else if (len <= vdm1_len[kind])
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
        for (i = 0, tier_len = 63; i < 4; i++, tier_len = 2 * tier_len + 1)
        {
            if (len <= tier_len && dim >= vdm2_dim[kind][i])
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

    if (dim >= matmul_dim[kind] && (len <= matmul_len[kind] || dim >= matmul_long_dim[kind]))
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

