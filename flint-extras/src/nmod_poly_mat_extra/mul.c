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

        2026-10-01 Thresholds fitted on square products on Zen 4, one
        thread, against FLINT-dev with the u32/u52/fp50 kernels of
        nmod_mat_mul (PR #2842), over six moduli: 21-, 30-, 40- and 50-bit
        FFT primes (a single transform, np == 1), a 30-bit prime (np == 2)
        and a 60-bit prime (np == 3). Against the previous thresholds, on
        the same measurements: average slowdown with respect to the best
        routine 1.38 -> 1.03, 95th percentile 4.2 -> 1.2.

        What changed with those kernels is that nmod_mat_mul is now fast at
        every modulus up to 52 bits, which is what the Vandermonde routines
        and sd_fft_matmul spend their time in:

        - short products go to vandermonde1 (result length <= 15) or
          vandermonde2 (up to 127 or 255), at every dimension from 8 on,
          by 1.2-3x over the next routine; geometric no longer wins
          anywhere by more than a few percent and is not selected
          (this would probably change for very large matrices of large degree);
        - sd_fft_matmul takes over from sd_fft_direct at large dimension
          for every modulus, not only in a 20-24-bit window as before; its
          pointwise products are modulo the 50-bit fft_small primes, or
          modulo p itself for an FFT prime, all within the new kernels.

        The modulus enters through two properties: whether fft_small uses
        a single transform (single_prime), and, when it does not, whether
        nmod_mat_mul modulo p itself is within the fast kernels (at most
        52 bits) -- which is what the Vandermonde routines need.
    */

    const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);

    /* the cheap part of what makes fft_small use a single transform
     * rather than several CRT primes, i.e. np == 1 (primality is
     * left to the plan, which checks it) */
    const int single_prime = (modbits <= 50)
        && ((slong) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) len + 3));

    /* 0: one transform; 1: several, p of at most 52 bits; 2: several, larger p */
    const int kind = single_prime ? 0 : (modbits <= 52 ? 1 : 2);
    static const slong vdm2_len[3] = { 127, 255, 127 };   /* vandermonde2 up to this length */
    static const slong vdm2_dim[3] = {  32,  16,  48 };   /* ... from this dimension on */
    static const slong matmul_dim[3] = { 192, 96, 192 };  /* sd_fft_matmul from this dimension */

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
    else if (len <= 15)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE1(modn, len) && n_is_prime(modn))
        {
            nmod_poly_mat_mul_vandermonde1(res, pmat1, pmat2);
            return;
        }
    }
    else if (len <= vdm2_len[kind] && dim >= vdm2_dim[kind])
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && n_is_prime(modn))
        {
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
            return;
        }
    }

    if (len >= 127 && dim >= matmul_dim[kind])
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

