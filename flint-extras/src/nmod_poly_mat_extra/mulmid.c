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

        The rule, by class of modulus:
        - the Vandermonde routines up to nhi = 16 at dimension at most 4,
          and above up to
            single transform:          min(max(24, dim), 96)
            several, p <= 40 bits:     min(8 dim, 384)
            several, p <= 52 bits:     min(16 dim, 1024)
            several, p > 52 bits:      min(max(16, dim/2), 48)
          vandermonde1 up to nhi = 24 (32 for p <= 40 bits, 9 above 52
          bits), vandermonde2 above; they need a prime modulus, and naive
          takes their place otherwise up to nhi = 16;
        - the windowed fft_small routines everything else: sd_fft_matmul
          from dimension 128 (96 with several transforms and p of at most
          52 bits), except for nhi above 2048 below dimension 192, and
          sd_fft_direct otherwise; they transform at length about nhi
          (a power of two N fitting the window, or the transposed product
          when it is shorter, see impl.h). These thresholds were fitted
          before the transposed product, which only makes the fft_small
          routines faster.
    */
    {
        const ulong modn = pmat1->modulus;
        const slong dim = n_cbrt(pmat1->r * pmat1->c * pmat2->c);
        const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);
        const int single_prime = (modbits <= 50)
            && ((flint_bitcnt_t) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) (len1 + len2 + 2)));
        slong vdm_nhi, vdm1_nhi, matmul_dim;

        if (single_prime)
        {
            vdm_nhi = FLINT_MIN(FLINT_MAX(24, dim), 96);
            vdm1_nhi = 24;
            matmul_dim = 128;
        }
        else if (modbits <= 40)
        {
            vdm_nhi = FLINT_MIN(8 * dim, 384);
            vdm1_nhi = 32;
            matmul_dim = 96;
        }
        else if (modbits <= 52)
        {
            vdm_nhi = FLINT_MIN(16 * dim, 1024);
            vdm1_nhi = 24;
            matmul_dim = 96;
        }
        else
        {
            vdm_nhi = FLINT_MIN(FLINT_MAX(16, dim / 2), 48);
            vdm1_nhi = 9;
            matmul_dim = 128;
        }
        if (dim <= 4)
            vdm_nhi = 16;

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

        if (dim >= matmul_dim && (nhi <= 2048 || dim >= 192))
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
