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

    const slong dim = n_cbrt(pmat1->r * pmat1->c * pmat2->c);
    const slong len = len1 + len2 - 1;
    const ulong modn = pmat1->modulus;

#if PML_HAVE_MACHINE_VECTORS
    /*
        Evaluation-interpolation, either at the roots of unity of
        fft_small (nmod_poly_mat_mul_sd_fft_direct) or at a geometric
        progression in Z/p itself (nmod_poly_mat_mul_geometric).

        NOTE
        The FFT route evaluates at np * ztrunc points, where np is from 1
        (fft_small can transform directly modulo p) to 3 or even 4 in rare
        cases. One thing is that ztrunc, the transform length, is rounded up
        which potentially multiplies the number of points (and of pointwise
        products) by up to 2. Also, the number of pointwise multiplications
        is multiplied by `np` since they have to be done for each prime. The
        geometric route evaluates at exactly len points.

        2026-09-13 Thresholds fitted on square products over three machines
        (Zen 4, Ice Lake, Apple M4), three moduli (a 50-bit FFT prime, a 30-bit
        and a 60-bit prime), at 1, 2 and 8 threads, against the current
        FLINT-dev built without BLAS.

        TODO fitted on square products; fat and thin ones are not covered.

        TODO measurements were focused on "small" matrices (but products
        already taking up to about 60s), with dimension <= 512. The fft_matmul
        strategy is not called currently but should be interesting in a corner
        case of large instances: with matrices of large dimensions (gaining
        from nmod_mat_mul) of large degree (otherwise the Vandermonde approach
        may be faster) and when there is a single FFT prime (otherwise
        geometric may be faster).

        TODO see if BLAS brings big changes to these thresholds

        NOTE: unlike the other variants, these two allocate the
        evaluations of the operands, up to a soft budget which is the
        larger of a fixed floor (currently 256MB) and twice the size of
        the operands and the result (see the MEM_FLOOR constants).
    */

    /* the cheap part of what makes fft_small transform directly modulo p
     * rather than modulo several CRT primes, i.e. np == 1 (primality is
     * left to the plan, which checks it) */
    /* TODO use some function already in fft_small for checking if prime is FFT of sufficient depth? */
    const int single_prime = (FLINT_BIT_COUNT(modn) <= 50)
        && ((slong) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) len + 3));
    int use_fft = 0, use_geometric = 0;

    if (single_prime)
    {
        if (len >= 128 || dim <= (len >= 32 ? 384 : 128))
            use_fft = 1;
        else if (len >= 32)
            use_geometric = 1;
        /* len < 32 and dim > 128: the tiers below */
    }
    else if (len >= 128)
    {
        if (dim >= 384)
            use_geometric = 1;
        else
            use_fft = 1;
    }
    else if (dim >= 8)   /* several primes, len < 128 */
    {
        if (len >= 63)
        {
            if (dim < 64)
                use_fft = 1;
            else if (dim <= 384)
                use_geometric = 1;
        }
        else if (dim >= 24 && dim <= 192)
            use_geometric = 1;
    }

    if (use_fft)
    {
        nmod_poly_mat_mul_sd_fft_direct(res, pmat1, pmat2);
        return;
    }
    if (use_geometric && NMOD_POLY_CAN_USE_GEOMETRIC(modn, len))
    {
        nmod_poly_mat_mul_geometric(res, pmat1, pmat2);
        return;
    }
#endif /* PML_HAVE_MACHINE_VECTORS */

    if (dim > 12)
    {
        if (NMOD_POLY_CAN_USE_GEOMETRIC(modn, len) && len > 300)
            nmod_poly_mat_mul_geometric(res, pmat1, pmat2);
        else if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len))
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
        /* FIXME small fields (cannot use geom/vdm): should probably call waksman in some cases... */
    }

    else if (dim > 10)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && len < 200)
            nmod_poly_mat_mul_vandermonde2(res, pmat1, pmat2);
        else if (NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn))
            nmod_poly_mat_mul_waksman(res, pmat1, pmat2);
        else
            nmod_poly_mat_mul(res, pmat1, pmat2);
    }

    else if (dim > 8)
    {
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) && len < 120)
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

