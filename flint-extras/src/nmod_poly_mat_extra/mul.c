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
        Evaluation-interpolation, either at the roots of unity of
        fft_small (nmod_poly_mat_mul_sd_fft_direct and
        nmod_poly_mat_mul_sd_fft_matmul) or at a geometric progression in
        Z/p itself (nmod_poly_mat_mul_geometric).

        NOTE
        The FFT routes evaluate at np * ztrunc points, where np is from 1
        (fft_small can transform directly modulo p) to 3 or even 4 in rare
        cases. One thing is that ztrunc, the transform length, is rounded up
        which increases the number of points (and of pointwise products),
        potentially by up to 2. Also, the number of pointwise multiplications
        is multiplied by `np` since they have to be done for each prime. The
        geometric route evaluates at exactly len points and does exactly
        len pointwise multiplications. The two FFT routes differ in the
        pointwise stage only: sd_fft_direct multiplies the transforms
        entry by entry with its own kernels, sd_fft_matmul hands each
        evaluation point to nmod_mat_mul.

        2026-09-13 Thresholds fitted on square products over three machines
        (Zen 4, Ice Lake, Apple M4), four moduli (a 21-bit and a 50-bit FFT
        prime, a 30-bit and a 60-bit prime), at 1 thread, against the
        current FLINT-dev built without an external BLAS.
    */

    const flint_bitcnt_t modbits = FLINT_BIT_COUNT(modn);

    /* the cheap part of what makes fft_small use a single transform
     * rather than several CRT primes, i.e. np == 1 (primality is
     * left to the plan, which checks it) */
    /* TODO use some function already in fft_small for checking if prime is FFT of sufficient depth? */
    const int single_prime = (modbits <= 50)
        && ((slong) flint_ctz(modn - 1) >= FLINT_BIT_COUNT((ulong) len + 3));

    /*
        Whether that single transform is modulo p itself -- fft_small
        does that only from 20 bits up, below which it uses one of its own
        50-bit primes (see _nmod_poly_should_directly_fft in
        fft_small/plan.c) -- and the resulting pointwise matrix products
        are then at a modulus small enough for nmod_mat_mul to take its
        fastest route: one gemm over doubles with no chinese remaindering,
        which it does when the smallest dimension is above 100 and
        FLINT_BIT_COUNT(k) + 2*bits < 58 (see nmod_mat/mul.c).

        This is what pays for the extra evaluations of the matmul variant,
        and sd_fft_direct cannot follow: its pointwise kernels work on the
        transforms themselves and do not get cheaper as p shrinks. The
        window is narrow -- 20 bits up to about 24 -- but inside it the
        matmul variant is the fastest route by a wide margin: at dimension
        512 and length 63 it is 1.8 times faster than sd_fft_direct on
        Zen 4 and 4 times on Apple M4.
    */
    const int fast_matmul = single_prime && modbits >= 20 && dim > 100
        && (FLINT_BIT_COUNT((ulong) pmat1->c) + 2 * modbits < 58);

    int use_fft = 0, use_matmul = 0, use_geometric = 0;

    if (fast_matmul)
    {
        if (len >= 32 || dim <= 128)
            use_matmul = 1;
        /* len < 32 and dim > 128: the tiers below */
    }
    else if (single_prime)
    {
        if (len >= 128 || dim <= (len >= 32 ? 256 : 128))
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
            /* the FFT route pays for the points it adds by rounding the
               transform up, and the more matrix entries there are to
               transform the less that costs relative to the pointwise
               stage: the crossover sits at about dim == len */
            if (dim < len)
                use_fft = 1;
            else if (dim <= 384)
                use_geometric = 1;
            /* dim > 384: the tiers below */
        }
        else if (dim >= 24 && dim <= 192)
            use_geometric = 1;
    }

    if (use_fft)
    {
        nmod_poly_mat_mul_sd_fft_direct(res, pmat1, pmat2);
        return;
    }
    if (use_matmul)
    {
        nmod_poly_mat_mul_sd_fft_matmul(res, pmat1, pmat2);
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

