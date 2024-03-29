#include "nmod_mat_poly.h"

void _nmod_mat_poly_shift_left(nmod_mat_struct * smatp,
                               const nmod_mat_struct * matp,
                               slong len,
                               slong n)
{
    if (smatp != matp)
    {
        for (slong i = 0; i < len; i++)
            nmod_mat_set(smatp + n + i, matp + i);
    }
    else
    {
        /* Copy in reverse to avoid writing over unshifted coefficients */
        for (slong i = len-1; i >= 0; i--)
            nmod_mat_swap(smatp + n + i, smatp + i);
    }

    for (slong i = 0; i < n; i++)
        nmod_mat_zero(smatp + i);
}

//void _nmod_mat_poly_shift_right(nmod_mat_poly_t smatp,
//                                const nmod_mat_poly_t matp,
//                                slong len,
//                                slong n)
//{
//}

void nmod_mat_poly_shift_left(nmod_mat_poly_t smatp,
                              const nmod_mat_poly_t matp,
                              slong n)
{
    if (n == 0)
    {
        nmod_mat_poly_set(smatp, matp);
        return;
    }

    if (matp->length == 0)
    {
        nmod_mat_poly_zero(smatp);
        return;
    }

    nmod_mat_poly_fit_length(smatp, matp->length + n);
    _nmod_mat_poly_set_length(smatp, matp->length + n);
    _nmod_mat_poly_shift_left(smatp->coeffs, matp->coeffs, matp->length - n, n);
}

//void nmod_mat_poly_shift_right(nmod_mat_poly_t smatp,
//                               const nmod_mat_poly_t matp,
//                               slong n);


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
