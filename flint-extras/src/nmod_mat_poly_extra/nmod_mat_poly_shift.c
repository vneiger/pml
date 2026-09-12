#include "nmod_mat_poly.h"

void _nmod_mat_poly_shift_left(nmod_mat_poly_t smatp,
                               const nmod_mat_poly_t matp,
                               slong len,
                               slong n)
{
    nmod_mat_t dst, src;

    if (smatp != matp)
    {
        for (slong i = 0; i < len; i++)
        {
            nmod_mat_poly_coeff_attach(dst, smatp, n + i);
            nmod_mat_poly_coeff_attach(src, matp, i);
            nmod_mat_set(dst, src);
        }
    }
    else
    {
        /* Exchange the entry arrays; in reverse, to avoid writing over
           unshifted coefficients */
        for (slong i = len-1; i >= 0; i--)
            FLINT_SWAP(nn_ptr, smatp->coeffs[n + i], smatp->coeffs[i]);
    }

    for (slong i = 0; i < n; i++)
    {
        nmod_mat_poly_coeff_attach(dst, smatp, i);
        nmod_mat_zero(dst);
    }
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

    /* read the length of the input before growing the output: when the two
       are aliased, _nmod_mat_poly_set_length below changes it */
    const slong len = matp->length;

    nmod_mat_poly_fit_length(smatp, len + n);
    _nmod_mat_poly_set_length(smatp, len + n);
    _nmod_mat_poly_shift_left(smatp, matp, len, n);
}

//void nmod_mat_poly_shift_right(nmod_mat_poly_t smatp,
//                               const nmod_mat_poly_t matp,
//                               slong n);
