#include "nmod_poly_mat_utils.h"  // truncate, shift
#include "nmod_poly_mat_approximant.h"

// computes x**(-d) (A*B mod x**h)
void middle_product(nmod_poly_mat_t res,
                    const nmod_poly_mat_t A,
                    const nmod_poly_mat_t B,
                    slong d,
                    slong h)
{
    nmod_poly_mat_mul(res, A, B);
    nmod_poly_mat_truncate(res, h);
    nmod_poly_mat_shift_right(res, res, d);
}


void pmbasis(nmod_poly_mat_t res, slong *res_shifts,
              const nmod_poly_mat_t F, ulong sigma, const slong *shifts)
{
    slong rdim = F->r, cdim = F->c;
    mp_limb_t prime = F->modulus;
    nmod_poly_mat_t Pl, Ph, F_prime;
    slong ul[rdim];

    if (sigma <= PMBASIS_THRES)
    {
        mbasisIII(res, res_shifts, F, sigma, shifts);
        return;
    }

    nmod_poly_mat_init(Pl, rdim, rdim, prime);
    pmbasis(Pl, ul, F, sigma / 2, shifts);

    nmod_poly_mat_init(F_prime, rdim, cdim, prime);
    middle_product(F_prime, Pl, F, sigma / 2, sigma);

    nmod_poly_mat_init(Ph, rdim, rdim, prime);
    if (sigma % 2 == 0)
        pmbasis(Ph, res_shifts, F_prime, sigma / 2, ul);
    else
        pmbasis(Ph, res_shifts, F_prime, sigma / 2 + 1, ul);

    nmod_poly_mat_mul(res, Ph, Pl);

    nmod_poly_mat_clear(Ph);
    nmod_poly_mat_clear(Pl);
    nmod_poly_mat_clear(F_prime);
}

void nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           slong order)
{
    if (order <= PMBASIS_THRES)
    {
        nmod_poly_mat_mbasis(appbas, shift, pmat, order);
        return;
    }

    const long order1 = order>>1;
    const long order2 = order - order1;
    nmod_poly_mat_t appbas2, residual;

    nmod_poly_mat_init(appbas2, pmat->r, pmat->r, pmat->modulus);
    nmod_poly_mat_init(residual, pmat->r, pmat->c, pmat->modulus);

    nmod_poly_mat_pmbasis(appbas, shift, pmat, order1);

    middle_product(residual, appbas, pmat, order1, order);

    nmod_poly_mat_pmbasis(appbas2, shift, residual, order2);

    nmod_poly_mat_mul(appbas, appbas2, appbas);

    nmod_poly_mat_clear(appbas2);
    nmod_poly_mat_clear(residual);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
