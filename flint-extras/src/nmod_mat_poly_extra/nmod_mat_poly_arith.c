#include <flint/nmod_mat.h>
#include "nmod_mat_poly.h"

void nmod_mat_poly_mul_coeff(nmod_mat_t coeff,
                             const nmod_mat_poly_t mat1,
                             const nmod_mat_poly_t mat2,
                             slong k)
{
    // only consider indices i such that:
    //     0 <= i <= k
    //     i < mat1->length
    //     k-i < mat2->length
    // so  i < min(k+1,mat1->length)
    // and i >= max(0, k + 1 - mat2->length)
    const slong ubound = FLINT_MIN(k+1, mat1->length);
    const slong lbound = FLINT_MAX(0, k+1 - mat2->length);

    // lbound >= ubound ==> coeff is zero, no term to consider
    if (lbound >= ubound)
    {
        nmod_mat_zero(coeff);
        return;
    }

    // now lbound < ubound
    // first handle i == lbound separately, to avoid wasting time zero-ing `coeff`
    nmod_mat_t cmat1, cmat2;
    nmod_mat_poly_coeff_attach(cmat1, mat1, lbound);
    nmod_mat_poly_coeff_attach(cmat2, mat2, k - lbound);
    nmod_mat_mul(coeff, cmat1, cmat2);

    // `if` just here to avoid initializing temp for nothing
    if (lbound + 1 < ubound)
    {
        nmod_mat_t temp;
        nmod_mat_init(temp, mat1->r, mat2->c, mat1->mod.n);
        for (slong i = lbound+1; i < ubound; i++)
        {
            nmod_mat_poly_coeff_attach(cmat1, mat1, i);
            nmod_mat_poly_coeff_attach(cmat2, mat2, k - i);
            nmod_mat_mul(temp, cmat1, cmat2);
            nmod_mat_add(coeff, coeff, temp);
        }
        nmod_mat_clear(temp);
    }
}

void nmod_mat_poly_evaluate_nmod(nmod_mat_t eval,
                                 const nmod_mat_poly_t matp,
                                 ulong pt)
{
    slong k = matp->length;
    nmod_mat_t cmat;

    if (k == 0)
    {
        nmod_mat_zero(eval);
        return;
    }

    if (k == 1 || pt == 0)
    {
        nmod_mat_poly_coeff_attach(cmat, matp, 0);
        nmod_mat_set(eval, cmat);
        return;
    }

    k--; // k == degree
    nmod_mat_poly_coeff_attach(cmat, matp, k);
    nmod_mat_set(eval, cmat);
    k--; // k == degree-1

    // Horner: eval = matp[k] + eval*pt, k = degree-1 ... 0
    for ( ; k >= 0; k--)
    {
        nmod_mat_poly_coeff_attach(cmat, matp, k);
        nmod_mat_scalar_addmul_ui(eval, cmat, eval, pt);
    }
}
