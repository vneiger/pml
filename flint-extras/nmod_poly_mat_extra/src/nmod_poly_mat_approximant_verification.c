#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"

int is_zero_mod_xk(const nmod_poly_mat_t mat, slong k)
{
    nmod_poly_t P;
    nmod_poly_init(P, mat->modulus);
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
        {
            nmod_poly_set(P, nmod_poly_mat_entry(mat, i, j));
            nmod_poly_shift_right(P, P, k);
            if (!nmod_poly_is_zero(P))
                return 0;
        }
    return 1;
}


/** TO FIX **/
int is_minimal_approximant_basis(const nmod_poly_mat_t base,
                                 const nmod_mat_t mat, slong order,
                                 const slong *shifts)
{
    slong rdim = mat->r, cdim = mat->c;
    mp_limb_t prime = mat->mod.n;
    nmod_poly_t constant;
    nmod_poly_mat_t mat_poly, res_mul;

    if (base->c != base->r)
    {
        printf("not basis: wrong shape");
        return 0;
    }

    nmod_poly_mat_init(mat_poly, rdim, cdim, prime);
    slong alloc;
    nmod_poly_init(constant, prime);
    for (slong i = 0; i < rdim; i++)
        for (slong j = 0; j < cdim; j++)
        {
            alloc = (slong) nmod_mat_get_entry(mat, i, j);
            nmod_poly_set_coeff_ui(constant, 0, alloc);
            nmod_poly_set(nmod_poly_mat_entry(mat_poly, i, j), constant);
        }
    nmod_poly_mat_init(res_mul, rdim, cdim, prime);
    nmod_poly_mat_mul(res_mul, base, mat_poly);

    if (! is_zero_mod_xk(res_mul, order))
    {
        printf("not zero");
        return 0;
    }
    slong lead_pos[rdim];
    leading_positions(lead_pos, base, shifts, ROW_WISE);
    printf("\nleading positions\n");
    for (slong i = 0; i < rdim; i++)
        printf("%lu ", lead_pos[i]);
    return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
