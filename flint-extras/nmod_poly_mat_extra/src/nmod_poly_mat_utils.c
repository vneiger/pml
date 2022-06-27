#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"

void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, slong degree)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            nmod_mat_set_entry(res, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j), degree));
}

void nmod_poly_mat_shift(nmod_poly_mat_t res, slong k)
{
    nmod_poly_struct *P;
    if (k > 0)
    {
        for (slong i = 0; i < res->r; i++)
            for (slong j = 0; j < res->c; j++)
            {
                P = nmod_poly_mat_entry(res, i, j);
                if (!nmod_poly_is_zero(P))
                    nmod_poly_shift_left(P, P, k);
            }
        return;
    }


    if (k < 0)
    {
        for (slong i = 0; i < res->r; i++)
            for (slong j = 0; j < res->c; j++)
            {
                P = nmod_poly_mat_entry(res, i, j);
                if (!nmod_poly_is_zero(P))
                    nmod_poly_shift_right(P, P, k);
            }
        return;
    }
}


void nmod_mat_to_poly_mat(nmod_poly_mat_t res, const nmod_mat_t M)
{
	slong r = M->r, c = M->c;
	nmod_poly_t P;
	nmod_poly_init(P, res->modulus);
	for (slong i = 0; i < r; i++)
		for (slong j = 0; j < c; j++)
		{
			nmod_poly_set_coeff_ui(P, 0, nmod_mat_get_entry(M, i, j));
			nmod_poly_set(nmod_poly_mat_entry(res, i, j), P);
		}
	nmod_poly_clear(P);
}

void nmod_mat_poly_to_poly_mat(nmod_poly_mat_t res, const nmod_mat_poly_t F)
{
	slong degree, r, c;
	mp_limb_t mod;
	degree = F->degree;
	r = F->r;
	c = F->c;
	mod = F->mod;

	if (res->modulus != mod)
	{
		printf("\nERROR! Wrong modulus: nmod_mat_poly_to_poly_mat\n");
		return;
	}
	if (res->r != r || res->c != c)
	{
		printf("\nERROR! Wrong shape: nmod_mat_poly_to_poly_mat\n");
		printf("shape res = (%ld, %ld), shape F = (%ld, %ld)\n", res->r, res->c, r, c);
		return;
	}

	nmod_poly_mat_zero(res);


	nmod_poly_mat_t mat;
	nmod_poly_mat_init(mat, r, c, mod);

	for (slong i = 0; i <= degree; i++)
	{

		nmod_mat_to_poly_mat(mat, F->mat + i);
		nmod_poly_mat_shift(mat, i);
		nmod_poly_mat_add(res, res, mat);
	}
	nmod_poly_mat_clear(mat);
}



/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
