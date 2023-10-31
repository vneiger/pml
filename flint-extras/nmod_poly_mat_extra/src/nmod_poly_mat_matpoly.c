#include "nmod_poly_mat_mat_poly.h"
#include "nmod_poly_mat_utils.h"

void nmod_mat_poly0_init(nmod_mat_poly0_t matp,
                        slong degree,
                        slong length,
                        slong r,
                        slong c,
                        mp_limb_t mod)
{

	slong i;

	if (r > 0)
	{
		matp->mat = (nmod_mat_struct *) flint_malloc(length * sizeof(nmod_mat_struct));
		for (i = 0; i < length; i++) // TODO STOP AT DEGREE?
			nmod_mat_init(matp->mat + i, r, c, mod);
	}
	else
		matp->mat = NULL;


	matp->length = length;
	matp->degree = degree;
	matp->mod = mod;
	matp->r = r;
	matp->c = c;
}

void nmod_mat_poly0_clear(nmod_mat_poly0_t A)
{
	if (A == NULL)
		return;
	if (A->mat != NULL)
	{
		for (slong i = 0; i < A->length; i++)
			nmod_mat_clear(A->mat + i);

		flint_free(A->mat);
	}
}

void nmod_mat_poly0_print(const nmod_mat_poly0_t A)
{
	for (slong i = 0; i < A->length; i++)
		nmod_mat_print_pretty(A->mat + i);
}

/** suppose that res->length >= degree + 1 **/
void nmod_mat_poly0_set(nmod_mat_poly0_t res, const nmod_poly_mat_t F)
{
	slong degree, r, c;
	mp_limb_t mod;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	mod = nmod_poly_mat_modulus(F);

	nmod_mat_t mat;
	nmod_mat_init(mat, r, c, mod);

	for (slong i = 0; i <= degree; i++)
	{
		nmod_poly_mat_get_coeff_mat(mat, F, i);
		nmod_mat_set(res->mat + i, mat);
	}
	nmod_mat_clear(mat);

	res->degree = degree;
}

void nmod_mat_poly0_init_set(nmod_mat_poly0_t res,
                            const nmod_poly_mat_t F)
{
	slong degree, r, c;
	mp_limb_t mod;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	mod = nmod_poly_mat_modulus(F);

	nmod_mat_poly0_init(res, degree, degree + 1, r, c, mod);

	nmod_poly_struct *P;
	for (slong i = 0; i < r; i++)
		for (slong j = 0; j < c; j++)
		{
			P = nmod_poly_mat_entry(F, i, j);
			for (slong k = 0; k <= nmod_poly_degree(P); k++)
				nmod_mat_set_entry(res->mat + k, i, j,nmod_poly_get_coeff_ui(P, k));
		}
}

void nmod_mat_poly0_init_setII(nmod_mat_poly0_t res,
                              const nmod_poly_mat_t F,
                              slong length)
{
	slong degree, r, c, min, d;
	mp_limb_t mod;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	mod = nmod_poly_mat_modulus(F);

	nmod_mat_poly0_init(res, degree, length, r, c, mod);

	nmod_poly_struct *P;
	for (slong i = 0; i < r; i++)
		for (slong j = 0; j < c; j++)
		{
			P = nmod_poly_mat_entry(F, i, j);
			d = nmod_poly_degree(P) + 1;
			min = (length < d) ? length : d;
			for (slong k = 0; k < min; k++)
				nmod_mat_set_entry(res->mat + k, i, j,nmod_poly_get_coeff_ui(P, k));
		}
}

void nmod_mat_poly0_init_setIII(nmod_mat_poly0_t res,
                               const nmod_poly_mat_t F,
                               slong length)
{
	slong degree, r, c;
	mp_limb_t mod;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	mod = nmod_poly_mat_modulus(F);

	nmod_mat_poly0_init(res, degree, degree + length, r, c, mod);

	nmod_poly_struct *P;
	for (slong i = 0; i < r; i++)
		for (slong j = 0; j < c; j++)
		{
			P = nmod_poly_mat_entry(F, i, j);
			for (slong k = 0; k <= nmod_poly_degree(P); k++)
				nmod_mat_set_entry(res->mat + k, i, j, nmod_poly_get_coeff_ui(P, k));
		}
}

void nmod_mat_poly0_get_coef(nmod_mat_t res,
                            const nmod_mat_poly0_t F,
                            slong k)
{
	nmod_mat_set(res, F->mat + k);
}

void nmod_mat_poly0_naive_mul_coef(nmod_mat_t res,
                                  const nmod_mat_poly0_t A,
                                  const nmod_mat_poly0_t B,
                                  slong k)
{
	nmod_mat_t temp;
	slong A_r, A_c, A_degree, B_r, B_c, B_degree;
	mp_limb_t A_mod, B_mod;

	A_r = A->r;
	A_c = A->c;
	A_mod = A->mod;
	A_degree = A->degree;

	B_r = B->r;
	B_c = B->c;
	B_mod = B->mod;
	B_degree = B->degree;

	if (A_mod != B_mod)
	{
		printf("\nERROR! Wrong modulus: nmod_mat_poly0_naive_mul_coef\n");
		return;
	}
	if (A_c != B_r)
	{
		printf("\nERROR! Wrong shape: nmod_mat_poly0_naive_mul_coef\n");
		printf("shape A = (%ld, %ld), shape B = (%ld, %ld)\n", A_r, A_c, B_r, B_c);
		return;
	}

	nmod_mat_zero(res);
	nmod_mat_init(temp, A_r, B_c, A_mod);
	for (slong i = 0; i <= k; i++)
	{
		if (i <= A_degree && (k-i) <= B_degree)
		{
			nmod_mat_mul(temp, A->mat + i, B->mat + (k - i));
			nmod_mat_add(res, res, temp);
		}
	}
	nmod_mat_clear(temp);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
