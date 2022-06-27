#include "nmod_poly_mat_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_approximant.h" // TODO for apply_perm ->  should move

void nmod_mat_poly_init(nmod_mat_poly_t res, slong degree, slong length,
		slong r, slong c, mp_limb_t modulus)
{

	slong i;

	if (r > 0)
	{
		res->mat = (nmod_mat_struct *) flint_malloc(length * sizeof(nmod_mat_struct));
		for (i = 0; i < length; i++)
			nmod_mat_init(res->mat + i, r, c, modulus);
	}
	else
		res->mat = NULL;


	res->length = length;
	res->degree = degree;
	res->modulus = modulus;
	res->r = r;
	res->c = c;
}

void nmod_mat_poly_clear(nmod_mat_poly_t A)
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

void nmod_mat_poly_print(const nmod_mat_poly_t A)
{
	for (slong i = 0; i < A->length; i++)
		nmod_mat_print_pretty(A->mat + i);
}

/** suppose that res->length >= degree + 1 **/
void nmod_mat_poly_set(nmod_mat_poly_t res, const nmod_poly_mat_t F)
{
	slong degree, r, c;
	mp_limb_t modulus;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	modulus = nmod_poly_mat_modulus(F);

	nmod_mat_t mat;
	nmod_mat_init(mat, r, c, modulus);

	for (slong i = 0; i <= degree; i++)
	{
		coefficient_matrix(mat, F, i);
		nmod_mat_set(res->mat + i, mat);
	}
	nmod_mat_clear(mat);

	res->degree = degree;
}

void nmod_mat_poly_init_set(nmod_mat_poly_t res, const nmod_poly_mat_t F)
{
	slong degree, r, c;
	mp_limb_t modulus;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	modulus = nmod_poly_mat_modulus(F);

	nmod_mat_poly_init(res, degree, degree + 1, r, c, modulus);

	nmod_poly_struct *P;
	for (slong i = 0; i < r; i++)
		for (slong j = 0; j < c; j++)
		{
			P = nmod_poly_mat_entry(F, i, j);
			for (slong k = 0; k <= nmod_poly_degree(P); k++)
				nmod_mat_set_entry(res->mat + k, i, j,nmod_poly_get_coeff_ui(P, k));
		}
}

void nmod_mat_poly_init_setII(nmod_mat_poly_t res,
		const nmod_poly_mat_t F, slong length)
{
	slong degree, r, c, min, d;
	mp_limb_t modulus;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	modulus = nmod_poly_mat_modulus(F);

	nmod_mat_poly_init(res, degree, length, r, c, modulus);

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

void nmod_mat_poly_init_setIII(nmod_mat_poly_t res,
		const nmod_poly_mat_t F, slong length)
{
	slong degree, r, c;
	mp_limb_t modulus;
	degree = nmod_poly_mat_degree(F);
	r = nmod_poly_mat_nrows(F);
	c = nmod_poly_mat_ncols(F);
	modulus = nmod_poly_mat_modulus(F);

	nmod_mat_poly_init(res, degree, degree + length, r, c, modulus);

	nmod_poly_struct *P;
	for (slong i = 0; i < r; i++)
		for (slong j = 0; j < c; j++)
		{
			P = nmod_poly_mat_entry(F, i, j);
			for (slong k = 0; k <= nmod_poly_degree(P); k++)
				nmod_mat_set_entry(res->mat + k, i, j, nmod_poly_get_coeff_ui(P, k));
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
	mp_limb_t modulus;
	degree = F->degree;
	r = F->r;
	c = F->c;
	modulus = F->modulus;

	if (res->modulus != modulus)
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
	nmod_poly_mat_init(mat, r, c, modulus);

	for (slong i = 0; i <= degree; i++)
	{

		nmod_mat_to_poly_mat(mat, F->mat + i);
		nmod_poly_mat_shift(mat, i);
		nmod_poly_mat_add(res, res, mat);
	}
	nmod_poly_mat_clear(mat);
}

void nmod_mat_poly_get_coef(nmod_mat_t res, const nmod_mat_poly_t F, slong k)
{
	nmod_mat_set(res, F->mat + k);
}

void nmod_mat_poly_naive_mul_coef(nmod_mat_t res, const nmod_mat_poly_t A,
		const nmod_mat_poly_t B, slong k)

{
	nmod_mat_t temp;
	slong A_r, A_c, A_degree, B_r, B_c, B_degree;
	mp_limb_t A_modulus, B_modulus;

	A_r = A->r;
	A_c = A->c;
	A_modulus = A->modulus;
	A_degree = A->degree;

	B_r = B->r;
	B_c = B->c;
	B_modulus = B->modulus;
	B_degree = B->degree;

	if (A_modulus != B_modulus)
	{
		printf("\nERROR! Wrong modulus: nmod_mat_poly_naive_mul_coef\n");
		return;
	}
	if (A_c != B_r)
	{
		printf("\nERROR! Wrong shape: nmod_mat_poly_naive_mul_coef\n");
		printf("shape A = (%ld, %ld), shape B = (%ld, %ld)\n", A_r, A_c, B_r, B_c);
		return;
	}

	nmod_mat_zero(res);
	nmod_mat_init(temp, A_r, B_c, A_modulus);
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

void structured_list_multiplication_blocks(nmod_mat_poly_t res, const nmod_mat_t A,
		const slong *perm, slong rank, slong deb, slong sigma)
{
	slong r = res->r, c = res->c;
	mp_limb_t modulus = res->modulus;
	slong i;
	nmod_mat_t previous_R1, R1, R2, R1_cp, R2_cp;
	nmod_mat_struct *r_i;
	slong *inv_perm = _perm_init(r);

	/** init **/
	nmod_mat_init(R1_cp, rank, c, modulus);
	nmod_mat_init(previous_R1, rank, c, modulus);
	nmod_mat_init(R2_cp, r - rank, c, modulus);

	for (i = deb; i < sigma; i++)
	{
		r_i = res->mat + i;
		apply_perm_rows_to_matrix(r_i, perm, r);
		nmod_mat_window_init(R1, r_i, 0, 0, rank, c);
		nmod_mat_set(R1_cp, R1);

		nmod_mat_set(R1, previous_R1); //it will set the top block of r_{i-1} on r_{i}

	nmod_mat_set(previous_R1, R1_cp);

	nmod_mat_window_init(R2, r_i, rank, 0, r, c);
	nmod_mat_set(R2_cp, R2);

	nmod_mat_mul(R2, A, R1_cp);
	nmod_mat_add(R2, R2, R2_cp);

	nmod_mat_window_clear(R1);
	nmod_mat_window_clear(R2);
	}

	/** apply perm^(-1) **/
	_perm_inv(inv_perm, perm, r);
	for (i = deb + 1; i < sigma; i++)
		apply_perm_rows_to_matrix(res->mat + i, inv_perm, r);

	/** clear **/
	nmod_mat_clear(R1_cp);
	nmod_mat_clear(R2_cp);
	nmod_mat_clear(previous_R1);

	_perm_clear(inv_perm);
}

void structured_list_multiplication_blocks_full(nmod_mat_poly_t res, const nmod_mat_t A,
		const slong *perm, slong rank)
{
	slong r = res->r, c = res->c;
	mp_limb_t modulus = res->modulus;
	slong i;
	nmod_mat_t previous_R1, R1, R2, R1_cp, R2_cp;
	nmod_mat_struct *r_i;
	slong *inv_perm = _perm_init(r);

	/** init **/
	nmod_mat_init(R1_cp, rank, c, modulus);
	nmod_mat_init(previous_R1, rank, c, modulus);
	nmod_mat_init(R2_cp, r - rank, c, modulus);

	res->degree += 1;

	for (i = 0; i <= res->degree; i++)
	{
		r_i = res->mat + i;
		apply_perm_rows_to_matrix(r_i, perm, r);
		nmod_mat_window_init(R1, r_i, 0, 0, rank, c);
		nmod_mat_set(R1_cp, R1);

		nmod_mat_set(R1, previous_R1); //it will set the top block of r_{i-1} on r_{i}

	nmod_mat_set(previous_R1, R1_cp);

	nmod_mat_window_init(R2, r_i, rank, 0, r, c);
	nmod_mat_set(R2_cp, R2);

	nmod_mat_mul(R2, A, R1_cp);
	nmod_mat_add(R2, R2, R2_cp);

	nmod_mat_window_clear(R1);
	nmod_mat_window_clear(R2);
	}

	/** apply perm^(-1) **/
	_perm_inv(inv_perm, perm, r);
	for (i = 0; i <= res->degree; i++)
		apply_perm_rows_to_matrix(res->mat + i, inv_perm, r);

	/** clear **/
	nmod_mat_clear(R1_cp);
	nmod_mat_clear(R2_cp);
	nmod_mat_clear(previous_R1);

	_perm_clear(inv_perm);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
