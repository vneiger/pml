#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

int main()
{
	// take example from SageMath doc
	// Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
	nmod_poly_mat_t mat;
	nmod_poly_mat_init(mat, 2, 3, 7);

	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);

    printf("Input matrix (flint printing) is:\n");
    nmod_poly_mat_print(mat, "x");


    printf("With SageMath printing, this is:\n");
    nmod_poly_mat_print_sagemath(mat, "x");

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
