#include <flint/nmod_poly.h>
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_forms.h"

int main()
{
	// take example from SageMath doc
	// Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
	nmod_poly_mat_t mat;
	nmod_poly_mat_init(mat, 2, 3, 7);

    slong rshift[] = {0,1,2};
    slong cshift[] = {-2,1};

	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    // sage: M.nmod_poly_mat_degree_matrix()
    // [ 1 -1  0]
    // [ 3 -1 -1]
    printf("Printing degree matrix:\n");
    nmod_poly_mat_degree_matrix_print_pretty(mat);
    printf("\n");

    // sage: M.nmod_poly_mat_degree_matrix(shifts=[0,1,2])
    // [ 1 -1  2]
    // [ 3 -1 -1]
    // Different convention for zero entries --> here
    // [ 1 0  2]
    // [ 3 0  1]
    printf("Printing row-wise [0,1,2]-shifted degree matrix:\n");
    nmod_poly_mat_degree_matrix_shifted_print_pretty(mat, rshift, ROW_LOWER);
    printf("\n");

    // sage: M.nmod_poly_mat_degree_matrix(shifts=[-2,1,2])
    // [-1 -3  2]
    // [ 1 -3 -3]
    // Different convention for zero entries --> here
    // [-1 0 2]
    // [ 1 0 1]
    rshift[0] = -2; rshift[1] = 1; rshift[2] = 2;
    printf("Printing row-wise [-2,1,2]-shifted degree matrix:\n");
    nmod_poly_mat_degree_matrix_shifted_print_pretty(mat, rshift, ROW_LOWER);
    printf("\n");

    // sage: M.nmod_poly_mat_degree_matrix(shifts=[-1,2], row_wise=False)
    // [ 0 -2 -1]
    // [ 5 -2 -2]
    // Different convention for zero entries --> here
    // [ 0 -2 -1]
    // [ 5  1  1]
    cshift[0] = -1; cshift[1] = 2;
    printf("Printing column-wise [-1,2]-shifted degree matrix:\n");
    nmod_poly_mat_degree_matrix_shifted_print_pretty(mat, cshift, COL_UPPER);
    printf("\n");

    nmod_poly_mat_clear(mat);

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
