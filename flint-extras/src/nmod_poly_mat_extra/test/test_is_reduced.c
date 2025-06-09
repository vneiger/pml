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
    slong cshift[] = {2,0};

	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    int isred;

    // sage: M.is_reduced()
    // False
    isred = nmod_poly_mat_is_reduced(mat, NULL, ROW_LOWER);
    printf("Checking is reduced row-wise: %d --> %s\n",
           isred, (isred==0) ? "ok" : "not ok");

    // sage: M.is_reduced(shifts=[0,1,2])
    // True
    isred = nmod_poly_mat_is_reduced(mat, rshift, ROW_LOWER);
    printf("Checking is (0,1,2)-reduced row-wise: %d --> %s\n",
           isred, (isred==1) ? "ok" : "not ok");

    isred = nmod_poly_mat_is_reduced(mat, NULL, COL_UPPER);
    printf("Checking is reduced column-wise: %d --> %s\n",
           isred, (isred==0) ? "ok" : "not ok");

    // sage: M.is_reduced(shifts=[2,0], row_wise=False,
    // ....:                           include_zero_vectors=False)
    // False
    isred = nmod_poly_mat_is_reduced(mat, cshift, COL_UPPER);
    printf("Checking is (2,0)-reduced column-wise: %d --> %s\n",
           isred, (isred==0) ? "ok" : "not ok");

    // sage: M = Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0], [0, 1, 0] ])
    // sage: M.is_reduced(shifts=[2,0,0], row_wise=False)
    // True
    nmod_poly_mat_clear(mat);
    nmod_poly_mat_init(mat, 3, 3, 7);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 0, 1);
    printf("\nInput matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    isred = nmod_poly_mat_is_reduced(mat, NULL, COL_UPPER);
    printf("Checking is reduced column-wise: %d --> %s\n",
           isred, (isred==1) ? "ok" : "not ok");

    isred = nmod_poly_mat_is_reduced(mat, cshift, COL_UPPER);
    printf("Checking is (2,0)-reduced column-wise: %d --> %s\n",
           isred, (isred==1) ? "ok" : "not ok");

    nmod_poly_mat_clear(mat);

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
