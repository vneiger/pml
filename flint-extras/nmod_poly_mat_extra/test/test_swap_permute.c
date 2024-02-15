#include <flint/nmod_poly.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"

int main()
{
	// build [ [3*x+1, 0, 1], [x^3+3, 0, 0] ]
	nmod_poly_mat_t mat;
	nmod_poly_mat_init(mat, 2, 3, 7);

	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);

    // build permutations
    slong perm_act_rows[] = {0,1};
    slong perm_store_rows[] = {0,1};
    slong perm_act_cols[] = {2,0,1};
    slong perm_store_cols[] = {0,1,2};

    printf("Initial perm_store_rows is [%ld,%ld]\n",
           perm_store_rows[0], perm_store_rows[1]);
    printf("Initial perm_act_rows is [%ld,%ld]\n",
           perm_act_rows[0], perm_act_rows[1]);

    printf("Initial perm_store_cols is [%ld,%ld,%ld]\n",
           perm_store_cols[0], perm_store_cols[1], perm_store_cols[2]);
    printf("Initial perm_act_cols is [%ld,%ld,%ld]\n",
           perm_act_cols[0], perm_act_cols[1], perm_act_cols[2]);


    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    nmod_poly_mat_swap_rows(mat, perm_store_rows, 0, 1);
    printf("After swapping rows 0-1, perm_store_rows is [%ld,%ld] and matrix is\n",
           perm_store_rows[0], perm_store_rows[1]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    nmod_poly_mat_swap_rows(mat, perm_store_rows, 0, 0);
    printf("After swapping rows 0-0, perm_store_rows is [%ld,%ld] and matrix is\n",
           perm_store_rows[0], perm_store_rows[1]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    nmod_poly_mat_swap_columns(mat, perm_store_cols, 1, 2);
    printf("After swapping coluns 1-2, perm_store_cols is [%ld,%ld,%ld] and matrix is\n",
           perm_store_cols[0], perm_store_cols[1], perm_store_cols[2]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");


    nmod_poly_mat_invert_rows(mat, perm_store_rows);
    printf("After inverting rows, perm_store_rows is [%ld,%ld] and matrix is\n",
           perm_store_rows[0], perm_store_rows[1]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    nmod_poly_mat_invert_columns(mat, perm_store_cols);
    printf("After inverting columns, perm_store_cols is [%ld,%ld,%ld] and matrix is\n",
           perm_store_cols[0], perm_store_cols[1], perm_store_cols[2]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");


    nmod_poly_mat_permute_rows(mat, perm_act_rows, perm_store_rows);
    printf("After permuting rows with permutation [0,1], perm_store_rows is [%ld,%ld] and matrix is\n",
           perm_store_rows[0], perm_store_rows[1]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    perm_act_rows[0] = 1; perm_act_rows[1] = 0; 
    nmod_poly_mat_permute_rows(mat, perm_act_rows, perm_store_rows);
    printf("After permuting rows with permutation [1,0], perm_store_rows is [%ld,%ld] and matrix is\n",
           perm_store_rows[0], perm_store_rows[1]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    nmod_poly_mat_permute_columns(mat, perm_act_cols, perm_store_cols);
    printf("After permuting columns with permutation [2,0,1], perm_store_cols is [%ld,%ld,%ld] and matrix is\n",
           perm_store_cols[0], perm_store_cols[1], perm_store_cols[2]);
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");


    // matrix [ [3*x+1, 0, 1], [x^3+3, 0, 0], [0, 1, 0] ]
    //nmod_poly_mat_clear(mat);
    //nmod_poly_mat_init(mat, 3, 3, 7);
	//nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
	//nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
	//nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
	//nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
	//nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);
	//nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 0, 1);
    //printf("\nInput matrix is:\n");
    //nmod_poly_mat_print_pretty(mat, "x");
    //printf("\n");

    nmod_poly_mat_clear(mat);

	return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
