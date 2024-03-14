#include <flint/nmod_poly.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h"

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

    // initialize rdeg
    slong rdeg[2], cdeg[3];

    // sage: M.row_degrees()
    // [1, 3]
    nmod_poly_mat_row_degree(rdeg, mat, NULL);
    printf("Checking row degree: [%ld, %ld] --> %s\n",
           rdeg[0], rdeg[1], (rdeg[0] == 1 && rdeg[1] == 3) ? "ok" : "not ok");

    // sage: M.row_degrees(shifts=[0,1,2])
    // [2, 3]
    nmod_poly_mat_row_degree(rdeg, mat, rshift);
    printf("Checking (0,1,2)-shifted row degree: [%ld, %ld] --> %s\n",
           rdeg[0], rdeg[1], (rdeg[0] == 2 && rdeg[1] == 3) ? "ok" : "not ok");

    // sage: M.column_degrees()
    // [3, -1, 0]
    nmod_poly_mat_column_degree(cdeg, mat, NULL);
    printf("Checking column degree: [%ld, %ld, %ld] --> %s\n",
           cdeg[0], cdeg[1], cdeg[2],
           (cdeg[0] == 3 && cdeg[1] == -1 && cdeg[2] == 0) ? "ok" : "not ok");

    // sage: M.column_degrees(shifts=[0,2])
    // [5, -1, 0]
    nmod_poly_mat_column_degree(cdeg, mat, cshift);
    printf("Checking (-2,1)-shifted column degree: [%ld, %ld, %ld] --> %s\n",
           cdeg[0], cdeg[1], cdeg[2],
           (cdeg[0] == 4 && cdeg[1] == -3 && cdeg[2] == -2) ? "ok" : "not ok");

    nmod_poly_mat_clear(mat);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
