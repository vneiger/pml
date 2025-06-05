#include <flint/nmod_poly.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // for print_pretty

int main()
{
    // take example from SageMath doc
    // Matrix(pR, [ [3*x+1, 0, 1], [x^3+3, 0, 0] ])
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, 2, 3, 7);

    slong rshift[] = {0,5,2};
    slong cshift[] = {1,2};

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 3, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    // initialize pivot index and pivot degree
    slong rpivind[2], cpivind[3];
    slong rpivdeg[2], cpivdeg[3];

    // sage: M.leading_positions()
    // [0, 0]
    nmod_poly_mat_pivot_profile(rpivind, rpivdeg, mat, NULL, ROW_LOWER);
    printf("Checking row-wise pivot profile [%ld, %ld], [%ld, %ld] --> %s\n",
           rpivind[0], rpivind[1],
           rpivdeg[0], rpivdeg[1],
           (rpivind[0] == 0 && rpivind[1] == 0
            && rpivdeg[0] == 1 && rpivdeg[1] == 3) ? "ok" : "not ok");

    // sage: M.leading_positions(shifts=[0,5,2], return_degree=True)
    // ([2, 0], [0, 3])
    nmod_poly_mat_pivot_profile(rpivind, rpivdeg, mat, rshift, ROW_LOWER);
    printf("Checking (0,5,2)-shifted row-wise pivot profile [%ld, %ld], [%ld, %ld] --> %s\n",
           rpivind[0], rpivind[1],
           rpivdeg[0], rpivdeg[1],
           (rpivind[0] == 2 && rpivind[1] == 0
            && rpivdeg[0] == 0 && rpivdeg[1] == 3) ? "ok" : "not ok");

    // sage: M.leading_positions(row_wise=False, return_degree=True)
    // ([1, -1, 0], [3, -1, 0])
    nmod_poly_mat_pivot_profile(cpivind, cpivdeg, mat, NULL, COL_UPPER);
    printf("Checking column-wise pivot profile [%ld, %ld, %ld], [%ld, %ld, %ld] --> %s\n",
           cpivind[0], cpivind[1], cpivind[2],
           cpivdeg[0], cpivdeg[1], cpivdeg[2],
           (cpivind[0] == 1 && cpivind[1] == -1 && cpivind[2] == 0
            && cpivdeg[0] == 3 && cpivdeg[1] == -1 && cpivdeg[2] == 0) ? "ok" : "not ok");

    // sage: M.leading_positions(shifts=[1,2], row_wise=False,
    //                           ....:   return_degree=True)
    // ([1, -1, 0], [3, -1, 0])
    nmod_poly_mat_pivot_profile(cpivind, cpivdeg, mat, cshift, COL_UPPER);
    printf("Checking (1,2)-shifted column-wise pivot profile [%ld, %ld, %ld], [%ld, %ld, %ld] --> %s\n",
           cpivind[0], cpivind[1], cpivind[2],
           cpivdeg[0], cpivdeg[1], cpivdeg[2],
           (cpivind[0] == 1 && cpivind[1] == -1 && cpivind[2] == 0
            && cpivdeg[0] == 3 && cpivdeg[1] == -1 && cpivdeg[2] == 0) ? "ok" : "not ok");

    // In case several entries in the row (resp. column) reach the shifted row
    // (resp. column) degree, the leading position is chosen as the rightmost
    // (resp. bottommost) such entry::

    // sage: M.leading_positions(shifts=[0,5,1],return_degree=True)
    // ([2, 0], [0, 3])
    rshift[2] = 1;
    nmod_poly_mat_pivot_profile(rpivind, rpivdeg, mat, rshift, ROW_LOWER);
    printf("Checking (0,5,1)-shifted row-wise pivot profile [%ld, %ld], [%ld, %ld] --> %s\n",
           rpivind[0], rpivind[1],
           rpivdeg[0], rpivdeg[1],
           (rpivind[0] == 2 && rpivind[1] == 0
            && rpivdeg[0] == 0 && rpivdeg[1] == 3) ? "ok" : "not ok");

    // sage: M.leading_positions(shifts=[2,0], row_wise=False,return_degree=True)
    // ([1, -1, 0], [3, -1, 0])
    cshift[0]=2; cshift[1]=0;
    nmod_poly_mat_pivot_profile(cpivind, cpivdeg, mat, cshift, COL_UPPER);
    printf("Checking (2,0)-shifted column-wise pivot profile [%ld, %ld, %ld], [%ld, %ld, %ld] --> %s\n",
           cpivind[0], cpivind[1], cpivind[2],
           cpivdeg[0], cpivdeg[1], cpivdeg[2],
           (cpivind[0] == 1 && cpivind[1] == -1 && cpivind[2] == 0
            && cpivdeg[0] == 3 && cpivdeg[1] == -1 && cpivdeg[2] == 0) ? "ok" : "not ok");

    nmod_poly_mat_clear(mat);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
