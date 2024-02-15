#include <flint/nmod_poly.h>
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_forms.h"

int main()
{
    // take example from SageMath doc
    // sage: pR.<x> = GF(7)[]
    // sage: M = Matrix([ [x^3+3*x^2+6*x+6, 3*x^2+3*x+6, 4*x^2+x+3],
    // ....:              [5,               1,           0        ],
    // ....:              [2*x^2+2,         2*x+5,       x^2+4*x+6] ])
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, 3, 3, 7);

    slong shift[3];

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 2, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 3, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 1, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 2, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 1, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 2, 4);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 0, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 0, 2);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 2, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 1, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 1, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 2, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    int iswp;

    // sage: M.is_weak_popov()
    // True
    iswp = nmod_poly_mat_is_weak_popov(mat, NULL, ROW_LOWER);
    printf("Checking is weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==1) ? "ok" : "not ok");

    // One can check whether the leading positions, in addition to being
    // pairwise distinct, are actually in increasing order::

    // sage: M.is_weak_popov(ordered=True)
    // True
    iswp = nmod_poly_mat_is_ordered_weak_popov(mat, NULL, ROW_LOWER);
    printf("Checking is ordered weak Popov row-wise: %d --> %s\n\n",
           iswp, (iswp==1) ? "ok" : "not ok");

    // sage: N = M.with_swapped_rows(1,2)
    // sage: N.is_weak_popov()
    // True
    // sage: N.is_weak_popov(ordered=True)
    // False
    nmod_poly_mat_t mat2;
    nmod_poly_mat_init_set(mat2, mat);
    nmod_poly_struct * tmp = mat2->rows[0];
    mat2->rows[0] = mat2->rows[1];
    mat2->rows[1] = tmp;
    printf("Consider input matrix with swapped rows:\n");
    nmod_poly_mat_print_pretty(mat2, "x");
    printf("\n");

    iswp = nmod_poly_mat_is_weak_popov(mat2, NULL, ROW_LOWER);
    printf("Checking is weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==1) ? "ok" : "not ok");
    iswp = nmod_poly_mat_is_ordered_weak_popov(mat2, NULL, ROW_LOWER);
    printf("Checking is ordered weak Popov row-wise: %d --> %s\n\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    nmod_poly_mat_clear(mat2);

    printf("Back to the original matrix:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    // sage: M.is_weak_popov(shifts=[2,3,1])
    // False
    shift[0] = 2; shift[1] = 3; shift[2] = 1;
    iswp = nmod_poly_mat_is_weak_popov(mat, shift, ROW_LOWER);
    printf("Checking is (2,3,1)-weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    // sage: M.is_weak_popov(shifts=[0,2,0],row_wise=False,ordered=True)
    // True
    shift[0] = 0; shift[1] = 2; shift[2] = 0;
    iswp = nmod_poly_mat_is_ordered_weak_popov(mat, shift, COL_UPPER);
    printf("Checking is (0,2,0)-ordered weak Popov column-wise: %d --> %s\n\n",
           iswp, (iswp==1) ? "ok" : "not ok");
    // TODO

    // Rectangular cases:
    // sage: M = Matrix([
    // ....:    [  x^3+5*x^2+5*x+1,       5,       6*x+4,         0],
    // ....:    [      6*x^2+3*x+1,       1,           2,         0],
    // ....:    [2*x^3+4*x^2+6*x+4, 5*x + 1, 2*x^2+5*x+5, x^2+5*x+6]
    // ....:     ])
    nmod_poly_mat_clear(mat);
    nmod_poly_mat_init(mat, 3, 4, 7);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 2, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 3, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 0, 5);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 1, 6);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 1, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 2, 6);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 0, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 2), 0, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 0, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 1, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 2, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 3, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 1, 5);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 2, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 3), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 3), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 3), 2, 1);

    printf("Now a rectangular matrix:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    iswp = nmod_poly_mat_is_weak_popov(mat, NULL, ROW_LOWER);
    printf("Checking is weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    iswp = nmod_poly_mat_is_ordered_weak_popov(mat, NULL, ROW_LOWER);
    printf("Checking is ordered weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    slong shift2[4] = {0,3,2,2};
    iswp = nmod_poly_mat_is_weak_popov(mat, shift2, ROW_LOWER);
    printf("Checking is (0,3,2,2)-weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==1) ? "ok" : "not ok");
    // TODO (todo what??)

    iswp = nmod_poly_mat_is_ordered_weak_popov(mat, shift2, ROW_LOWER);
    printf("Checking is (0,3,2,2)-ordered weak Popov row-wise: %d --> %s\n\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    // sage: M.is_weak_popov(shifts=[0,2,1,3])
    // True
    // sage: M.is_weak_popov(shifts=[0,2,1,3],ordered=True)
    // True
    slong shift3[4] = {0,2,1,3};
    iswp = nmod_poly_mat_is_weak_popov(mat, shift3, ROW_LOWER);
    printf("Checking is (0,2,1,3)-weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==1) ? "ok" : "not ok");

    iswp = nmod_poly_mat_is_ordered_weak_popov(mat, shift3, ROW_LOWER);
    printf("Checking is (0,2,1,3)-ordered weak Popov row-wise: %d --> %s\n\n",
           iswp, (iswp==1) ? "ok" : "not ok");

    // Zero rows (resp. columns if column-wise) cannot yield weak Popov matrices:
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 0));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 1));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 2));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 3));

    printf("Finally same matrix with zero second row:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    iswp = nmod_poly_mat_is_weak_popov(mat, NULL, ROW_LOWER);
    printf("Checking is weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    iswp = nmod_poly_mat_is_weak_popov(mat, shift3, ROW_LOWER);
    printf("Checking is (0,2,1,3)-weak Popov row-wise: %d --> %s\n",
           iswp, (iswp==0) ? "ok" : "not ok");

    iswp = nmod_poly_mat_is_weak_popov(mat, NULL, COL_UPPER);
    printf("Checking is weak Popov column-wise: %d --> %s\n",
           iswp, (iswp==0) ? "ok" : "not ok");


    nmod_poly_mat_clear(mat);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
