#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"

int main()
{
    // take example from SageMath doc
    //pR.<x> = GF(7)[]
    //M = Matrix(pR, [ [x^4+6*x^3+4*x+4, 3*x+6,     3  ],
    //                 [x^2+6*x+6,       x^2+5*x+5, 2  ],
    //                 [3*x,             6*x+5,     x+5] ])
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, 3, 3, 7);

    slong shift[3];

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 3, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 4, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 1, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 1, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 2, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 2, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 2), 0, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 1, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 1), 1, 6);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 1, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    int isp;

    //M.is_popov()
    //True
    isp = nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER);
    printf("Checking is Popov row-wise: %d --> %s\n",
           isp, (isp==1) ? "ok" : "not ok");

    //M.is_popov(shifts=[0,1,2])
    //True
    shift[0] = 0; shift[1] = 1; shift[2] = 2;
    isp = nmod_poly_mat_is_popov(mat, shift, ROW_LOWER);
    printf("Checking is (0,1,2)-Popov row-wise: %d --> %s\n",
           isp, (isp==1) ? "ok" : "not ok");

    //M[:,:2].is_popov()
    //False
    nmod_poly_mat_t mat2;
    nmod_poly_mat_init(mat2, 3, 2, mat->modulus);
    nmod_poly_set(nmod_poly_mat_entry(mat2, 0, 0), nmod_poly_mat_entry(mat, 0, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 0, 1), nmod_poly_mat_entry(mat, 0, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 1, 0), nmod_poly_mat_entry(mat, 1, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 1, 1), nmod_poly_mat_entry(mat, 1, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 2, 0), nmod_poly_mat_entry(mat, 2, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 2, 1), nmod_poly_mat_entry(mat, 2, 1));

    isp = nmod_poly_mat_is_popov(mat2, NULL, ROW_LOWER);
    printf("Checking if :,0:2 submatrix is Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    //M[:2,:].is_popov(shifts=[0,1,2])
    //True
    nmod_poly_mat_clear(mat2);
    nmod_poly_mat_init(mat2, 2, 3, mat->modulus);
    nmod_poly_set(nmod_poly_mat_entry(mat2, 0, 0), nmod_poly_mat_entry(mat, 0, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 0, 1), nmod_poly_mat_entry(mat, 0, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 0, 2), nmod_poly_mat_entry(mat, 0, 2));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 1, 0), nmod_poly_mat_entry(mat, 1, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 1, 1), nmod_poly_mat_entry(mat, 1, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat2, 1, 2), nmod_poly_mat_entry(mat, 1, 2));

    isp = nmod_poly_mat_is_popov(mat2, shift, ROW_LOWER);
    printf("Checking if 0:2,: submatrix is (0,1,2)-Popov row-wise: %d --> %s\n",
           isp, (isp==1) ? "ok" : "not ok");

    //M = Matrix(pR, [ [x^4+3*x^3+x^2+2*x+6, x^3+5*x^2+5*x+1],
    //                 [6*x+1,               x^2+4*x+1      ],
    //                 [6,                   6              ] ])
    nmod_poly_mat_clear(mat);
    nmod_poly_mat_init(mat, 3, 2, 7);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 2);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 2, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 3, 3);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 4, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 2, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 3, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 0), 1, 6);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 1, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 2, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 0, 6);

    printf("\n");
    printf("Now, input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");

    //M.is_popov(row_wise=False)
    //False
    isp = nmod_poly_mat_is_popov(mat, NULL, COL_UPPER);
    printf("Checking is Popov column-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    //M.is_popov(shifts=[0,2,3], row_wise=False)
    //True
    shift[0] = 0; shift[1] = 2; shift[2] = 3;
    isp = nmod_poly_mat_is_popov(mat, shift, COL_UPPER);
    printf("Checking is (0,2,3)-Popov column-wise: %d --> %s\n",
           isp, (isp==1) ? "ok" : "not ok");

    //One can forbid zero rows (or columns if not working row-wise):
    // --> NOTE in flint zero rows ARE forbidden always
    //
    //N = Matrix(pR, [ [x^4+3*x^3+x^2+2*x+6, 6*x+1     ],
    //                 [5*x^2+5*x+1,         x^2+4*x+1 ],
    //                 [0,                   0         ] ])
    nmod_poly_mat_init_set(mat2, mat);
    nmod_poly_zero(nmod_poly_mat_entry(mat2, 2, 0));
    nmod_poly_zero(nmod_poly_mat_entry(mat2, 2, 1));
    nmod_poly_zero(nmod_poly_mat_entry(mat2, 0, 1));
    nmod_poly_zero(nmod_poly_mat_entry(mat2, 1, 0));

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 0, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 0, 1), 1, 6);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 0), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 0), 2, 5);

    //N.is_popov()
    //True
    //N.is_popov(include_zero_vectors=False)
    //False
    isp = nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER);
    printf("Checking is Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok"); // here NOT Popov (there is zero row)

    //One can verify Popov form up to row permutation (or column permutation if not working row-wise):
    //M.swap_columns(0,1)
    //M.is_popov(shifts=[0,2,3], row_wise=False)
    //False
    nmod_poly_mat_swap_columns(mat, NULL, 0, 1);
    isp = nmod_poly_mat_is_popov(mat, shift, ROW_LOWER);
    printf("Checking on swapped columns if is (0,2,3)-Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    //M.is_popov(shifts=[0,2,3], row_wise=False, up_to_permutation=True)
    //True
    // TODO up_to_permutation Not supported in Flint at the moment

    //N.swap_rows(0,2)
    //N.is_popov()
    //False
    nmod_poly_mat_swap_rows(mat2, NULL, 0, 2);
    isp = nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER);
    printf("Checking on swapped rows if is Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    //N.is_popov(up_to_permutation=True)
    //True
    // TODO up_to_permutation Not supported in Flint at the moment

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

    printf("\n");
    printf("Now a rectangular matrix:\n");
    nmod_poly_mat_print_pretty(mat, "x");

    isp = nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER);
    printf("Checking is Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    slong shift2[4] = {0,3,2,2};
    isp = nmod_poly_mat_is_popov(mat, shift2, ROW_LOWER);
    printf("Checking is (0,3,2,2)-Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    slong shift3[4] = {0,2,1,3};
    isp = nmod_poly_mat_is_popov(mat, shift3, ROW_LOWER);
    printf("Checking is (0,2,1,3)-Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    printf("-> it has non-reduced entries above/below pivots; let's reduce them\n");
    nmod_poly_zero(nmod_poly_mat_entry(mat, 0, 1));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 2, 1));
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 0), 3, 0);

    isp = nmod_poly_mat_is_popov(mat, shift3, ROW_LOWER);
    printf("Checking is (0,2,1,3)-Popov row-wise: %d --> %s\n",
           isp, (isp==1) ? "ok" : "not ok");


    // Zero rows (resp. columns if column-wise) cannot yield Popov matrices:
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 0));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 1));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 2));
    nmod_poly_zero(nmod_poly_mat_entry(mat, 1, 3));

    printf("\n");
    printf("Finally same matrix with zero second row:\n");
    nmod_poly_mat_print_pretty(mat, "x");

    isp = nmod_poly_mat_is_popov(mat, NULL, ROW_LOWER);
    printf("Checking is Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    isp = nmod_poly_mat_is_popov(mat, shift3, ROW_LOWER);
    printf("Checking is (0,2,1,3)-Popov row-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    isp = nmod_poly_mat_is_popov(mat, NULL, COL_UPPER);
    printf("Checking is Popov column-wise: %d --> %s\n",
           isp, (isp==0) ? "ok" : "not ok");

    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(mat2);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
