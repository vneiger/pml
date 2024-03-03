#include <flint/nmod_poly.h>
#include "nmod_poly_mat_io.h"
#include "nmod_poly_mat_forms.h"

int main()
{
    // take example from SageMath doc
    // pR.<x> = GF(7)[]
    // M = Matrix(pR, [ [x^4+6*x^3+4*x+4, 3*x+6,     3  ],
    //                  [0,               x^2+5*x+5, 2  ],
    //                  [0,               0,         x+5] ])
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, 3, 3, 7);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 0, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 1, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 3, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 0), 4, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 1), 1, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 0, 2), 0, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 1), 2, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 1, 2), 0, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat, 2, 2), 1, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat, "x");
    printf("\n");

    int ish;


    // M.is_hermite()
    // True
    ish = nmod_poly_mat_is_hermite(mat, ROW_UPPER); // sage's default is upper
    printf("Checking is upper Hermite row-wise: %d --> %s\n",
           ish, (ish==1) ? "ok" : "not ok");

    // M.is_hermite(row_wise=False)
    // True
    ish = nmod_poly_mat_is_hermite(mat, COL_UPPER); // sage's default is upper
    printf("Checking is upper Hermite column-wise: %d --> %s\n",
           ish, (ish==1) ? "ok" : "not ok");

    // M.is_hermite(row_wise=False, lower_echelon=True)
    // False
    ish = nmod_poly_mat_is_hermite(mat, ROW_LOWER);
    printf("Checking is lower Hermite row-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");

    // N = Matrix(pR, [ [x+5, 0,               0        ],
    //                  [2,   x^4+6*x^3+4*x+4, 0        ],
    //                  [3,   3*x^3+6,         x^2+5*x+5] ])
    nmod_poly_mat_t mat2;
    nmod_poly_mat_init(mat2, 3, 3, 7);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 0, 0), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 0, 0), 1, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 0), 0, 2);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 1), 0, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 1), 1, 4);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 1), 3, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 1, 1), 4, 1);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 2, 0), 0, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 2, 1), 0, 6);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 2, 1), 3, 3);

    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 2, 2), 0, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 2, 2), 1, 5);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(mat2, 2, 2), 2, 1);

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat2, "x");
    printf("\n");

    // N.is_hermite()
    // False
    ish = nmod_poly_mat_is_hermite(mat2, ROW_UPPER);
    printf("Checking is upper Hermite row-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");

    // N.is_hermite(lower_echelon=True)
    // True
    ish = nmod_poly_mat_is_hermite(mat2, ROW_LOWER);
    printf("Checking is lower Hermite row-wise: %d --> %s\n",
           ish, (ish==1) ? "ok" : "not ok");

    // N.is_hermite(row_wise=False)
    // False
    ish = nmod_poly_mat_is_hermite(mat2, COL_UPPER);
    printf("Checking is upper Hermite column-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");

    // N.is_hermite(row_wise=False, lower_echelon=True)
    // False 
    ish = nmod_poly_mat_is_hermite(mat2, COL_LOWER);
    printf("Checking is lower Hermite row-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");

    // Rectangular matrices with zero rows are supported. Zero rows (resp. columns) can be forbidden, and otherwise they should be at the bottom (resp. the right-hand side) of the matrix:
    // NOTE: zero vectors not accepted in FLINT

    // N[:,1:].is_hermite(lower_echelon=True)
    // False
    nmod_poly_mat_t mat3;
    nmod_poly_mat_init(mat3, 3, 2, 7);
    nmod_poly_set(nmod_poly_mat_entry(mat3, 0, 0), nmod_poly_mat_entry(mat2, 0, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 1, 0), nmod_poly_mat_entry(mat2, 1, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 2, 0), nmod_poly_mat_entry(mat2, 2, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 0, 1), nmod_poly_mat_entry(mat2, 0, 2));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 1, 1), nmod_poly_mat_entry(mat2, 1, 2));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 2, 1), nmod_poly_mat_entry(mat2, 2, 2));

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat3, "x");
    printf("\n");

    ish = nmod_poly_mat_is_hermite(mat3, ROW_LOWER);
    printf("Checking is lower Hermite row-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");

    ish = nmod_poly_mat_is_echelon(mat3, COL_LOWER);
    printf("Checking is lower echelon column-wise: %d --> %s\n",
           ish, (ish==1) ? "ok" : "not ok");

    ish = nmod_poly_mat_is_hermite(mat3, COL_LOWER);
    printf("Checking is lower Hermite column-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");

    // N[[1,2,0],1:].is_hermite(lower_echelon=True)
    // True
    // --> let's not try this one here, zero row again

    // N[:2,:].is_hermite(row_wise=False, lower_echelon=True)
    // True
    nmod_poly_mat_clear(mat3);
    nmod_poly_mat_init(mat3, 2, 3, 7);
    nmod_poly_set(nmod_poly_mat_entry(mat3, 0, 0), nmod_poly_mat_entry(mat2, 0, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 0, 1), nmod_poly_mat_entry(mat2, 0, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 0, 2), nmod_poly_mat_entry(mat2, 0, 2));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 1, 0), nmod_poly_mat_entry(mat2, 1, 0));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 1, 1), nmod_poly_mat_entry(mat2, 1, 1));
    nmod_poly_set(nmod_poly_mat_entry(mat3, 1, 2), nmod_poly_mat_entry(mat2, 1, 2));

    printf("Input matrix is:\n");
    nmod_poly_mat_print_pretty(mat3, "x");
    printf("\n");

    ish = nmod_poly_mat_is_hermite(mat3, ROW_LOWER);
    printf("Checking is lower Hermite row-wise: %d --> %s\n",
           ish, (ish==1) ? "ok" : "not ok");

    ish = nmod_poly_mat_is_echelon(mat3, COL_LOWER);
    printf("Checking is lower echelon column-wise: %d --> %s\n",
           ish, (ish==1) ? "ok" : "not ok");

    ish = nmod_poly_mat_is_hermite(mat3, COL_LOWER);
    printf("Checking is lower Hermite column-wise: %d --> %s\n",
           ish, (ish==0) ? "ok" : "not ok");


    // N[:2,:].is_hermite(row_wise=False,
    //                    lower_echelon=True,
    //                    include_zero_vectors=False)
    // False
    // --> let's not try this one here, no corresponding option in flint

    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(mat2);
    nmod_poly_mat_clear(mat3);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
