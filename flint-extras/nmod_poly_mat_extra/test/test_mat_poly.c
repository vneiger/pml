#include <time.h>
#include <stdlib.h>
#include <flint/nmod_types.h>
#include <flint/fmpz_mat.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

// TODO make random choice for given prime (or for random prime of given size)
#define PRIME_30_BITS 536870923
#define PRIME_60_BITS 576460752303423619

/** Test
 * Verify if nmod_poly_mat_mul + coefficient_matrix and
 * nmod_mat_poly0_naive_mul_coef
 * return the same result.
 * Verify if nmod_poly_mat -> nmod_mat_poly0 -> nmod_poly_mat
 * unchanged the initial matrix
 */
int test_nmod_mat_poly0(void)
{
    nmod_mat_t prod, prod_2;
    nmod_poly_mat_t mat, mat_2, res, res_2, poly_prod;
    nmod_mat_poly0_t mat_repr, mat_repr_2;
    slong rdim = 6, cdim = 2, prime = 7,
          len = 10, len_2 = 20, k = 28;
    char *x = malloc(1);

    /** init **/
    nmod_poly_mat_init(mat, rdim, cdim, prime);
    nmod_poly_mat_init(mat_2, cdim, rdim, prime);

    nmod_poly_mat_init(res, rdim, cdim, prime);
    nmod_poly_mat_init(res_2, cdim, rdim, prime);

    nmod_poly_mat_init(poly_prod, rdim, rdim, prime);

    nmod_mat_init(prod, rdim, rdim, prime);
    nmod_mat_init(prod_2, rdim, rdim, prime);

    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    /** test poly_mat <-> mat_poly0 **/
    nmod_poly_mat_randtest(mat, state, len);
    nmod_poly_mat_randtest(mat_2, state, len_2);

    nmod_mat_poly0_init_setII(mat_repr, mat, 10);
    nmod_mat_poly0_init_setII(mat_repr_2, mat_2, 9);

    nmod_poly_mat_set_from_mat_poly0(res, mat_repr);
    nmod_poly_mat_set_from_mat_poly0(res_2, mat_repr_2);

    if (!nmod_poly_mat_equal(mat, res))
        printf("error mat not equals");
    if (!nmod_poly_mat_equal(mat_2, res_2))
        printf("\nerror mat_2 not equals\n");

    /** test mul: poly_mat <-> mat_poly0 **/
    nmod_poly_mat_mul(poly_prod, res, res_2);
    nmod_poly_mat_get_coeff_mat(prod, poly_prod, k);

    nmod_poly_mat_print(poly_prod, x);

    nmod_mat_poly0_naive_mul_coef(prod_2, mat_repr, mat_repr_2, k);


    //print
    nmod_mat_print_pretty(prod);
    nmod_mat_print_pretty(prod_2);

    if (!nmod_mat_equal(prod_2, prod))
        printf("error mul not equals");

    /** clear **/
    nmod_mat_poly0_clear(mat_repr);
    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(res);

    nmod_mat_poly0_clear(mat_repr_2);
    nmod_poly_mat_clear(mat_2);
    nmod_poly_mat_clear(res_2);

    nmod_poly_mat_clear(poly_prod);
    nmod_mat_clear(prod);
    nmod_mat_clear(prod_2);

    flint_rand_clear(state);
    free(x);

    return 0;
}

/** Test
 * many function of matpol, verify if there are leaks
 **/
int test_matpol(void)
{
    nmod_poly_mat_t A;
    nmod_poly_mat_init(A,3,3,4);

    slong cdim = A->c, rdim = A->r;
    printf("number of columns: %lu, and rows: %lu\n", cdim, rdim);

    flint_rand_t seed;
    flint_rand_init(seed);
    srand(time(NULL));
    flint_randseed(seed, rand(), rand());


    nmod_poly_mat_randtest(A, seed, 5);

    printf("A\n");
    nmod_poly_mat_print_pretty(A,"x");



    nmod_mat_t B;
    nmod_mat_init(B, rdim, cdim, 4);
    nmod_poly_mat_get_coeff_mat(B, A, 2);
    printf("\ncoefficient matrix for degree 2 of A\n");
    nmod_mat_print_pretty(B);

    slong shifts[cdim];

    shifts[0] = 1; shifts[1] = 2; shifts[2] = 3;
    printf("shifts: ");
    slongvec_print_sagemath(shifts, cdim);

    slong cols_deg[cdim];
    nmod_poly_mat_column_degree(cols_deg, A, shifts);
    slong rows_deg[rdim];
    nmod_poly_mat_row_degree(rows_deg, A, shifts);

    printf("\ncolumn degree of A with shift \n");
    slongvec_print_sagemath(cols_deg, rdim);
    printf("\nrow degree of A with shift \n");
    slongvec_print_sagemath(rows_deg, cdim);

    slong deg_A = nmod_poly_mat_degree(A);
    printf("\nA's degree: %ld\n", deg_A);


    fmpz_mat_t mat_deg;
    fmpz_mat_init(mat_deg, rdim, cdim);

    nmod_poly_mat_degree_matrix_shifted(mat_deg, A, shifts, ROW_LOWER);
    printf("\n");
    fmpz_mat_print_pretty(mat_deg);

    slong lead_pos[cdim];
    nmod_poly_mat_pivot_index(lead_pos, A, shifts, COL_UPPER);
    printf("leading position: ");
    slongvec_print_sagemath(lead_pos, cdim);

    nmod_poly_mat_leading_matrix(B, A, shifts, ROW_LOWER);
    printf("leading matrix for shift shifts of A\n");
    nmod_mat_print_sagemath(B);

    nmod_poly_mat_clear(A);
    nmod_mat_clear(B);
    flint_rand_clear(seed);
    fmpz_mat_clear(mat_deg);

    return 1;
}

int main(void)
{
    //test_nmod_mat_poly0();
    test_matpol();
    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
