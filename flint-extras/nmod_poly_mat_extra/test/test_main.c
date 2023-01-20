#include <flint/fmpz_mat.h>
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_extra.h"
#include "sagemath_extra.h"

#include <time.h>
#include <flint/profiler.h>

#define PRIME_30_BITS 536870923
#define PRIME_60_BITS 576460752303423619

int int64_equal(int64_t *shifts_1, int64_t *shifts_2, slong length)
{
    for (slong i = 0; i < length; i++)
        if (shifts_1[i] != shifts_2[i])
            return 0;
    return 1;
}

int test_pmbasis(void)
{
    nmod_poly_mat_t mat, first_res, second_res;
    slong rdim = 20, cdim = 10, prime = 101 , sigma = 400, len = 400,
          shifts[rdim], first_shifts[rdim], second_shifts[rdim];
    flint_rand_t state;
    timeit_t t0;

    /** init **/
    nmod_poly_mat_init(mat, rdim, cdim, prime);

    nmod_poly_mat_init(first_res, rdim, rdim, prime);
    nmod_poly_mat_init(second_res, rdim, rdim, prime);

    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    for (slong i = 0; i < rdim; i++)
        shifts[i] = rand() % 10 - 5;

    nmod_poly_mat_randtest(mat, state, len);

    timeit_start(t0);
    M_basisIV(first_res, first_shifts, mat, sigma, shifts);
    timeit_stop(t0);
    flint_printf("M_basisIV: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

    timeit_start(t0);
    PM_basis(second_res, second_shifts, mat, sigma, shifts);
    timeit_stop(t0);
    flint_printf("PM_basis: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

    if (!nmod_poly_mat_equal(first_res, second_res))
        printf("M_basis and M_basisIV don't return the same result\n");

    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(first_res);
    nmod_poly_mat_clear(second_res);

    flint_randclear(state);
    return 0;
}

/** Test
 * Verify if all versions of mbasis give the same results, test time
 */
int test_mbasis(void)
{
    nmod_poly_mat_t mat, first_res, second_res, third_res, fourth_res, fifth_res;
    slong rdim = 128, cdim = 64, prime = 3 , sigma = 32, len = 32,
          shifts[rdim], first_shifts[rdim], second_shifts[rdim],
          third_shifts[rdim], fourth_shifts[rdim], fifth_shifts[rdim];
          timeit_t t0;
          flint_rand_t state;

          /** init **/
          nmod_poly_mat_init(mat, rdim, cdim, prime);

          nmod_poly_mat_init(first_res, rdim, rdim, prime);
          nmod_poly_mat_init(second_res, rdim, rdim, prime);
          nmod_poly_mat_init(third_res, rdim, rdim, prime);
          nmod_poly_mat_init(fourth_res, rdim, rdim, prime);
          nmod_poly_mat_init(fifth_res, rdim, rdim, prime);

          flint_randinit(state);
          srand(time(NULL));
          flint_randseed(state, rand(), rand());

          nmod_poly_mat_randtest(mat, state, len);

          for (slong i = 0; i < 3; i++)
          {
              nmod_poly_mat_randtest(mat, state, len);
              //_perm_randtest(shifts, rdim, state);
              //for (slong i = 0; i < rdim; i++)
              //    shifts[i] = rand() % len - 10;

              timeit_start(t0);
              M_basis(first_res, first_shifts, mat, sigma, shifts);
              timeit_stop(t0);
              flint_printf("M_basis: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

              timeit_start(t0);
              M_basisII(second_res, second_shifts, mat, sigma, shifts);
              timeit_stop(t0);
              flint_printf("M_basisII: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

              timeit_start(t0);
              M_basisIII(third_res, third_shifts, mat, sigma, shifts);
              timeit_stop(t0);
              flint_printf("M_basisIII: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

              timeit_start(t0);
              M_basisIV(fourth_res, fourth_shifts, mat, sigma, shifts);
              timeit_stop(t0);
              flint_printf("M_basisIV: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

              timeit_start(t0);
              M_basisV(fifth_res, fifth_shifts, mat, sigma, shifts);
              timeit_stop(t0);
              flint_printf("M_basisV: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

              /** test **/
              if (!nmod_poly_mat_equal(first_res, second_res))
              {
                  printf("M_basis and M_basisII don't return the same result\n");
                  nmod_poly_mat_print(first_res,"X");
                  nmod_poly_mat_print(second_res,"X");
              }

              if (!nmod_poly_mat_equal(third_res, first_res))
              {
                  printf("M_basis and M_basisIII don't return the same result\n");
                  nmod_poly_mat_print(first_res,"X");
                  nmod_poly_mat_print(third_res,"X");
              }

              if (!nmod_poly_mat_equal(fourth_res, first_res))
              {
                  printf("M_basis and M_basisIV don't return the same result\n");
                  nmod_poly_mat_print(first_res,"X");
                  nmod_poly_mat_print(fourth_res,"X");
              }

              if (!nmod_poly_mat_equal(fifth_res, first_res))
              {
                  printf("M_basis and M_basisV don't return the same result\n");
                  nmod_poly_mat_print(first_res,"X");
                  nmod_poly_mat_print(fifth_res,"X");
              }

              if (!int64_equal(first_shifts, second_shifts, rdim))
                  printf("M_basis and M_basisII don't return the same shifts result\n");

              if (!int64_equal(first_shifts, third_shifts, rdim))
                  printf("M_basis and M_basisIII don't return the same shifts result\n");

              if (!int64_equal(first_shifts, fourth_shifts, rdim))
                  printf("M_basis and M_basisIV don't return the same shifts result\n");

              if (!int64_equal(first_shifts, fifth_shifts, rdim))
                  printf("M_basis and M_basisV don't return the same shifts result\n");

              printf("\n");
          }
          /** clear **/
          nmod_poly_mat_clear(mat);
          nmod_poly_mat_clear(first_res);
          nmod_poly_mat_clear(second_res);
          nmod_poly_mat_clear(third_res);
          nmod_poly_mat_clear(fourth_res);
          nmod_poly_mat_clear(fifth_res);

          flint_randclear(state);
          return 0;
}

/** Test
 * Verify if nmod_poly_mat_mul + coefficient_matrix and
 * nmod_mat_poly_naive_mul_coef
 * return the same result.
 * Verify if nmod_poly_mat -> nmod_mat_poly -> nmod_poly_mat
 * unchanged the initial matrix
 */
int test_nmod_mat_poly(void)
{
    nmod_mat_t prod, prod_2;
    nmod_poly_mat_t mat, mat_2, res, res_2, poly_prod;
    nmod_mat_poly_t mat_repr, mat_repr_2;
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
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    /** test poly_mat <-> mat_poly **/
    nmod_poly_mat_randtest(mat, state, len);
    nmod_poly_mat_randtest(mat_2, state, len_2);

    nmod_mat_poly_init_setII(mat_repr, mat, 10);
    nmod_mat_poly_init_setII(mat_repr_2, mat_2, 9);

    nmod_poly_mat_set_from_mat_poly(res, mat_repr);
    nmod_poly_mat_set_from_mat_poly(res_2, mat_repr_2);

    if (!nmod_poly_mat_equal(mat, res))
        printf("error mat not equals");
    if (!nmod_poly_mat_equal(mat_2, res_2))
        printf("\nerror mat_2 not equals\n");

    /** test mul: poly_mat <-> mat_poly **/
    nmod_poly_mat_mul(poly_prod, res, res_2);
    nmod_poly_mat_coefficient_matrix(prod, poly_prod, k);

    nmod_poly_mat_print(poly_prod, x);

    nmod_mat_poly_naive_mul_coef(prod_2, mat_repr, mat_repr_2, k);


    //print
    nmod_mat_print_pretty(prod);
    nmod_mat_print_pretty(prod_2);

    if (!nmod_mat_equal(prod_2, prod))
        printf("error mul not equals");

    /** clear **/
    nmod_mat_poly_clear(mat_repr);
    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(res);

    nmod_mat_poly_clear(mat_repr_2);
    nmod_poly_mat_clear(mat_2);
    nmod_poly_mat_clear(res_2);

    nmod_poly_mat_clear(poly_prod);
    nmod_mat_clear(prod);
    nmod_mat_clear(prod_2);

    flint_randclear(state);
    free(x);

    return 0;
}

/** Test
 * Verify if structured_multiplication blocks and
 * list_structured_multiplication_blocks and basis + nmod_poly_mat_mul
 * return the same result
 */
//int test_structured_mul_blocks(void)
//{
//    nmod_poly_mat_t mat, res, res_2, resid, resid_2, resid_3;
//    nmod_mat_poly_t mat_repr;
//    nmod_mat_t A, constant_mat;
//    slong rank, rdim = 8, cdim = 4, prime = 101, len = 10, sigma = 30,
//          shifts[rdim], res_shifts[rdim],
//          *perm = _perm_init(rdim);
//
//          /** init **/
//          nmod_poly_mat_init(mat, rdim, cdim, prime);
//
//          nmod_poly_mat_init(resid, rdim, cdim, prime);
//          nmod_poly_mat_init(resid_2, rdim, cdim, prime);
//          nmod_poly_mat_init(resid_3, rdim, cdim, prime);
//
//          nmod_poly_mat_init(res, rdim, rdim, prime);
//          nmod_poly_mat_init(res_2, rdim, rdim, prime);
//
//          flint_rand_t state;
//          flint_randinit(state);
//          srand(time(NULL));
//          flint_randseed(state, rand(), rand());
//
//          nmod_poly_mat_randtest(mat, state, len);
//          nmod_poly_mat_set(resid, mat);
//          nmod_poly_mat_set(resid_2, mat);
//          nmod_mat_poly_init_setII(mat_repr, mat, sigma);
//
//
//          for (slong i = 0; i < rdim; i++)
//              shifts[i] =  rand() % 10 - 5;
//
//          nmod_mat_init(constant_mat, rdim, cdim, prime);
//
//          nmod_poly_mat_coefficient_matrix(constant_mat, mat, 0);
//          rank = mbasis1_for_mbasis(A, res_shifts, perm, constant_mat, shifts);
//          mbasis1(res_2, res_shifts, constant_mat, shifts);
//          nmod_poly_mat_mul(resid_2, res_2, resid_2);
//
//          nmod_poly_mat_one(res);
//          structured_multiplication_blocks(res, A, perm, rank);
//
//          structured_multiplication_blocks(resid, A, perm, rank);
//
//          structured_list_multiplication_blocks(mat_repr, A, perm, rank, 0, sigma);
//
//          nmod_mat_poly_to_poly_mat(resid_3, mat_repr);
//
//          // print
//          char *x = malloc(1);
//          nmod_poly_mat_print(res, x);
//          nmod_poly_mat_print(res_2, x);
//          nmod_poly_mat_print(resid, x);
//          nmod_poly_mat_print(resid_2, x);
//          nmod_poly_mat_print(resid_3, x);
//          free(x);
//
//          // test
//          if (!nmod_poly_mat_equal(res, res_2))
//              printf("mbasis1 and mbasis1_for_mbasis not same result");
//
//          if (!nmod_poly_mat_equal(resid, resid_2))
//              printf("mbasis1 and mbasis1_for_mbasis not same result");
//
//          if (!nmod_poly_mat_equal(resid, resid_3))
//              printf("mbasis1 and list_struct_mul not same result");
//
//          /** clear **/
//          nmod_mat_clear(A);
//          nmod_mat_clear(constant_mat);
//
//          nmod_poly_mat_clear(mat);
//          nmod_poly_mat_clear(res);
//          nmod_poly_mat_clear(resid);
//          nmod_poly_mat_clear(res_2);
//          nmod_poly_mat_clear(resid_2);
//          nmod_poly_mat_clear(resid_3);
//
//          nmod_mat_poly_clear(mat_repr);
//
//          _perm_clear(perm);
//          flint_randclear(state);
//          return 1;
//}


/** Test
 * verify if there are leaks in basis.c functions
 **/
int test_basis(void)
{
    nmod_mat_t mat;
    slong rdim = 9, cdim = 3, prime = 7;

    nmod_mat_init(mat, rdim, cdim, prime);


    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    nmod_mat_randtest(mat, state);


    nmod_poly_mat_t res;
    nmod_poly_mat_init(res, rdim, rdim, prime);
    nmod_poly_mat_zero(res);

    int64_t res_shift[rdim];

    int64_t shift[rdim];

    for (slong i = 0; i < rdim; i++)
        shift[i] =  rand() % 10 - 5;

    mbasis1(res, res_shift, mat, shift);


    printf("Matrix\n");
    nmod_mat_print_pretty(mat);
    printf("Result mbasis1\n");
    nmod_poly_mat_print_pretty(res,"x");
    nmod_poly_mat_clear(res);
    flint_randclear(state);

    nmod_mat_t res2;
    int64_t res_shift2[rdim];
    slong res_perm2[rdim];
    slong rank = mbasis1_for_mbasis(res2, res_shift2, res_perm2, mat, shift);

    nmod_mat_clear(mat);

    printf("\nMatrix A from basis = [[x,0],[A,1]]\n");
    nmod_mat_print_pretty(res2);

    printf("\nshifts = ");
    for(slong i = 0; i < rdim; i++)
        printf("%ld ", res_shift2[i]);

    printf("\npermutation = ");
    for(slong i = 0; i < rdim; i++)
        printf("%ld ", res_perm2[i]);


    printf("\nrank = %ld\n", rank);

    nmod_mat_clear(res2);
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
    flint_randinit(seed);
    srand(time(NULL));
    flint_randseed(seed, rand(), rand());


    nmod_poly_mat_randtest(A, seed, 5);

    printf("A\n");
    nmod_poly_mat_print_pretty(A,"x");



    nmod_mat_t B;
    nmod_mat_init(B, rdim, cdim, 4);
    nmod_poly_mat_coefficient_matrix(B, A, 2);
    printf("\ncoefficient matrix for degree 2 of A\n");
    nmod_mat_print_pretty(B);

    int64_t shifts[cdim];

    shifts[0] = 1; shifts[1] = 2; shifts[2] = 3;
    printf("shifts: ");
    slongvec_print_sagemath(shifts, cdim);

    int64_t cols_deg[cdim];
    nmod_poly_mat_column_degree_shifted(cols_deg, A, shifts);
    int64_t rows_deg[rdim];
    nmod_poly_mat_row_degree_shifted(rows_deg, A, shifts);

    printf("\ncolumn degree of A with shift \n");
    slongvec_print_sagemath(cols_deg, rdim);
    printf("\nrow degree of A with shift \n");
    slongvec_print_sagemath(rows_deg, cdim);

    slong deg_A = nmod_poly_mat_degree(A);
    printf("\nA's degree: %ld\n", deg_A);


    orientation_t row_wise = 0;
    fmpz_mat_t mat_deg;
    fmpz_mat_init(mat_deg, rdim, cdim);

    nmod_poly_mat_degree_matrix_shifted(mat_deg, A, shifts, row_wise);
    printf("\n");
    fmpz_mat_print_pretty(mat_deg);

    int64_t lead_pos[cdim];
    nmod_poly_mat_pivot_index_shifted_columnwise(lead_pos, A, shifts);
    printf("leading position: ");
    slongvec_print_sagemath(lead_pos, cdim);

    nmod_poly_mat_leading_matrix_shifted(B, A, shifts, row_wise);
    printf("leading matrix for shift shifts of A\n");
    nmod_mat_print_sagemath(B);

    nmod_poly_mat_clear(A);
    nmod_mat_clear(B);
    flint_randclear(seed);
    fmpz_mat_clear(mat_deg);

    return 1;
}


int main(void)
{
    //test_nmod_mat_poly();
    test_mbasis();
    //test_matpol();
    //test_basis();
    //test_structured_mul_blocks();
    test_pmbasis();
    return EXIT_SUCCESS;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
