#include <assert.h>
#include <flint/flint.h>

#include "nmod_mat_extra.h"

/** checks X is in the nullspace of A, i.e X*A == 0 */
int is_in_nullspace(nmod_mat_t X, nmod_mat_t A)
{
    nmod_mat_t Z;
    nmod_mat_init(Z, X->r, A->c, A->mod.n);
    nmod_mat_mul(Z, X, A);
    if (nmod_mat_is_zero(Z) != 1)
    {
        printf("not in nullspace\n");
        return 0;
    }
    nmod_mat_clear(Z);
    return 1;
}

/** checks X generates the nullspace of A,
 *  and has full row rank, assuming X*A == 0 */
int basis_of_nullspace(nmod_mat_t X, nmod_mat_t A)
{
    slong rank = nmod_mat_rank(A);
    slong nullity = nmod_mat_rank(X);
    if (nullity != A->r - rank)
    {
        printf("nullspace rank not equal to nullity \n");
        return 0;
    }
    if (nullity != X->r)
    {
        printf("nullspace not full row rank\n");
        return 0;
    }
    return 1;
}

/** checks X is in reduced row echelon form, defining pivot as rightmost
 * nonzero entry in a row; assumes X has full row rank */
int is_rref(nmod_mat_t X)
{
    slong * pivot = malloc(X->r * sizeof(slong));
    for (slong i = 0; i < X->r; ++i)
    {
        pivot[i] = X->c - 1;
        while (pivot[i] >= 0 && nmod_mat_entry(X, i, pivot[i]) == 0)
            --pivot[i];
        if (pivot[i] < 0)
        {
            printf("found zero row in nullspace\n");
            return 0;
        }
        if (i>0 && pivot[i] <= pivot[i-1])
        {
            printf("pivots not increasing in nullspace\n");
            return 0;
        }
        for (slong j = i+1; j < X->r; ++j)
        {
            if (nmod_mat_entry(X,j,pivot[i]))
            {
                printf("entries not zero below pivot in nullspace\n");
                return 0;
            }
        }
    }
    return 1;
}

int check(slong field_prime, slong iterations, flint_rand_t state, slong nrows, slong ncols)
{
    nmod_mat_t A;
    nmod_mat_init(A, nrows, ncols, field_prime);
    for (slong k = 0; k < iterations; ++k)
    {
        //printf("%ld -->\n",k);
        nmod_mat_randfull(A, state);
        nmod_mat_t X;
        //printf("-- %ld:kernel in progress --\n",k);
        nmod_mat_left_nullspace(X, A);
        //printf("-- %ld: kernel OK --\n",k);
        //printf("%ldOK\n",k);
        //is_in_nullspace(X, A);
        //printf("%ldOK\n",k);
        //basis_of_nullspace(X, A);
        //printf("%ldOK\n",k);
        //is_rref(X);
        //printf("%ldOK\n",k);
        if (
               is_in_nullspace(X, A) == 0 
            || basis_of_nullspace(X, A) == 0
            || is_rref(X) == 0
           )
            return 0;
        nmod_mat_clear(X);
        //printf("<-- %ld\n",k);
    }
    return 1;
}

int main(int argc, char *argv[])
{
    if (argc > 3)
    {
        printf("Usage: %s OR %s field_prime OR %s field_prime number_of_tests\n", argv[0], argv[0], argv[0]);
        return 1;
    }
    
    srand(time(NULL));
    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, rand(), rand());

    slong field_prime = (argc>=2) ? atol(argv[1]) : 3;
    slong iterations = (argc==3) ? atol(argv[2]) : 1000;

    if (check(field_prime, iterations, state, 10, 4) == 0)
        printf("BUG1\n");
    if (check(field_prime, iterations, state, 6, 6) == 0)
        printf("BUG2\n");
    if (check(field_prime, iterations, state, 4, 10) == 0)
        printf("BUG3\n");

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
