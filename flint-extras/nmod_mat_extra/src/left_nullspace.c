#include "nmod_mat_extra.h"

slong nmod_mat_left_nullspace(nmod_mat_t X, nmod_mat_t A)
{
    // helper lists of pivot|nonpivot columns of X
    slong * permutation = malloc(A->r * sizeof(slong));

    // compute compact form of left nullspace
    nmod_mat_t Y;
    const slong nullity = nmod_mat_left_nullspace_compact(Y,permutation,A);

    // deduce nullspace
    nmod_mat_init(X, nullity, A->r, A->mod.n);
    for (slong i = 0; i < nullity; ++i)
        nmod_mat_entry(X, i, permutation[i]) = 1;
    for (slong j = 0; j < A->r - nullity; ++j)
        for (slong i = 0; i < nullity; ++i)
            nmod_mat_entry(X, i, permutation[nullity+j]) = nmod_mat_entry(Y, i, j);

    free(permutation);

    return nullity;
}

slong nmod_mat_left_nullspace_compact(
                                      nmod_mat_t X,
                                      slong * permutation,
                                      nmod_mat_t A
                                     )
{
    // At <- transpose A
    nmod_mat_t At;
    nmod_mat_init(At, A->c, A->r, A->mod.n);
    nmod_mat_transpose(At,A);

    // Xt <- right kernel of At
    nmod_mat_t Xt;
    nmod_mat_init(Xt, At->c, At->c, At->mod.n);
    slong nullity = nmod_mat_nullspace(Xt, At);

    // find row rank profile of At and its complement,
    // store them concatenated in permutation (first complement, then rrp, each
    // in increasing order)
    // -> recall: row rank profile = indices of nonpivots in Xt which is in
    // reduced column echelon form
    slong rank = At->c - nullity;
    for (slong j = 0; j < nullity; ++j)
    {
        permutation[j] = Xt->r - 1;
        while (permutation[j] >= 0 && nmod_mat_entry(Xt, permutation[j], j) == 0)
            --permutation[j];
        if (permutation[j] < 0)
        {
            printf("BUG!\n");
            return -1; // should never happen!
        }
    }

    slong r = nullity;
    for (slong i = 0; i < permutation[0]; ++i)
    {
        permutation[r] = i;
        ++r;
    }
    for (slong j = 0; j < nullity-1; ++j)
    {
        for (slong i = permutation[j]+1; i < permutation[j+1]; ++i)
        {
            permutation[r] = i;
            ++r;
        }
    }
    for (slong i = permutation[nullity-1]+1; i < Xt->r; ++i)
    {
        permutation[r] = i;
        ++r;
    }

    // extract dense part of kernel
    nmod_mat_init(X, nullity, rank, A->mod.n);
    for (slong i = 0; i < nullity; ++i)
        for (slong j = 0; j < rank; ++j)
            nmod_mat_entry(X, i, j) = nmod_mat_entry(Xt, permutation[nullity+j], i);

    // clean
    nmod_mat_clear(At);
    nmod_mat_clear(Xt);

    return nullity;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
