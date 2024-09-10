#include <flint/ulong_extras.h>

#include "nmod_poly_fft.h"

// this version assumes d >= 5
void _n_geometric_sequence_with_precomp(ulong * seq, ulong a, ulong d, ulong n)
{
    ulong a_pr_quo, a_pr_rem;
    n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, n);
    seq[0] = UWORD(1);
    seq[1] = n_mulmod_precomp_shoup(UWORD(1), n);
    seq[2] = a;
    seq[3] = a_pr_quo;
    // a**2 == a * a
    n_mulmod_and_precomp_shoup(seq+4, seq+5, a, a, a_pr_quo, a_pr_rem, a_pr_quo, n);
    // a**3 == a**2 * a
    n_mulmod_and_precomp_shoup(seq+6, seq+7, a, seq[4], a_pr_quo, a_pr_rem, seq[5], n);

    // a**4 = a**3 * a
    n_mulmod_and_precomp_shoup(seq+8, seq+9, a, seq[6], a_pr_quo, a_pr_rem, seq[7], n);
    // a <- a**4 and precomp for a**4
    a = seq[8];
    a_pr_quo = seq[9];
    a_pr_rem = n_mulmod_precomp_shoup_rem_from_quo(a_pr_quo, n);
    ulong i;
    for (i = 5; i+3 < d; i+=4)
    {
        n_mulmod_and_precomp_shoup(seq+(2*(i+0)), seq+(2*(i+0)+1), a, seq[2*(i-4)], a_pr_quo, a_pr_rem, seq[2*(i-4)+1], n);
        n_mulmod_and_precomp_shoup(seq+(2*(i+1)), seq+(2*(i+1)+1), a, seq[2*(i-3)], a_pr_quo, a_pr_rem, seq[2*(i-3)+1], n);
        n_mulmod_and_precomp_shoup(seq+(2*(i+2)), seq+(2*(i+2)+1), a, seq[2*(i-2)], a_pr_quo, a_pr_rem, seq[2*(i-2)+1], n);
        n_mulmod_and_precomp_shoup(seq+(2*(i+3)), seq+(2*(i+3)+1), a, seq[2*(i-1)], a_pr_quo, a_pr_rem, seq[2*(i-1)+1], n);
    }
    for ( ; i < d; i++)
    {
        n_mulmod_and_precomp_shoup(seq+(2*(i+0)), seq+(2*(i+0)+1), a, seq[2*(i-4)], a_pr_quo, a_pr_rem, seq[2*(i-4)+1], n);
    }
}

// general, any d
void n_geometric_sequence_with_precomp(ulong * seq, ulong a, ulong d, ulong n)
{
    // note: unrolling helps
    // note: computing with a**4 to reduce dependencies helps
    if (d == 1)
    {
        seq[0] = UWORD(1);
        seq[1] = n_mulmod_precomp_shoup(UWORD(1), n);
    }
    else if (d == 2)
    {
        seq[0] = UWORD(1);
        seq[1] = n_mulmod_precomp_shoup(UWORD(1), n);
        seq[2] = a;
        seq[3] = n_mulmod_precomp_shoup(a, n);
    }
    else if (d == 3)
    {
        ulong a_pr_quo, a_pr_rem;
        n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, n);
        seq[0] = UWORD(1);  // a**0
        seq[1] = n_mulmod_precomp_shoup(UWORD(1), n);
        seq[2] = a;         // a**1
        seq[3] = a_pr_quo;
        // a**2 == a * a
        n_mulmod_and_precomp_shoup(seq+4, seq+5, a, a, a_pr_quo, a_pr_rem, a_pr_quo, n);
    }
    else if (d == 4)
    {
        ulong a_pr_quo, a_pr_rem;
        n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, n);
        seq[0] = UWORD(1);  // a**0
        seq[1] = n_mulmod_precomp_shoup(UWORD(1), n);
        seq[2] = a;         // a**1
        seq[3] = a_pr_quo;
        // a**2 == a * a
        n_mulmod_and_precomp_shoup(seq+4, seq+5, a, a, a_pr_quo, a_pr_rem, a_pr_quo, n);
        // a**3 == a**2 * a
        n_mulmod_and_precomp_shoup(seq+6, seq+7, a, seq[4], a_pr_quo, a_pr_rem, seq[5], n);
    }
    else // d >= 5
        _n_geometric_sequence_with_precomp(seq, a, d, n);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
