/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"
#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* compute A s.t. A mod primes[j] = m[j] for all j              */
/* ------------------------------------------------------------ */
void 
fmpz_CRT_naive_CRT(fmpz_t A, nn_srcptr m, const fmpz_CRT_naive_t mCRT)
{
    fmpz_t comb;
    nn_ptr m_premul;
    ulong i;

    fmpz_init(comb);
    m_premul = _nmod_vec_init(mCRT->num_primes);
    for (i = 0; i < mCRT->num_primes; i++)
	m_premul[i] = nmod_mul(m[i], mCRT->inverses[i], mCRT->mod[i]);

    fmpz_CRT_naive_combine(comb, m_premul, mCRT);
    fmpz_mod(A, comb, mCRT->prod);

    fmpz_clear(comb);
    _nmod_vec_clear(m_premul);
}
