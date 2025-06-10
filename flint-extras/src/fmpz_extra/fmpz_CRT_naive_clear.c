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

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* clears all memory                                            */
/* ------------------------------------------------------------ */
void 
fmpz_CRT_naive_clear(fmpz_CRT_naive_t mCRT)
{
    ulong i;

    flint_free(mCRT->primes);
    flint_free(mCRT->inverses);
    for (i = 0; i < mCRT->num_limbs; i++)
	flint_free(mCRT->coefficients[i]);

    flint_free(mCRT->coefficients);
    flint_free(mCRT->mod);
    fmpz_clear(mCRT->prod);
}
