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
/* clears all memory used by mmod                               */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_clear(fmpz_multimod_naive_t mmod)
{
    ulong i;

    flint_free(mmod->primes);
    flint_free(mmod->mod);
    
    for (i = 0; i < mmod->num_primes; i++)
    {
        flint_free(mmod->powers_of_two[i]);
    }
    
    flint_free(mmod->powers_of_two);
    fmpz_clear(mmod->prod);
}
