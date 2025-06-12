/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

#include "fmpz_extra.h"

/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_reduce(ulong * out, const fmpz_t A, const fmpz_multimod_CRT_t mmod)
{
    ulong i;
    fmpz * tmp;

    if (mmod->num_leaves == 1)
        fmpz_multimod_naive_reduce(out, A, mmod->leaves_mod[0]);
    else
    {
        tmp = (fmpz *) malloc(mmod->num_leaves * sizeof(fmpz));
        for (i = 0; i < mmod->num_leaves; i++)
            fmpz_init(tmp + i);
        
        fmpz_multi_mod_precomp(tmp, mmod->top_mod, A, 0);
        
        for (i = 0; i < mmod->num_leaves; i++)
        {
            fmpz_multimod_naive_reduce(out + mmod->size_leaves * i, tmp + i, mmod->leaves_mod[i]);
            fmpz_clear(tmp + i);
        }
        
        free(tmp);
    }
}
