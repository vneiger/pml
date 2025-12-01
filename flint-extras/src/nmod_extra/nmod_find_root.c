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

#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* finds an element of order at least n                       */
/* returns 0 if not found                                     */
/*------------------------------------------------------------*/
ulong nmod_find_root(ulong n, nmod_t mod)
{
    ulong q;
    for(q = 2; q < mod.n; q++)
    {
        ulong k = 1;
        ulong qk = q;
        while (qk != 1 && k < n)
        {
            qk = nmod_mul(qk, q, mod);
            k++;
        }
        if (qk != 1)
        {
            return q;
        }
    }
    return 0; 
}

