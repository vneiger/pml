/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* clear
/*------------------------------------------------------------*/
void nmod_multimod_CRT_clear(nmod_multimod_CRT_t C)
{
#if FLINT_HAVE_FFT_SMALL
    if (C->p >= (1L << 50)) // large modulus
#endif  // FLINT_HAVE_FFT_SMALL
        if (C->num_primes > 1)
            _nmod_vec_clear(C->data);
}
