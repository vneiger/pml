/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef __NMOD_POLY_EXTRA__H
#define __NMOD_POLY_EXTRA__H

#include <flint/nmod_poly.h>

/* #include "pml.h" */
/* #include "machine_vectors.h" */

#ifdef __cplusplus
extern "C" {
#endif


/** Test whether the prime modulus `modn` allows for evaluation, interpolation,
 * and extrapolation at `len` points in geometric progression
 *
 * For this, the functions currently in FLINT require to find an element of
 * order >= 2*len, this is feasible as soon as modn > 2*len
 *
 * TODO For the moment we require modn > 10*len just to make sure we find
 * such an element easily with nmod_find_root
 * Ideally, the constraint would be lowered to 2*len, and for example one could
 * store in some nmod context a primitive root of maximal order, so that there
 * is nothing to compute when we run some geometric progression functions
 * */
#define NMOD_CAN_USE_GEOMETRIC(modn, len) ((modn) >= UWORD(10) * (len))

/** Generates random polynomial `pol` of length up to `len` with uniformly
 * random coefficients. If `len` is nonpositive, `pol` is set to zero. */
void nmod_poly_rand(nmod_poly_t pol,
                    flint_rand_t state,
                    slong len);


/** Generates random monic polynomial `pol` of length exactly `len` with
 * uniformly random coefficients. If `len` is nonpositive, `pol` is set to
 * zero. */
void nmod_poly_rand_monic(nmod_poly_t pol,
                          flint_rand_t state,
                          slong len);

#ifdef __cplusplus
}
#endif

#endif /* __NMOD_POLY_EXTRA__H */
