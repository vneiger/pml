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

#include "pml.h"
#include "machine_vectors.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#endif // __NMOD_POLY_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
