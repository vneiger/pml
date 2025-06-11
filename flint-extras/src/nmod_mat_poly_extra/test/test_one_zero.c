/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <time.h>
#include <stdlib.h>

#include <flint/flint.h>

#include "nmod_mat_poly.h"

int main(int argc, char *argv[])
{
    flint_rand_t state;
    flint_rand_init(state);
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), rand());

    printf("building random matrix over GF(7)...\n");
    nmod_mat_poly_t matp;
    nmod_mat_poly_init(matp, 4, 5, 7);
    nmod_mat_poly_rand(matp, state, 3);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if zero: %d\n", nmod_mat_poly_is_zero(matp));

    printf("\nnow making it zero\n");
    nmod_mat_poly_zero(matp);
    nmod_mat_poly_print_pretty(matp);

    printf("testing if zero: %d\n", nmod_mat_poly_is_zero(matp));
    printf("testing if one: %d\n", nmod_mat_poly_is_one(matp));

    printf("\nnow making it one\n");
    nmod_mat_poly_one(matp);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if one: %d\n", nmod_mat_poly_is_one(matp));

    printf("\nmaking it random of larger degree again\n");
    nmod_mat_poly_rand(matp, state, 3);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if one: %d\n", nmod_mat_poly_is_zero(matp));

    printf("\nnow making it one\n");
    nmod_mat_poly_one(matp);
    nmod_mat_poly_print_pretty(matp);
    printf("testing if one: %d\n", nmod_mat_poly_is_one(matp));

    // clears memory
    nmod_mat_poly_clear(matp);
    flint_rand_clear(state);

	return 0;
}
