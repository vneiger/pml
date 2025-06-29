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
    nmod_mat_poly_init(matp, 2, 3, 7);
    nmod_mat_poly_rand(matp, state, 4);
    printf("\nPretty print:\n");
    nmod_mat_poly_print_pretty(matp);

    printf("\nGetting info: nrows: %ld; ncols: %ld; deg: %ld; len: %ld\n",
           nmod_mat_poly_nrows(matp),
           nmod_mat_poly_ncols(matp),
           nmod_mat_poly_degree(matp),
           nmod_mat_poly_length(matp));

    printf("\nGetting coefficient (k,i,j) = (1,1,1) with macro:  %ld\n",
           nmod_mat_poly_entry(matp, 1,1,1));
    printf("Getting coefficient (k,i,j) = (1,1,1) with func:  %ld\n",
           nmod_mat_poly_get_entry(matp, 1,1,1));

    nmod_mat_poly_entry(matp, 1,1,1) = 5;
    printf("Getting coefficient (1,1,1) after setting to 5 with macro:  %ld\n",
           nmod_mat_poly_entry(matp, 1,1,1));
    nmod_mat_poly_set_entry(matp,1,1,1,4);
    printf("Getting coefficient (1,1,1) after setting to 4 with func:  %ld\n",
           nmod_mat_poly_entry(matp, 1,1,1));

    printf("\nUsing nmod_mat_print in conjunction with nmod_mat_poly_lead:\n");
    nmod_mat_print(nmod_mat_poly_lead(matp));

    printf("\nTruncating at order 2:\n");
    nmod_mat_poly_truncate(matp, 2);
    nmod_mat_poly_print_pretty(matp);

    printf("\nUsing nmod_mat_print in conjunction with nmod_mat_poly_lead:\n");
    nmod_mat_print(nmod_mat_poly_lead(matp));

    // clears memory
    nmod_mat_poly_clear(matp);
    flint_rand_clear(state);

	return 0;
}
