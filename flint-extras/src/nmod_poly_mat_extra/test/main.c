/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost, Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "testing_collection.h"

#include "t-det.c"
#include "t-dixon.c"
#include "t-hermite_normal_form.c"
#include "t-middle_product_geometric.c"
#include "t-mul_geometric.c"
#include "t-mbasis.c"
#include "t-pmbasis.c"
#include "t-mul_waksman.c"
#include "t-rand.c"
#include "t-weak_popov_form.c"
#include "t-kernel.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_poly_mat_middle_product_geometric),
    TEST_FUNCTION(nmod_poly_mat_mul_geometric),
    TEST_FUNCTION(nmod_poly_mat_det),
    TEST_FUNCTION(nmod_poly_mat_dixon),
    TEST_FUNCTION(nmod_poly_mat_hnf),
    TEST_FUNCTION(nmod_poly_mat_mbasis),
    TEST_FUNCTION(nmod_poly_mat_pmbasis),
    /* TEST_FUNCTION(nmod_poly_mat_mul_waksman), */  /* TODO */
    TEST_FUNCTION(nmod_poly_mat_rand),
    TEST_FUNCTION(nmod_poly_mat_weak_popov_form),
    TEST_FUNCTION(nmod_poly_mat_kernel)
};

/* main function *************************************************************/

TEST_MAIN(tests)
