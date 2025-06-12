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

#include "t-det.c"
#include "t-dixon.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_poly_mat_dixon),
    TEST_FUNCTION(nmod_poly_mat_det)
};

/* main function *************************************************************/

TEST_MAIN(tests)
