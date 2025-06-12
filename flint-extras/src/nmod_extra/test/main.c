/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML. It is adapted from the files src/flint.h.in
    in FLINT (GNU LGPL version 3 or later).

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "pml.h"
#if PML_HAVE_AVX2
# include "t-multimod_CRT_reduce.c"
#endif  /* PML_HAVE_AVX2 */
#include "t-multimod_CRT_CRT.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
#if PML_HAVE_AVX2
    TEST_FUNCTION(nmod_multimod_CRT_reduce),
#endif  /* PML_HAVE_AVX2 */
    TEST_FUNCTION(nmod_multimod_CRT_CRT),  /* TODO failing test for non-PML_HAVE_AVX2, see C file */
};

/* main function *************************************************************/

TEST_MAIN(tests)

