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

#include "t-add_shift_left.c"
#include "t-det.c"
#include "t-hermite_diagonal.c"
#include "t-hermite_normal_form.c"
#include "t-inv_series.c"
#include "t-kernel.c"
#include "t-linsolve_geometric.c"
#include "t-linsolve_via_kernel.c"
#include "t-linsolve_via_series.c"
#include "t-degree_det_bound.c"
#include "t-degree_sum.c"
#include "t-mul_sd_fft_direct.c"
#include "t-mul_sd_fft_matmul.c"
#include "t-mul_geometric.c"
#include "t-mul_vandermonde.c"
#include "t-mul_waksman.c"
#include "t-mulmid.c"
#include "t-mulmid_precomp_right.c"
#include "t-rank_at_point.c"
#include "t-set_from_mat_poly.c"
#include "t-shift_right_trunc.c"
#include "t-solve_series_dixon.c"
#include "t-solve_series_high_order_lifting.c"
#include "t-mbasis.c"
#include "t-pmbasis.c"
#include "t-pmintbasis.c"
#include "t-pmintbasis_nondistinct.c"
#include "t-pmintbasis_geometric.c"
#include "t-pmintbasis_geometric_auto.c"
#include "t-reverse.c"
#include "t-rand.c"
#include "t-row_basis.c"
#include "t-row_basis_left_factor.c"
#include "t-weak_popov_form.c"
#include "t-krylov_kernel_basis.c"
#include "t-relation_basis_krylov.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(nmod_poly_mat_add_shift_left),
    TEST_FUNCTION(nmod_poly_mat_det),
    TEST_FUNCTION(nmod_poly_mat_hnf),
    TEST_FUNCTION(nmod_poly_mat_hnf_diagonal),
    TEST_FUNCTION(nmod_poly_mat_inv_series),
    TEST_FUNCTION(nmod_poly_mat_kernel),
    TEST_FUNCTION(nmod_poly_mat_linsolve_geometric),
    TEST_FUNCTION(nmod_poly_mat_linsolve_via_kernel),
    TEST_FUNCTION(nmod_poly_mat_linsolve_via_series),
    TEST_FUNCTION(nmod_poly_mat_degree_det_bound),
    TEST_FUNCTION(nmod_poly_mat_degree_sum),
    TEST_FUNCTION(nmod_poly_mat_mbasis),
    TEST_FUNCTION(nmod_poly_mat_mul_sd_fft_direct),
    TEST_FUNCTION(nmod_poly_mat_mul_sd_fft_matmul),
    TEST_FUNCTION(nmod_poly_mat_mul_geometric),
    TEST_FUNCTION(nmod_poly_mat_mul_waksman),
    TEST_FUNCTION(nmod_poly_mat_mul_vandermonde),
    TEST_FUNCTION(nmod_poly_mat_mulmid),
    TEST_FUNCTION(nmod_poly_mat_mulmid_precomp_right),
    TEST_FUNCTION(nmod_poly_mat_rank_at_point),
    TEST_FUNCTION(nmod_poly_mat_set_from_mat_poly),
    TEST_FUNCTION(nmod_poly_mat_shift_right_trunc),
    TEST_FUNCTION(nmod_poly_mat_solve_series_dixon),
    TEST_FUNCTION(nmod_poly_mat_solve_series_high_order_lifting),
    TEST_FUNCTION(nmod_poly_mat_pmbasis),
    TEST_FUNCTION(nmod_poly_mat_pmintbasis),
    TEST_FUNCTION(nmod_poly_mat_pmintbasis_nondistinct),
    TEST_FUNCTION(nmod_poly_mat_pmintbasis_geometric),
    TEST_FUNCTION(nmod_poly_mat_pmintbasis_geometric_auto),
    TEST_FUNCTION(nmod_poly_mat_reverse),
    TEST_FUNCTION(nmod_poly_mat_rand),
    TEST_FUNCTION(nmod_poly_mat_row_basis),
    TEST_FUNCTION(nmod_poly_mat_row_basis_left_factor),
    TEST_FUNCTION(nmod_poly_mat_weak_popov_form),
    TEST_FUNCTION(nmod_poly_mat_krylov_kernel_basis),
    TEST_FUNCTION(nmod_poly_mat_relation_basis_krylov),
};

/* main function *************************************************************/

TEST_MAIN(tests)
