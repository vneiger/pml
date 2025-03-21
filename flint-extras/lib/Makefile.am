
# 
#   Copyright (C) 2025 PML
#
#   This file is part of pml is free software: you
#   can redistribute it and/or modify it under the terms of the GNU Lesser
#   General Public License as published by the Free Software Foundation,
#   either version 2.1 of the License, or (at your option) any later version.
#
#   mapml is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with mapml. If not, see <http://www.gnu.org/licenses/>. 


ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig

include_pmldir=$(includedir)/pml

VPATH= 	../nmod_extra/src \
		../fmpz_extra/src \
		../nmod_poly_extra/src \
		../nmod_vec_extra/src \
		../nmod_mat_extra/src \
		../sagemath_extra/src \
		../nmod_poly_mat_extra/src \
		../nmod_mat_poly_extra/src

lib_LTLIBRARIES = libpml.la

libpml_la_SOURCES=  nmod_find_root.c    \
			        nmod_multimod_CRT_clear.c   \
			        nmod_multimod_CRT_reduce.c \
			        nmod_extra.c \
			        nmod_multimod_CRT_CRT.c \
			        nmod_multimod_CRT_init.c  

libpml_la_SOURCES+= fmpz_multi_CRT_combine.c \
					fmpz_multimod_naive_clear.c \
					fmpz_CRT_naive_CRT.c	\
					fmpz_multimod_CRT_CRT.c	\
					fmpz_multimod_naive_init.c \
					fmpz_CRT_naive_clear.c	\
					fmpz_multimod_CRT_clear.c	\
					fmpz_multimod_naive_reduce.c \
					fmpz_CRT_naive_combine.c	\
					fmpz_multimod_CRT_init.c \
					fmpz_CRT_naive_init.c	\
					fmpz_multimod_CRT_reduce.c


libpml_la_SOURCES+= nmod_geometric_progression_init_set.c \
					nmod_geometric_progression_interpolate.c \
					nmod_geometric_progression_clear.c		\
					nmod_poly_rand.c \
					nmod_geometric_progression_evaluate.c


libpml_la_SOURCES+= nmod_vec_dot_product.c	\
					nmod_vec_primes.c \
					nmod_vec_dot_product_multi.c\
					nmod_vec_rand.c \
					nmod_vec_dot2_small_modulus.c \
					nmod_vec_integer_dot_product.c

#nmod_mat_mul_small_modulus.c	

libpml_la_SOURCES+= nmod_mat_mul_newdot.c	\
					nmod_mat_mul_small_modulus.c \
                 	nmod_mat_rand.c \
					nmod_mat_rotations.c \
					left_nullspace.c	\
					nmod_mat_pluq.c

libpml_la_SOURCES+= sagemath_extras.c
					

libpml_la_SOURCES+= nmod_poly_mat_middle_product_geometric.c \
					nmod_poly_mat_approximant_mbasis.c		\
					nmod_poly_mat_middle_product_naive.c \
					nmod_poly_mat_approximant_mbasis1.c	\
					nmod_poly_mat_approximant_pmbasis.c		\
					nmod_poly_mat_approximant_verification.c	\
					nmod_poly_mat_degree_matrix.c	\
					nmod_poly_mat_det.c			\
					nmod_poly_mat_hermite.c	\
					nmod_poly_mat_mul_waksman.c \
					nmod_poly_mat_io.c		\
					nmod_poly_mat_pivots.c \
					nmod_poly_mat_is_form.c		\
					nmod_poly_mat_rand.c \
					nmod_poly_mat_is_unimodular.c	\
					nmod_poly_mat_rdeg_cdeg.c \
					nmod_poly_mat_leading_matrix.c \
					nmod_poly_mat_utils.c \
					nmod_poly_mat_weak_popov.c

libpml_la_SOURCES+= nmod_mat_poly_arith.c	\
					nmod_mat_poly_mem.c		\
					nmod_mat_poly_shift.c  \
					nmod_mat_poly_io.c		\
					nmod_mat_poly_rand.c \
					nmod_mat_poly_mbasis.c		\
					nmod_mat_poly_set_from.c



