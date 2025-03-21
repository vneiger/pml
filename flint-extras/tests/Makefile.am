
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


AUTOMAKE_OPTIONS = gnu

EXTRA_DIST=  perfpublisher.sh 

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig

include_pmldir=$(includedir)/pml


AM_CPPFLAGS = -I../../include 

LDADD =  ../lib/libpml.la

PERFPUBLISHERFILE=tests-report.xml


VPATH= 	../fmpz_extra/tests \
		../nmod_vec_extra/tests \
		../nmod_mat_extra/tests


#       nmod_extra  TODO 
# ______________________________________





#       fmpz_extra  
# ______________________________________

TESTS = test_CRT_naive_combine \
		test_CRT_naive_CRT \
		test_multimod_CRT_CRT 


test_CRT_naive_combine_SOURCES=test_CRT_naive_combine.c 
test_CRT_naive_CRT_SOURCES=test_CRT_naive_CRT.c
test_multimod_CRT_CRT_SOURCES=test_multimod_CRT_CRT.c

#       nmod_vec_extra  
# ______________________________________

TESTS += test_nmod_vec_integer_dot_product

test_nmod_vec_integer_dot_product_SOURCES=test_nmod_vec_integer_dot_product.c 


#       nmod_mat_extra  
# ______________________________________

TESTS += test_pluq test_rand test_nullspace 

test_pluq_SOURCES=test_pluq.c 
test_rand_SOURCES=test_rand.c 
test_nullspace_SOURCES=test_nullspace.c 






check_PROGRAMS = $(TESTS)

EXTRA_PROGRAMS = 

# Perfpublisher script interaction - AC 2016/06/09
perfpublisher:
	+chmod 777 ./perfpublisher.sh
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(TESTS)" "$(CXX)"

