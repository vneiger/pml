
# 
# Copyright (C) 2025  PML
#
# Created from https://github.com/fplll/fplll/blob/master/configure.ac
#
#   This file is part of pml. pml is free software: you
#   can redistribute it and/or modify it under the terms of the GNU Lesser
#   General Public License as published by the Free Software Foundation,
#   either version 2.1 of the License, or (at your option) any later version.
#
#   pml is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with pml. If not, see <http://www.gnu.org/licenses/>. 


# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pml],[0.0.1])
AC_CONFIG_SRCDIR([nmod_extra/src/nmod_extra.c])
AM_MAINTAINER_MODE([enable])
AC_CANONICAL_TARGET


AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_HEADERS([config.h])

#GV do it cleaner for local include, compilation in lib (temporary) 
DEFAULT_CPPFLAGS="-O3 -Wall -I/usr/local/include -I../include" 

CPPFLAGS=${CPPFLAGS:-$DEFAULT_CPPFLAGS}

DEFAULT_LDFLAGS="-L/usr/local/lib"

LDFLAGS=${LDFLAGS:-$DEFAULT_LDFLAGS}

#Defaults 


LT_INIT


AC_PROG_INSTALL

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])     # GV outdated ? 


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strtol])


CHECK_ATOMIC()

#################################################
AS_BOX([          PML config               ])
#################################################


# Check gmp
#----------

AC_ARG_WITH(gmp,
      AS_HELP_STRING([--with-gmp=@<:@=DIR@:>@], [GMP install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
   [GMP version >= 4.2.0 needed, see http://gmplib.org])])


# Check flint
#------------

AC_ARG_WITH(flint,
      AS_HELP_STRING([--with-flint=@<:@=DIR@:>@], [Flint install directory]), [
         CPPFLAGS="$CPPFLAGS -I$withval/include"
         LDFLAGS="$LDFLAGS -L$withval/lib"
      ])

AC_CHECK_LIB(flint, __flint_get_memory_functions, , [AC_MSG_ERROR(
   [ERROR FLINT])])

LIBS="$LIBS -lm -lpthread"

################   PML FILES   #################################


AS_ECHO(["##--------------------------------------##"])

AC_CONFIG_LINKS([include/nmod_extra.h:nmod_extra/nmod_extra.h
include/fmpz_extra.h:fmpz_extra/fmpz_extra.h
include/nmod_vec_extra.h:nmod_vec_extra/nmod_vec_extra.h
include/nmod_mat_extra.h:nmod_mat_extra/nmod_mat_extra.h
])

AC_CONFIG_FILES([
Makefile
lib/Makefile
tests/Makefile
])

AC_OUTPUT

