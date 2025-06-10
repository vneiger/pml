################################################################################
# preamble
################################################################################

## 
# this file is adapted from FLINT's configure.ac (June 2025)
##

define(PML_COPYRIGHT,[[
Copyright (C) 2019, 2025 Vincent Neiger, Eric Schost, Gilles Villard

This file is part of PML.

PML is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License version 2.0 (GPL-2.0-or-later) as published by the Free
Software Foundation; either version 2.0 of the License, or (at your option) any
later version. See <https://www.gnu.org/licenses/>.

]])
define(PML_COPYRIGHT_C,[[/*
    Copyright (C) 2019 Seung Gyu Hyun, Vincent Neiger, Eric Schost
    Copyright (C) 2025 Vincent Neiger, Eric Schost, Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 3 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/]])

AC_COPYRIGHT(PML_COPYRIGHT)
AH_TOP(PML_COPYRIGHT_C)

AC_PREREQ(2.62)

m4_define(PML_VERSION,m4_include([VERSION]))

AC_INIT(PML, PML_VERSION, [https://github.com/vneiger/pml/issues/], pml)
AC_CONFIG_AUX_DIR([config])
AC_LANG(C)

# NOTE: This must be after AC_INIT
PML_MAJOR_SO=1
PML_MINOR_SO=0
PML_PATCH_SO=0

# Make sure that we are in the right directory
AC_CONFIG_SRCDIR(src/nmod_extra/nmod_extra.c)

# Ensure that source directory is not already configured, so there are no
# clashes between different configurations.
AX_INIT

################################################################################
# initialize libtool
################################################################################

# Later we want to overwrite the default CFLAGS. Check now whether they are set
# or not as Libtool overwrites CFLAGS.
if test -n "${CFLAGS+x}";
then
    cflags_set="yes"
else
    cflags_set="no"
fi

if test -n "${TESTCFLAGS+x}";
then
    testcflags_set="yes"
else
    testcflags_set="no"
fi

if test -n "${CXXFLAGS+x}";
then
    cxxflags_set="yes"
else
    cxxflags_set="no"
fi

LT_INIT([disable-static])

if test "$cflags_set" = "no";
then
    CFLAGS=
fi

if test "$cxxflags_set" = "no";
then
    CXXFLAGS=
fi

################################################################################
# build system
################################################################################

if test -n "$target_alias";
then
    AC_MSG_ERROR([PML does not support the use of specifying --target.])
fi

dnl Get system triplet
dnl NOTE: This is already invoked from LT_INIT
dnl AC_CANONICAL_HOST

################################################################################
# configure headers
################################################################################

AC_CONFIG_HEADERS(src/config.h src/pml-config.h)

################################################################################
# features
################################################################################

AC_ARG_ENABLE(pthread,
[AS_HELP_STRING([--enable-pthread],[Enable pthread [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-pthread. Need yes or no.])
    ;;
esac],
enable_pthread="yes")

AC_ARG_ENABLE(reentrant,
[AS_HELP_STRING([--enable-reentrant],[Build reentrant version of library [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-reentrant. Need yes or no.])
    ;;
esac],
enable_reentrant="no")

AC_ARG_ENABLE(thread-safe,
[AS_HELP_STRING([--enable-thread-safe],[Enable thread-local storage [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-tls. Need yes or no.])
    ;;
esac],
enable_thread_safe="yes")

if test -n "$enable_tls";
then
    AC_MSG_ERROR([--enable-tls is an old feature. Please use --enable-thread-safe instead.])
fi

AC_ARG_ENABLE(assert,
[AS_HELP_STRING([--enable-assert],[Enable use of asserts [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-assert. Need yes or no.])
    ;;
esac],
enable_assert="no")

AC_ARG_ENABLE(coverage,
[AS_HELP_STRING([--enable-coverage],[Enable test coverage [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-coverage. Need yes or no.])
    ;;
esac],
enable_coverage="no")

AC_ARG_ENABLE(debug,
[AS_HELP_STRING([--enable-debug],[Compile PML with debug information [default=yes]])],
[case $enableval in
yes|no)
    if test "$enableval" = "no" && test "$enable_coverage" = "yes";
    then
        AC_MSG_ERROR([Debug information is required by coverage.])
    fi
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-debug. Need yes or no.])
    ;;
esac],
enable_debug="yes")

AC_ARG_ENABLE(dependency-tracking,
[AS_HELP_STRING([--enable-dependency-tracking],[Enable GCC automated dependency tracking [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-dependency-tracking. Need yes or no.])
    ;;
esac],
enable_dependency_tracking="yes")

AC_ARG_ENABLE(pretty-tests,
[AS_HELP_STRING([--enable-pretty-tests],[Enable pretty printing for tests [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-pretty-tests. Need yes or no.])
    ;;
esac],
enable_pretty_tests="yes")

AC_ARG_ENABLE(gmp-internals,
[AS_HELP_STRING([--enable-gmp-internals],[Enable calling GMP internals directly [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-gmp-internals. Need yes or no.])
    ;;
esac],
enable_gmp_internals="yes")

AC_ARG_ENABLE(assembly,
[AS_HELP_STRING([--enable-assembly],[Enable assembly routines (if available) [default=yes]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-assembly. Need yes or no.])
    ;;
esac],
enable_assembly="yes")

AC_ARG_ENABLE(avx2,
[AS_HELP_STRING([--enable-avx2],[Enable AVX2 instructions [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-avx2. Need yes or no.])
    ;;
esac],
enable_avx2="no")

AC_ARG_ENABLE(avx512,
[AS_HELP_STRING([--enable-avx512],[Enable AVX512 instructions [default=no]])],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-avx512. Need yes or no.])
    ;;
esac],
enable_avx512="no")

# NOTE: Maintainer level only. Avoids any sort of optimization to speed up
# builds.
AC_ARG_ENABLE(fast-build,[],
[case $enableval in
yes|no)
    ;;
*)
    AC_MSG_ERROR([Bad value $enableval for --enable-fast-build. Need yes or no.])
    ;;
esac],
enable_fast_build="no")

################################################################################
# packages
################################################################################

AC_ARG_WITH(gmp_include,
[AS_HELP_STRING([--with-gmp-include=DIR],[GMP include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
gmp_include_path="$withval"
)

AC_ARG_WITH(gmp_lib,
[AS_HELP_STRING([--with-gmp-lib=DIR],[GMP library directory])],
LDFLAGS="$LDFLAGS -L$withval"
gmp_lib_path="$withval"
)

AC_ARG_WITH(gmp,
[AS_HELP_STRING([--with-gmp[[=DIR]]],[GMP install directory])],
if test "$with_gmp" = "no";
then
    AC_MSG_FAILURE([The use of GMP is mandatory.])
elif test "$with_gmp" = "yes";
then
    :
elif test -z "$gmp_include_path" && test -z "$gmp_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    gmp_include_path="$withval/include"
    gmp_lib_path="$withval/lib"
    with_gmp="yes"
else
    AC_MSG_FAILURE([Cannot use --with-gmp along with --with-gmp-include or --with-gmp-lib.])
fi,
with_gmp="yes"
)

AC_ARG_WITH(mpfr_include,
[AS_HELP_STRING([--with-mpfr-include=DIR],[MPFR include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
mpfr_include_path="$withval"
)

AC_ARG_WITH(mpfr_lib,
[AS_HELP_STRING([--with-mpfr-lib=DIR],[MPFR library directory])],
LDFLAGS="$LDFLAGS -L$withval"
mpfr_lib_path="$withval"
)

AC_ARG_WITH(mpfr,
[AS_HELP_STRING([--with-mpfr[[=DIR]]],[MPFR install directory])],
if test "$with_mpfr" = "no";
then
    AC_MSG_FAILURE([The use of MPFR is mandatory.])
elif test "$with_mpfr" = "yes";
then
    :
elif test -z "$mpfr_include_path" && test -z "$mpfr_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    mpfr_include_path="$withval/include"
    mpfr_lib_path="$withval/lib"
    with_mpfr="yes"
else
    AC_MSG_FAILURE([Cannot use --with-mpfr along with --with-mpfr-include or --with-mpfr-lib.])
fi,
with_mpfr="yes"
)

AC_ARG_WITH(flint_include,
[AS_HELP_STRING([--with-flint-include=DIR],[FLINT include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
flint_include_path="$withval"
)

AC_ARG_WITH(flint_lib,
[AS_HELP_STRING([--with-flint-lib=DIR],[FLINT library directory])],
LDFLAGS="$LDFLAGS -L$withval"
flint_lib_path="$withval"
)

AC_ARG_WITH(flint,
[AS_HELP_STRING([--with-flint[[=DIR]]],[FLINT install directory])],
if test "$with_flint" = "no";
then
    AC_MSG_FAILURE([The use of FLINT is mandatory.])
elif test "$with_flint" = "yes";
then
    :
elif test -z "$flint_include_path" && test -z "$flint_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    flint_include_path="$withval/include"
    flint_lib_path="$withval/lib"
    with_flint="yes"
else
    AC_MSG_FAILURE([Cannot use --with-flint along with --with-flint-include or --with-flint-lib.])
fi,
with_flint="yes"
)


AC_ARG_WITH(blas_include,
[AS_HELP_STRING([--with-blas-include=DIR],[Use BLAS and specify its include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
blas_include_path="$withval"
)

AC_ARG_WITH(blas_lib,
[AS_HELP_STRING([--with-blas-lib=DIR],[Use BLAS and specify its library directory])],
LDFLAGS="$LDFLAGS -L$withval"
blas_lib_path="$withval"
)

AC_ARG_WITH(blas,
[AS_HELP_STRING([--with-blas[[=DIR]]],[Use BLAS and specify its install directory])],
if test "$with_blas" = "yes";
then
    :
elif test "$with_blas" = "no";
then
    if test -n "$blas_include_path" || test -n "$blas_lib_path";
    then
        AC_MSG_FAILURE([--with-blas-include or --with-blas-lib was specified but --with-blas=no.])
    fi
elif test -z "$blas_include_path" && test -z "$blas_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    blas_include_path="$withval/include"
    blas_lib_path="$withval/lib"
    with_blas="yes"
else
    AC_MSG_FAILURE([Cannot use --with-blas along with --with-blas-include or --with-blas-lib.])
fi,
with_blas="no"
)

if test "$with_blas" != "no" || test -n "$blas_include_path" || test -n "$blas_lib_path";
then
    with_blas="yes"
fi

AC_ARG_WITH(gc_include,
[AS_HELP_STRING([--with-gc-include=DIR],[Use GC and specify its include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
gc_include_path="$withval",
with_gc_include="no"
)

AC_ARG_WITH(gc_lib,
[AS_HELP_STRING([--with-gc-lib=DIR],[Use GC and specify its library directory])],
LDFLAGS="$LDFLAGS -L$withval"
gc_lib_path="$withval",
with_gc_lib="no"
)

AC_ARG_WITH(gc,
[AS_HELP_STRING([--with-gc[[=DIR]]],[Use GC and optionally specify its install directory])],
if test "$with_gc" = "yes";
then
    :
elif test "$with_gc" = "no";
then
    if test -n "$gc_include_path" || test -n "$gc_lib_path";
    then
        AC_MSG_FAILURE([--with-gc-include or --with-gc-lib was specified but --with-gc=no.])
    fi
elif test -z "$gc_include_path" && test -z "$gc_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    gc_include_path="$withval/include"
    gc_lib_path="$withval/lib"
    with_gc="yes"
else
    AC_MSG_FAILURE([Cannot use --with-gc along with --with-gc-include or --with-gc-lib.])
fi,
with_gc="no"
)

if test "$with_gc" != "no" || test -n "$gc_include_path" || test -n "$gc_lib_path";
then
    with_gc="yes"
fi

AC_ARG_WITH(ntl_include,
[AS_HELP_STRING([--with-ntl-include=DIR],[Use NTL (in tests) and specify its include directory])],
CPPFLAGS="$CPPFLAGS -I$withval"
ntl_include_path="$withval"
)

AC_ARG_WITH(ntl_lib,
[AS_HELP_STRING([--with-ntl-lib=DIR],[Use NTL (in tests) and specify its library directory])],
LDFLAGS="$LDFLAGS -L$withval"
ntl_lib_path="$withval"
)

AC_ARG_WITH(ntl,
[AS_HELP_STRING([--with-ntl[[=DIR]]],[Use NTL (in tests) and optionally specify its install directory])],
if test "$with_ntl" = "yes";
then
    :
elif test "$with_ntl" = "no";
then
    if test -n "$ntl_include_path" || test -n "$ntl_lib_path";
    then
        AC_MSG_FAILURE([--with-ntl-include or --with-ntl-lib was specified but --with-ntl=no.])
    fi
elif test -z "$ntl_include_path" && test -z "$ntl_lib_path";
then
    CPPFLAGS="$CPPFLAGS -I$withval/include"
    LDFLAGS="$LDFLAGS -L$withval/lib"
    ntl_include_path="$withval/include"
    ntl_lib_path="$withval/lib"
    with_ntl="yes"
else
    AC_MSG_FAILURE([Cannot use --with-ntl along with --with-ntl-include or --with-ntl-lib.])
fi,
with_ntl="no"
)

if test "$with_ntl" != "no" || test -n "$ntl_include_path" || test -n "$ntl_lib_path";
then
    with_ntl="yes"
fi


# maple interface 

AC_MSG_CHECKING([use of maple])

AC_ARG_WITH(maplec,
[AS_HELP_STRING([--with-maplec=DIR], [Maple include directory])],
AC_MSG_RESULT([yes])
CPPFLAGS="$CPPFLAGS -I$withval/include"
maple_include_path="$withval", AC_MSG_RESULT([no])
)


AC_ARG_WITH(maplelib,
[AS_HELP_STRING([--with-maplelib=DIR], [Maple library directory])], 
LDFLAGS="$LDFLAGS -L$withval"
maple_lib_path="$withval"
)

if test "$with_maplec" != "no" && test -n "$maple_include_path" && test -n "$maple_lib_path";
then
    with_maple="yes"
fi


################################################################################
# architecture specifics
################################################################################

if test "$enable_fast_build" = "no";
then
    gcc_cflags="-O3 "
fi
gcc_cflags="$gcc_cflags-pedantic -std=c11"
gcc_warnings="-Werror=implicit-function-declaration -Werror=newline-eof -Wall -Wno-stringop-overread -Wno-stringop-overflow"

dnl from FLINT:
dnl We only try to provide specifics for those systems that currently supports
dnl our assembly routines. If more combinations are wished for than what is
dnl specified, please open up an issue at
dnl <https://github.com/vneiger/pml/issues/> and we will consider it.
dnl
dnl For these systems we aim to provide:
dnl
dnl   gcc_cflags     flags for GCC-compatible compilers
dnl
dnl   asm_path       directory for assembly, relative to src/mpn_extras
dnl   param_path     directory for flint-mparam.h, relative to src/mpn_extras
dnl
dnl For x86_64 systems, we can also set:
dnl
dnl   have_avx512    system has AVX512F and AVX512DQ (we assume this implies have_avx2)
dnl   have_avx2      system has AVX2 and FMA

case $host in

  ARM64_PATTERN)
    gcc_cflags="$gcc_cflags"
    gcc_cflags_optlist="arch tune"
    ##asm_path="arm64"
    ##param_path="arm64"
    flint_cv_have_fft_small_arm_i="yes"

    # NOTE: Cortex values where taken from
    # https://developer.arm.com/Processors/Cortex-XXX
    case $host_cpu in
      armcortexa35 | armcortexa35neon)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=cortex-a35"
        ;;
      armcortexa53 | armcortexa53neon)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=cortex-a53"
        ;;
      armcortexa55 | armcortexa55neon)
        gcc_cflags_arch="-march=armv8.2-a"
        gcc_cflags_tune="-mtune=cortex-a55"
        ;;
      armcortexa57 | armcortexa57neon)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=cortex-a57"
        ;;
      armcortexa72 | armcortexa72neon)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=cortex-a72"
        ;;
      armcortexa73 | armcortexa73neon)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=cortex-a73"
        ;;
      armcortexa75 | armcortexa75neon)
        gcc_cflags_arch="-march=armv8.2-a"
        gcc_cflags_tune="-mtune=cortex-a75"
        ;;
      armcortexa76 | armcortexa76neon)
        gcc_cflags_arch="-march=armv8.2-a"
        gcc_cflags_tune="-mtune=cortex-a76"
        ;;
      armcortexa77 | armcortexa77neon)
        gcc_cflags_arch="-march=armv8.2-a"
        gcc_cflags_tune="-mtune=cortex-a77"
        ;;
      armcortexa65 | armcortexa65neon)
        gcc_cflags_arch="-march=armv8.2-a"
        gcc_cflags_tune="-mtune=cortex-a65"
        ;;
      armcortexa34 | armcortexa34neon)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=cortex-a34"
        ;;
      armcortexa78 | armcortexa78neon)
        gcc_cflags_arch="-march=armv8.2-a"
        gcc_cflags_tune="-mtune=cortex-a78"
        ;;
      armexynosm1)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=exynos-m1"
        ;;
      armthunderx)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=thunderx"
        ;;
      armxgene1)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune="-mtune=xgene1"
        ;;
      applem1*)
        gcc_cflags_arch="-march=armv8.5-a"
        gcc_cflags_tune=""
        ;;
      applem2*)
        gcc_cflags_arch="-march=armv8.6-a"
        gcc_cflags_tune=""
        ;;
      applem3*)
        gcc_cflags_arch="-march=armv8.6-a"
        gcc_cflags_tune=""
        ;;
      applem4*)
        gcc_cflags_arch="-march=armv8.7-a"
        gcc_cflags_tune=""
        ;;
      aarch64* | armv8*)
        gcc_cflags_arch="-march=armv8-a"
        gcc_cflags_tune=""
        ;;
      *)
        AC_MSG_ERROR([Error in deciding flags for $host.
Please report at <https://github.com/vneiger/pml/issues/>])
        ;;
    esac
    ;;


  X86_64_ADX_PATTERN)
    gcc_cflags="$gcc_cflags"
    gcc_cflags_optlist="arch"
    #asm_path="x86_64/broadwell"
    #param_path="x86_64/broadwell"

    have_avx512="no"
    have_avx2="yes"

    case $host_cpu in
      zen)
        gcc_cflags_arch="-march=znver1"
        #param_path="x86_64/zen3"
        ;;
      zen2)
        gcc_cflags_arch="-march=znver2"
        #param_path="x86_64/zen3"
        ;;
      zen3)
        gcc_cflags_arch="-march=znver3"
        #param_path="x86_64/zen3"
        ;;
      zen4)
        gcc_cflags_arch="-march=znver4"
        #param_path="x86_64/zen4"
        have_avx512="yes"
        ;;
      zen5)
        gcc_cflags_arch="-march=znver5 -march=znver4"
        #param_path="x86_64/zen4"  # fallback, to be refined
        have_avx512="yes"
        ;;
      broadwell)
        gcc_cflags_arch="-march=broadwell"
        #param_path="x86_64/broadwell"
        ;;
      skylake)
        gcc_cflags_arch="-march=skylake"
        #param_path="x86_64/skylake"
        ;;
      skylake_server)
        gcc_cflags_arch="-march=skylake-avx512"
        #param_path="x86_64/skylake"
        have_avx512="yes"
        ;;
      kabylake)
        gcc_cflags_arch="-march=skylake"
        #param_path="x86_64/skylake"
        ;;
      cometlake)
        gcc_cflags_arch="-march=skylake"
        #param_path="x86_64/skylake"
        ;;
      knightslanding)
        gcc_cflags_arch="-march=knl"
        #param_path="x86_64/skylake"
        ;;
      cannonlake)
        gcc_cflags_arch="-march=cannonlake"
        #param_path="x86_64/icelake"  # fallback, to be refined
        have_avx512="yes"
        ;;
      icelake)
        gcc_cflags_arch="-march=icelake-client"
        #param_path="x86_64/icelake"
        have_avx512="yes"
        ;;
      icelake_server)
        gcc_cflags_arch="-march=icelake-server"
        #param_path="x86_64/icelake"
        have_avx512="yes"
        ;;
      rocketlake)
        gcc_cflags_arch="-march=rocketlake"
        #param_path="x86_64/icelake"  # fallback, to be refined
        have_avx512="yes"
        ;;
      tigerlake)
        gcc_cflags_arch="-march=tigerlake"
        #param_path="x86_64/tigerlake"
        have_avx512="yes"
        ;;
      alderlake)
        gcc_cflags_arch="-march=alderlake"
        #param_path="x86_64/icelake"  # fallback, to be refined
        ;;
      raptorlake)
        gcc_cflags_arch="-march=raptorlake -march=alderlake"
        #param_path="x86_64/raptorlake"
        ;;
      meteorlake)
        gcc_cflags_arch="-march=meteorlake -march=alderlake"
        #param_path="x86_64/meteorlake"
        ;;
      sapphirerapids)
        gcc_cflags_arch="-march=sapphirerapids"
        #param_path="x86_64/icelake"  # fallback, to be refined
        have_avx512="yes"
        ;;
      *)
        AC_MSG_ERROR([Error in deciding flags for $host.
Please report at <https://github.com/vneiger/pml/issues/>])
        ;;
    esac

    case $host in
	  *-*-mingw* | *-*-msys | *-*-cygwin)
        flint_nonstd_abi="yes"
        #asm_path=""
        #param_path="x86_64"
	    ;;
	esac

    if test "$cflags_set" = "yes";
    then
        # We cannot be sure that AVX2 and AVX512 are enabled.
        have_avx2="no"
        have_avx512="no"
    fi

    if test "$have_avx512" = "yes";
    then
        have_avx2="yes"
    fi

    if test "$have_avx2" = "yes";
    then
        flint_cv_have_fft_small_x86_i="yes"
    fi
	;;

  X86_64_PATTERN)
    # NOTE: We do not care if the user specifies noavx!
    gcc_cflags="$gcc_cflags"
    gcc_cflags_optlist="arch"
    asm_path=""
    param_path="x86_64"

    case $host_cpu in
      coreiibr | ivybridge)
        gcc_cflags_arch="-march=ivybridge"
        ;;
      coreihwl | haswell)
        gcc_cflags_arch="-march=haswell"
        ;;
      piledriver | bd2)
        gcc_cflags_arch="-march=bdver2"
        ;;
      steamroller | bd3)
        gcc_cflags_arch="-march=bdver3"
        ;;
      excavator | bd4)
        gcc_cflags_arch="-march=bdver4"
        have_avx2="yes"
        ;;
      x86_64v3)
        gcc_cflags_arch="-march=x86-64-v3"
        have_avx2="yes"
        ;;
      x86_64v4)
        gcc_cflags_arch="-march=x86-64-v4"
        have_avx512="yes"
        ;;
    esac

    if test "$cflags_set" = "yes";
    then
        # We cannot be sure that AVX2 and AVX512 are enabled.
        have_avx2="no"
        have_avx512="no"
    fi

    if test "$have_avx512" = "yes";
    then
        have_avx2="yes"
    fi

    if test "$have_avx2" = "yes";
    then
        flint_cv_have_fft_small_x86_i="yes"
    fi
	;;

  *)
    asm_path=""
    param_path="generic"
    gcc_cflags="$gcc_cflags"
    ;;

esac

################################################################################
# programs
################################################################################

dnl NOTE: The following line is invoked from LT_INIT
dnl AC_PROG_CC

AC_PROG_MKDIR_P

FLINT_CC_IS_GCC
FLINT_CC_IS_CLANG

################################################################################
# environment variables
################################################################################

AC_ARG_VAR([TESTCFLAGS], [Choose different C compiler flags for tests [default=CFLAGS]. Can be useful to compile library at a specific optimization level while compiling tests fast.])

################################################################################
# check programs and system
################################################################################

AC_C_INLINE

if test "$compiler_c_o" = "no";
then
    AC_MSG_ERROR([Compiler must support both -c and -o simultaneously!])
fi

AC_C_BIGENDIAN(flint_big_endian=yes)

case "$host_os" in
    darwin*)
        pml_lib="libpml.dylib"
        pml_lib_major="libpml.$PML_MAJOR_SO.dylib"
        pml_lib_full="libpml.$PML_MAJOR_SO.$PML_MINOR_SO.dylib"
        extra_shared_flags="-install_name '`pwd`/$pml_lib_full' -compatibility_version $PML_MAJOR_SO.$PML_MINOR_SO -current_version $PML_MAJOR_SO.$PML_MINOR_SO.$PML_PATCH_SO"
        pml_dylib="yes"
        ;;
    cygwin|mingw*|msys)
        pml_lib="libpml.dll"
        pml_lib_major="libpml.$PML_MAJOR_SO.dll"
        pml_lib_full="libpml.$PML_MAJOR_SO.$PML_MINOR_SO.$PML_PATCH_SO.dll"
        pml_implib="$pml_lib.$libext"
        extra_shared_flags="-static-libgcc $wl--export-all-symbols,--out-implib,$pml_implib"
        pml_dlllib="yes"
        ;;
    *)
        pml_lib="libpml.so"
        pml_lib_major="libpml.so.$PML_MAJOR_SO"
        pml_lib_full="libpml.so.$PML_MAJOR_SO.$PML_MINOR_SO.$PML_PATCH_SO"
        extra_shared_flags="$wl-soname,$pml_lib_major"
        pml_solib="yes"
        ;;
esac

pml_lib_static="libpml.$libext"

if test "`uname -s`" = "android";
then
    extra_shared_flags="$wl-soname,$pml_lib"
fi

AC_SUBST(WL,$wl)

AC_SUBST(PML_LIB,$pml_lib)
AC_SUBST(PML_LIB_MAJOR,$pml_lib_major)
AC_SUBST(PML_LIB_FULL,$pml_lib_full)
AC_SUBST(PML_IMPLIB,$pml_implib)
AC_SUBST(PML_LIB_STATIC,$pml_lib_static)
AC_SUBST(EXTRA_SHARED_FLAGS,$extra_shared_flags)

if test "$pml_dlllib" = "yes";
then
    AC_SUBST(PML_DLLLIB,1)
else
    AC_SUBST(PML_DLLLIB,0)
fi

if test "$pml_dylib" = "yes";
then
    AC_SUBST(PML_DYLIB,1)
else
    AC_SUBST(PML_DYLIB,0)
fi

if test "$pml_solib" = "yes";
then
    AC_SUBST(PML_SOLIB,1)
else
    AC_SUBST(PML_SOLIB,0)
fi

case "$host_os" in
    cygwin|mingw*|msys)
        if test "$DLLTOOL" = "false";
        then
            AC_MSG_ERROR(["Couldn't find any dlltool that is required to create import libraries for Windows-type systems."])
        fi
        ;;
esac

################################################################################
# check GMP
################################################################################

FLINT_CHECK_GMP_H(6,2,1)

################################################################################
# check MPFR
################################################################################

FLINT_CHECK_MPFR_H(4,1,0)

################################################################################
# check FLINT
################################################################################

PML_CHECK_FLINT_H(3,0,1)

################################################################################
# check headers
################################################################################

# Mandatory headers
# The following headers are checked previously:
#   stdio.h stdlib.h stdint.h string.h
AC_CHECK_HEADERS([stdarg.h math.h float.h errno.h],,
                 AC_MSG_ERROR([Could not find a mandatory header!]))

# Optional headers
# The following headers are checked previously:
#   unistd.h
AC_CHECK_HEADERS([fenv.h alloca.h malloc.h sys/param.h windows.h pthread_np.h])

if test "$enable_pthread" = "yes";
then
    AC_CHECK_HEADER([pthread.h],,
                    AC_MSG_ERROR([Could not find pthread.h!]))
fi

if test "$enable_assert" = "yes";
then
    AC_CHECK_HEADER([assert.h],,
                    AC_MSG_ERROR([Could not find assert.h!]))
fi

if test "$with_blas" = "yes";
then
    AC_CHECK_HEADER([cblas.h],,
                    AC_MSG_ERROR([Could not find cblas.h!]))
fi

if test "$with_gc" = "yes";
then
    AC_CHECK_HEADER([gc.h],,
                    AC_MSG_ERROR([Could not find gc.h!]))
fi

if test "$with_ntl" = "yes";
then
    AC_PROG_CXX
    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([NTL/ZZ.h],,
                    AC_MSG_ERROR([Could not find NTL/ZZ.h!]))
    AC_LANG_POP([C++])
fi

if test "$with_maple" = "yes";
then
    AC_CHECK_HEADER([maplec.h],,
                    AC_MSG_ERROR([Could not find maplec.h!]))
fi


################################################################################
# check libraries
################################################################################

# With certain CFLAGS, the compiler may return an error, thinking that the
# compiler has not found the function in the library.
save_CFLAGS="$CFLAGS"
CFLAGS=""

AC_SEARCH_LIBS([atan2],[m],
    [],
    [AC_MSG_ERROR(["The C math library was not found!"])])

AC_SEARCH_LIBS([__gmpz_init],[gmp],
    [],
    [AC_MSG_ERROR(["GMP library was not found.  If you indeed have GMP installed, please
submit a bug report to <https://github.com/vneiger/pml/issues/> so
that we can either fix the issue or give a more proper error message."])])

AC_SEARCH_LIBS([__gmpn_mul_basecase],[gmp],
    [],
    [AC_MSG_ERROR([`mpn_mul_basecase' was not found in the GMP library, which is required by FLINT.
Please submit a bug report to <https://github.com/vneiger/pml/issues/>.])])

if test "$enable_gmp_internals" = "yes";
then
    dnl The following functions should exist in the GMP library regardless.
    AC_SEARCH_LIBS([__gmpn_gcd_11],[gmp],
        [],
        [AC_MSG_ERROR(["`mpn_gcd_11' was not found in the GMP library. It is needed to enable GMP
internals."])])

    AC_SEARCH_LIBS([__gmpn_div_q],[gmp],
        [],
        [AC_MSG_ERROR([`mpn_div_q' was not found in the GMP library. It is needed to enable GMP
internals.])])

fi

if test "$enable_mpfr_check" = "yes";
then
    AC_SEARCH_LIBS([mpfr_init],[mpfr],
        [],
        [AC_MSG_ERROR([MPFR library was not found.  If you indeed have MPFR installed, please
submit a bug report to <https://github.com/vneiger/pml/issues/> so
that we can either fix the issue or give a more proper error message.])])

    AC_SEARCH_LIBS([mpfr_round_p],[mpfr],
        [],
        [AC_MSG_ERROR(["`mpfr_round_p' was not found in the MPFR library. It is required for PML/FLINT."])])

    AC_SEARCH_LIBS([mpfr_mulhigh_n],[mpfr],
        [],
        [AC_MSG_ERROR(["`mpfr_mulhigh_n' was not found in the MPFR library. It is required for PML/FLINT."])])

    AC_SEARCH_LIBS([mpfr_sqrhigh_n],[mpfr],
        [],
        [AC_MSG_ERROR(["`mpfr_sqrhigh_n' was not found in the MPFR library. It is required for PML/FLINT."])])
else
    LIBS="$LIBS -lmpfr"
fi

if test "$enable_flint_check" = "yes";
then
    AC_SEARCH_LIBS([flint_malloc],[flint],
        [],
        [AC_MSG_ERROR([FLINT library was not found.  If you indeed have FLINT installed, please
submit a bug report to <https://github.com/vneiger/pml/issues/> so
that we can either fix the issue or give a more proper error message.])])

else
    LIBS="$LIBS -lflint"
fi


if test "$enable_pthread" = "yes";
then
    # check if compiler accepts the preferred -pthread
    AX_CHECK_COMPILE_FLAG([-pthread],
        [LIBS="-pthread $LIBS"],
        # check if the compiler can do -lpthread
        AC_SEARCH_LIBS([pthread_create],[pthread],
            [],
            AC_MSG_ERROR([Couldn't link against POSIX threads via -pthread or -lpthread while
--enable-pthread was specified.])
        )
    )
fi

if test "$with_blas" = "yes";
then
    AC_SEARCH_LIBS([cblas_dgemm],[cblas openblas blas],[],
        [AC_MSG_ERROR([BLAS library was not found.  If you indeed have BLAS installed, please
submit a bug report to <https://github.com/vneiger/pml/issues/> so
that we can either fix the issue or give a more proper error message.])]
    )
fi

if test "$with_ntl" = "yes";
then
    FLINT_CHECK_NTL(,[AC_MSG_ERROR([NTL library was not found.  If you indeed have NTL installed, please
submit a bug report to <https://github.com/vneiger/pml/issues/> so
that we can either fix the issue or give a more proper error message.])])
fi

if test "$with_maple" = "yes";
then 
    AC_SEARCH_LIBS([MapleToInteger64],[maplec],,
    [AC_MSG_ERROR(["MAPLE library was not found.  If you indeed have maple installed, please
submit a bug report to <https://github.com/vneiger/pml/issues/> so
that we can either fix the issue or give a more proper error message.])])
fi

if test "$with_maple" = "yes";
then AC_SUBST(WANT_MAPLE,1)
else AC_SUBST(WANT_MAPLE,0)
fi

################################################################################
# check settings and environment
################################################################################

AC_SUBST(PIC_FLAG,$pic_flag)

FLINT_CHECK_PRAGMAS

# ensure alloca works
AC_MSG_CHECKING([if alloca works])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#  ifdef _MSC_VER
#   include <malloc.h>
#   define alloca _alloca
#  else
#   include <alloca.h> /* We assume then that you have alloca.h */
#  endif
# endif
#endif],
[char * p = alloca(128);])],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])
AC_MSG_ERROR([Couldn't find alloca, which is required for FLINT. Please submit a bug
report to <https://github.com/vneiger/pml/issues/> and specify your
operating system.])])

AC_CHECK_FUNCS([aligned_alloc _aligned_malloc])

################################################################################
# CFLAGS
################################################################################

if test "$flint_cv_cc_is_clang" = "yes";
then
    # The following is needed for Clang to check for unknown options.
    AX_CHECK_COMPILE_FLAG([-Werror=unknown-warning-option],
                          [CFLAGS="-Werror=unknown-warning-option"])
fi

if test "$enable_coverage" = "yes";
then
    AX_CHECK_COMPILE_FLAG([--coverage],
        [save_CFLAGS="--coverage $save_CFLAGS"],
        [AC_MSG_ERROR([$CC does not to support test coverage flags])])
    LDFLAGS="--coverage $LDFLAGS"

    if test "$enable_pthread" = "yes";
    then
        AX_CHECK_COMPILE_FLAG([-fprofile-update=atomic],
            [save_CFLAGS="-fprofile-update=atomic $save_CFLAGS"],
            [AC_MSG_ERROR([$CC does not to support multi-threaded test coverage flags])])
    fi
fi

if test "$enable_debug" = "yes";
then
    if test "$ac_cv_prog_cc_g" = "yes";
    then
        save_CFLAGS="-g $save_CFLAGS"
    else
        AC_MSG_ERROR([$CC does not support the flag -g needed to generate debug information])
    fi
fi

if test -n "$ABI";
then
    AX_CHECK_COMPILE_FLAG([-m$ABI],
        [save_CFLAGS="-m$ABI $save_CFLAGS"],
        [AC_MSG_ERROR([$CC does not support the ABI flag -m$ABI])]
    )
fi

if test "$enable_avx2" = "yes";
then
    AX_CHECK_COMPILE_FLAG([-mavx2],
    [save_CFLAGS="-mavx2 $save_CFLAGS"],
    [AC_MSG_ERROR([$CC does not support the flag -mavx2 needed for AVX2 instructions])])
    AX_CHECK_COMPILE_FLAG([-mfma], [save_CFLAGS="-mfma $save_CFLAGS"])
fi

if test "$enable_avx512" = "yes";
then
    AX_CHECK_COMPILE_FLAG([-mavx512f],
    [save_CFLAGS="-mavx512f $save_CFLAGS"],
    [AC_MSG_ERROR([$CC does not support the flag -mavx512f needed for AVX512 instructions])])
fi

if test "$cflags_set" = "no";
then
    for flag in $gcc_cflags $gcc_warnings;
    do
        AX_CHECK_COMPILE_FLAG([$flag],[save_CFLAGS="$flag $save_CFLAGS"])
    done
    for opt in $gcc_cflags_optlist;
    do
        eval tmp=\"\$gcc_cflags_${opt}\" # " Just to fix highlighting
        # found_single is for checking if we have found a flag that should be
        # alone in its type. Example: On Skylake, -march=skylake should not be
        # accompanied by -march=broadwell; instead, if we find that
        # -march=skylake is accepted, we do not want to check -march=broadwell.
        found_single="no"
        for flag in $tmp;
        do
            if test "$found_single" = "yes";
            then
                found_single="no"
                break
            elif test "$opt" = "arch";
            then
                AX_CHECK_COMPILE_FLAG([$flag],
                                      [save_CFLAGS="$flag $save_CFLAGS"
                                       found_single="yes"],
                                      [unset flint_cv_have_fft_small_x86_i
                                       unset flint_cv_have_fft_small_arm_i])
            else
                AX_CHECK_COMPILE_FLAG([$flag],[save_CFLAGS="$flag $save_CFLAGS"])
            fi
        done
    done
    if test "$enable_fast_build" = no;
    then
        AX_CHECK_COMPILE_FLAG([-funroll-loops],[unroll_loops="yes"])
    fi
    if test "$unroll_loops" = "yes";
    then
        AC_SUBST(UNROLL_LOOPS,1)
        AC_DEFINE(FLINT_UNROLL_LOOPS, 1, [Define to locally unroll some loops])
    else
        AC_SUBST(UNROLL_LOOPS,0)
    fi
    CFLAGS="$save_CFLAGS"
else
    CFLAGS="$save_CFLAGS"
    AX_CHECK_COMPILE_FLAG([$CFLAGS],[],AC_MSG_ERROR(["Couldn't compile with given CFLAGS!"]))
    AC_SUBST(UNROLL_LOOPS,0)
fi

if test "$testcflags_set" = "no";
then
    TESTCFLAGS="$CFLAGS"
else
    TESTCFLAGS="$TESTCFLAGS $save_CFLAGS"
fi

################################################################################
# Assembly
################################################################################

#if test "$enable_assembly" = "yes";
#then
#    if test "x$asm_path" = "x";
#    then
#        AC_MSG_WARN([Currently no assembly available for $host. Disabling assembly...])
#        enable_assembly="no"
#    elif test "$flint_cv_abi" != "64";
#    then
#        AC_MSG_WARN([Currently no assembly available for 32-bit systems. Disabling assembly...])
#        enable_assembly="no"
#    fi
#fi
#
#if test "$enable_assembly" = "yes";
#then
#    GMP_INIT(config.m4)
#
#    GMP_PROG_M4
#    GMP_PATH_NM
#
#    CL_ASM_NOEXECSTACK
#    GMP_ASM_TEXT
#    GMP_ASM_DATA
#    GMP_ASM_LABEL_SUFFIX
#    GMP_ASM_GLOBL
#    GMP_ASM_GLOBL_ATTR
#    GMP_ASM_UNDERSCORE
#    GMP_ASM_RODATA
#    GMP_ASM_TYPE
#    GMP_ASM_SIZE
#    GMP_ASM_LSYM_PREFIX
#    GMP_ASM_ALIGN_LOG
#
#    case $host in
#        ARM64_PATTERN)
#            AC_DEFINE([FLINT_HAVE_ASSEMBLY_armv8],1,[Define if Arm v8 assembly is available])
#            case $host in
#                *-*-darwin*)
#                    GMP_INCLUDE_MPN(src/mpn_extras/arm64/darwin.m4)
#                    ;;
#                *)
#                    GMP_INCLUDE_MPN(src/mpn_extras/arm64/arm64-defs.m4)
#                    ;;
#            esac
#            ;;
#        X86_64_ADX_PATTERN)
#            AC_DEFINE([FLINT_HAVE_ASSEMBLY_x86_64_adx],1,[Define if x86_64 ADX assembly is available])
#            GMP_ASM_ALIGN_FILL_0x90
#            GMP_ASM_COFF_TYPE
#            GMP_INCLUDE_MPN(src/mpn_extras/x86_64/x86_64-defs.m4)
#            case $host in
#                *-*-darwin*)
#                    GMP_INCLUDE_MPN(src/mpn_extras/x86_64/darwin.m4)
#                    ;;
#            esac
#            ;;
#        *)
#            AC_MSG_ERROR([Error in assembly. Report to <https://github.com/vneiger/pml/issues/>.])
#    esac
#
#    GMP_FINISH
#else
#    # LSYM_PREFIX is still needed for x86-64 and ARM64 in longlong_asm_gnu.h
#    case $host in
#        ARM64_PATTERN | X86_64_PATTERN)
#            GMP_INIT(tmpconf)
#            GMP_ASM_TEXT
#            GMP_ASM_LABEL_SUFFIX
#            GMP_ASM_LSYM_PREFIX
#            rm -f $gmp_tmpconfigm4 $gmp_tmpconfigm4i $gmp_tmpconfigm4p
#            ;;
#    esac
#fi
#
#if test "$enable_assembly" = "yes";
#then
#    AC_SUBST(WANT_ASSEMBLY,1)
#    AC_SUBST(ASM_PATH,$asm_path)
#else
#    AC_SUBST(WANT_ASSEMBLY,0)
#    AC_SUBST(ASM_PATH,)
#fi

################################################################################
# parameters
################################################################################

#AC_CONFIG_LINKS([src/flint-mparam.h:src/mpn_extras/$param_path/flint-mparam.h],[],[param_path="$param_path"])
#AC_SUBST(PARAM_PATH,$param_path)

################################################################################
# substitutions and definitions
################################################################################

AC_SUBST(PML_VERSION_FULL, PML_VERSION)

# split version into major/minor/patch using POSIX variable substitutions
tail=PML_VERSION
PML_MAJOR=${tail%%.*}
tail=${tail#*.}
PML_MINOR=${tail%%.*}
tail=${tail#*.}
PML_PATCH=${tail%-*}

AC_SUBST(PML_MAJOR)
AC_SUBST(PML_MINOR)
AC_SUBST(PML_PATCH)
AC_SUBST(PML_MAJOR_SO)
AC_SUBST(PML_MINOR_SO)
AC_SUBST(PML_PATCH_SO)

AC_SUBST(GMP_LIB_PATH,$gmp_lib_path)
AC_SUBST(MPFR_LIB_PATH,$mpfr_lib_path)
AC_SUBST(FLINT_LIB_PATH,$flint_lib_path)
AC_SUBST(BLAS_LIB_PATH,$blas_lib_path)
AC_SUBST(GC_LIB_PATH,$gc_lib_path)
AC_SUBST(NTL_LIB_PATH,$ntl_lib_path)
AC_SUBST(MAPLE_LIB_PATH,$maple_lib_path)

AC_SUBST(LIB_CPPFLAGS,["-DPML_NOSTDIO -DPML_NOSTDARG"])

if test "$enable_shared" = "yes";
then
    AC_SUBST(SHARED,1)
else
    AC_SUBST(SHARED,0)
fi

if test "$enable_static" = "yes";
then
    AC_SUBST(STATIC,1)
else
    AC_SUBST(STATIC,0)
fi

if test "$enable_coverage" = "yes";
then
    AC_SUBST(COVERAGE,1)
    AC_DEFINE(PML_COVERAGE,1,[Define to enable coverage.])
else
    AC_SUBST(COVERAGE,0)
fi

if test "$enable_dependency_tracking" = "yes";
then
    AC_SUBST(WANT_DEPS,1)
else
    AC_SUBST(WANT_DEPS,0)
fi

if test "$enable_pretty_tests" = "yes";
then
    AC_DEFINE(PML_WANT_PRETTY_TESTS,1,[Define to enable pretty printing for tests.])
fi

if test "$with_ntl" = "yes";
then
    AC_SUBST(WANT_NTL,1)
else
    AC_SUBST(WANT_NTL,0)
fi

if test "$with_blas" = "yes";
then
    AC_DEFINE(PML_USES_BLAS,1,[Define to enable BLAS.])
fi

if test "$flint_big_endian" = "yes";
then
    AC_DEFINE(PML_BIG_ENDIAN,1,[Define if system is big endian.])
fi

if test "$enable_pthread" = "yes";
then
    AC_DEFINE(PML_USES_PTHREAD,1,[Define to enable the use of pthread.])
fi

if test "$enable_thread_safe" = "yes";
then
    AC_DEFINE(PML_USES_TLS,1,[Define to enable thread-local storage.])
fi

if test "$enable_reentrant" = "yes";
then
    AC_DEFINE(PML_REENTRANT,1,[Define to enable reentrant.])
fi

if test "$enable_assert" = "yes";
then
    AC_DEFINE(PML_WANT_ASSERT,1,[Define to enable use of asserts.])
fi

if test "$enable_gmp_internals" = "yes";
then
    AC_DEFINE(PML_WANT_GMP_INTERNALS,1,[Define to enable use of GMP internals.])
fi

if test "$ac_abs_confdir" = "`pwd`";
then
    AC_SUBST(IS_OUT_OF_TREE,0)
else
    AC_SUBST(IS_OUT_OF_TREE,1)
fi

################################################################################
# epilog
################################################################################

AC_CONFIG_FILES([Makefile pml.pc src/pml.h])
AC_OUTPUT

dnl Shorten the original help message
m4_cleardivert([HELP_BEGIN])dnl
m4_divert_push([HELP_BEGIN])dnl

if test "$ac_init_help" = "long";
then
    cat <<_ACEOF
\`configure' configures m4_ifset([AC_PACKAGE_STRING],
			[AC_PACKAGE_STRING],
			[this package]) to adapt to many kinds of systems.

Usage: $[0] [[OPTION]]... [[VAR=VALUE]]...

[To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print \`checking ...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for \`--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or \`..']

Installation directories:
]AS_HELP_STRING([--prefix=PREFIX],
  [install architecture-independent files in PREFIX [$ac_default_prefix]])
AS_HELP_STRING([--exec-prefix=EPREFIX],
  [install architecture-dependent files in EPREFIX [PREFIX]])[

By default, \`make install' will install all the files in
\`$ac_default_prefix/include' and \`$ac_default_prefix/lib' etc.
You can specify an installation prefix other than \`$ac_default_prefix'
using \`--prefix', for instance \`--prefix=\$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR            user executables [EPREFIX/bin]
  --libdir=DIR            object code libraries [EPREFIX/lib]
  --includedir=DIR        C header files [PREFIX/include]
_ACEOF

    cat <<\_ACEOF]
m4_divert_pop([HELP_BEGIN])dnl

dnl NOTE: There is supposed to be a missing `fi' here.
