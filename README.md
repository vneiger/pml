# PML (Polynomial Matrix Library)

Additions to the NTL and FLINT libraries, with a focus on univariate polynomial
matrices, structured matrices, and their applications.

Authors of NTL version: Seung Gyu Hyun, Vincent Neiger, Eric Schost

Authors of FLINT version, and current maintainers: Vincent Neiger, Eric Schost

Version 0.3

## Licensing

PML v0.3 is distributed under the GNU General Public License version 2.0
(GPL-2.0-or-later). This applies to all files in the directory of this README
as well as in all subdirectories. See the file COPYING for a copy of the
license.

PML v0.3 is heavily based on [NTL](https://libntl.org/) and
[FLINT](https://flintlib.org/). See the file `ntl-extras/COPYING_NTL` for NTL's
copyright notice. FLINT is distributed under LGPL 2.1 (GNU Lesser General
Public License), see `flint-extras/COPYING_FLINT` for the license.

## Installation

If compiling the NTL version, NTL should be installed, version >=11.3.1
required.

The FLINT version is work in progress, and currently compiles with the latest
git version of FLINT. Currently, some parts of flint-extras may not compile on
processors without avx512 instructions.

Building PML has mostly been tested on linux distributions. Installation relies
on "make", and the documentation relies on Doxygen.

Each directory should contain one or more .h file and subdirectories src/,
test/, timings/ (sometimes also tune/).

Running "make" at the root (of either ntl-extras or flint-extras) (re)builds
the entire library from scratch. Running "make doc" (after "make") builds a
Doxygen documentation that can be found in ROOT/include/html/index.html
(building the documentation requires having doxygen and graphviz installed).

The NTL version uses C++ code and .cpp files, whereas the FLINT version uses C
and .c files.

In `src/`,
 - "make" -> compile and install the header / object files.
 - "make clean" -> removes the object files.

In `test/`, test files should be called test-something.\{c,cpp\}
 - "make" or "make clean" -> removes all executables and check files (.chk)
 - "make all" -> compiles all test files
 - "make run" -> runs all executables, outputs results to something.chk
 - "make something.exe" -> compiles only test-something.\{c,cpp\}
 - "make something.chk" -> runs only test-something

In `timings/`, timing files should be called time-something.\{c,cpp\}
 - "make" or "make clean" -> removes all executables and data files (.dat)
 - "make all" -> compiles all timing files
 - "make run" -> runs all executables, outputs results to something.dat
 - "make something.exe" -> compiles only time-something.\{c,cpp\}
 - "make something.chk" -> runs only time-something

## For developers

Code style:
  - with vim and emacs, simply respect the modelines
  - no hard tabulation
  - soft tabulations = 4 spaces
  - scope delimiter \{ \} on their own line
  - scope delimiters are not required here when they are not required by C++
    (e.g. for one-line if/for)
  - when creating a new file, add the modelines at the bottom of the file:

```
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
```
