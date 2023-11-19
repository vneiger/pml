# PML (Polynomial Matrix Library)

Additions to Victor Shoup's NTL, with a focus on univariate polynomial matrices, structured matrices, and their applications.

Authors: Seung Gyu Hyun, Vincent Neiger, Éric Schost

Version 0.2

# Licensing

PML v0.2 is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). This applies to all files in the directory of this README as well as in all subdirectories. See the file COPYING for a copy of the license.

PML v0.2 is heavily based on [NTL](https://www.shoup.net/ntl). See the file COPYING_NTL for NTL's copyright notice.

# Installation

NTL should be installed, version >=11.3.1 required.

For the moment, this has only been tested on some linux distributions, with NTL installed in default location. Installation relies on "make" and documentation is built via Doxygen.

Each directory should contain one or more .h file and subdirectories src/, test/, timings/ (sometimes also tune/).

Running "make" at the root, where this README is, (re)builds the entire library from scratch. Running "make doc" (after "make") builds a Doxygen documentation that can be found in ROOT/include/html/index.html  (building the documentation requires having doxygen and graphviz installed).

In `src/`,
 - "make" -> compile and install the header / object files.
 - "make clean" -> removes the object files.

In `test/`, test files should be called test-something.cpp
 - "make" or "make clean" -> removes all executables and check files (.chk)
 - "make all" -> compiles all test files
 - "make run" -> runs all executables, outputs results to something.chk
 - "make something.exe" -> compiles only test-something.cpp
 - "make something.chk" -> runs only test-something

In `timings/`, timing files should be called time-something.cpp
 - "make" or "make clean" -> removes all executables and data files (.dat)
 - "make all" -> compiles all timing files
 - "make run" -> runs all executables, outputs results to something.dat
 - "make something.exe" -> compiles only time-something.cpp
 - "make something.chk" -> runs only time-something

# For developers

Code style:
  - with vim and emacs, simply respect the modelines
  - no hard tabulation
  - soft tabulations = 4 spaces
  - scope delimiter { } on their own line
  - scope delimiters are not required here when they are not required by C++ (e.g. for one-line if/for)
  - when creating a new file, add the modelines at the bottom of the file:
```
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
```
