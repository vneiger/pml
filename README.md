# PML (Polynomial Matrix Library)

Additions to the NTL and FLINT libraries, with a focus on univariate polynomial
matrices, structured matrices, and their applications.

Version 0.5

Warning: the FLINT-based part of PML (flint-extras folder) is work in progress.
Currently, it requires FLINT version 3.4.0 (November 2025) or later.
 
## Authors

Current maintainers:

 - Vincent Neiger
 - Eric Schost
 - Gilles Villard 

Contributors:

 - Seung Gyu Hyun
 - Kevin Tran

## Licensing

PML v0.5 is distributed under the GNU General Public License version 2.0
(GPL-2.0-or-later). This applies to all files in the directory of this README
as well as in all subdirectories. See the file COPYING for a copy of the
license.

PML v0.5 is heavily based on [NTL](https://libntl.org/) and
[FLINT](https://flintlib.org/). See the file `ntl-extras/COPYING_NTL` for
NTL's copyright notice. FLINT is distributed under LGPL 3.0 (GNU Lesser General
Public License), see `flint-extras/COPYING_FLINT` for the license.

## Citing PML

```
@inproceedings{HyunNeigerSchost2019,
  author = {Hyun, S. G. and Neiger, V. and Schost, {\'E}.},
  title = {Implementations of Efficient Univariate Polynomial Matrix Algorithms and Application to Bivariate Resultants},
  year = {2019},
  isbn = {9781450360845},
  publisher = {ACM},
  doi = {10.1145/3326229.3326272},
  booktitle = {Proceedings ISSAC 2019},
  pages = {235--242},
  numpages = {8},
  keywords = {algorithms, implementation, resultant, polynomial matrices},
  note = {\url{https://github.com/vneiger/pml}},
}

@manual{pml,
    key = {PML},
    author = {The {PML} team},
    title = {{PML}: {P}olynomial {M}atrix {L}ibrary},
    year = {2025},
    note = {Version 0.5, \url{https://github.com/vneiger/pml}}
}
```

## Installation (FLINT-based version)

The FLINT-based version is work in progress, and currently compiles with FLINT
3.4.0 or later versions, including the current development version.

The build system is similar to (and directly derived from) FLINT's one. In
short, for a standard configuration, the following steps should work:
```
./bootstrap.sh
./configure
make
```

It is also recommended to run `make check` to verify that all functionalities
work as expected. If not, please submit a bug report.

## Installation (NTL-based version)

If compiling the NTL version, NTL should be installed, version at least 11.3.1
required. Building PML/NTL has mostly been tested on linux distributions.
Installation relies on "make", and the documentation relies on Doxygen. Each
directory should contain one or more .h file and subdirectories `src`, `test`,
`timings` (sometimes also `tune`).

Running "make" at the root (of either ntl-extras or flint-extras) (re)builds
the entire library from scratch. Running "make doc" (after "make") builds a
Doxygen documentation that can be found in `ROOT/include/html/index.html`
(building the documentation requires having doxygen and graphviz installed).

The NTL version uses C++ code and .cpp files, whereas the FLINT version uses C
and .c files.

In `src`,

 - "make": compile the object files.
 - "make clean": clean compilation products.

In `test`, test files should be called test-something.\{c,cpp\}

 - "make" or "make clean": removes all executables and check files (.chk)
 - "make all": compiles all test files
 - "make run": runs all executables, outputs results to something.chk
 - "make something.exe": compiles only test-something.\{c,cpp\}
 - "make something.chk": runs only test-something

In `timings`, timing files should be called time-something.\{c,cpp\}

 - "make" or "make clean": removes all executables and data files (.dat)
 - "make all": compiles all timing files
 - "make run": runs all executables, outputs results to something.dat
 - "make something.exe": compiles only time-something.\{c,cpp\}
 - "make something.chk": runs only time-something

## For developers

Code style:

  - no hard tabulation
  - soft tabulations, 4 spaces
  - scope delimiter \{ \} on their own line
  - scope delimiters are not required here when they are not required by C/C++
    (e.g. for one-line if or for)
