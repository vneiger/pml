# ntl-extra
Additions to Shoup's NTL, with a focus on univariate polynomial matrices, constant structured matrices, and their applications.

Version 0.0

# Installation

NTL should be installed.

For the moment, this has only been tested on some linux distributions, with NTL installed in default location. Installation relies on "make" and documentation is built via Doxygen.

Each directory should contain one or more .h file and subdirectories src/, test/, timings/ (sometimes also tune/).

Running "make" at the root, where this README is, (re)builds the entire library from scratch. Running "make doc" (after "make") builds a Doxygen documentation that can be found in ROOT/include/html/index.html

In src/, "make" -> compile and install the header / object files.
         "make clean" -> removes the object files.

In test/, test files should be called test-something.cpp
          "make" or "make clean" -> removes all executables and check files (.chk)
          "make all" -> compiles all test files
          "make run" -> runs all executables, outputs results to something.chk
                        something.chk should be empty if all goes well
          "make something.exe" -> compiles only test-something.cpp
          "make something.chk" -> runs only test-something

In timings/, timing files should be called time-something.cpp
          "make" or "make clean" -> removes all executables and data files (.dat)
          "make all" -> compiles all timing files
          "make run" -> runs all executables, outputs results to something.dat
          "make something.exe" -> compiles only time-something.cpp
          "make something.chk" -> runs only time-something

# For developers

Code style:
  - with vim and emacs, simply respect the modelines
  - no hard tabulation
  - soft tabulations = 4 spaces
  - scope delimiter { } on their own line
  - scope delimiters are not required here when they are not required by C++ (e.g. for one-line if/for)
