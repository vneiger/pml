
Each directory should contain a .h file and subdirectories src/, test/, timings/.

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

Code style:
  - with vim and emacs, simply respect the modelines
  - no hard tabulation
  - soft tabulations = 4 spaces
  - scope delimiter { } on their own line
  - may not insert scope delimiters when not necessary (e.g. for one-line if/for)
