STRUCTURE OF FLINT-EXTRAS
=========================

# misc -- compilation, documentation, compiled files

COPYING_FLINT
doxygen.conf
Makefile
include
lib

---

# fmpz -- multiprecision integers

- fmpz_extra
some CRT tools

- fmpz_mat_extra
multimodular matrix multiplication

- fmpz_poly_mat_extra
header for polynomial matrix multiplication (multimod | Waksman);
**no implementation yet**

---

# fmpz_mod -- arithmetic modulo multiprecision integers

**no implementation yet**

- fmpz_mod_mat_extra
matrices over Z/nZ, large n
no implementation yet

- fmpz_mod_poly_mat
univariate polynomial matrices over Z/nZ,
no implementation yet

---

# nmod -- arithmetic modulo word-size integers

- nmod_extra
  some small utility functions,
  some additions to machine_vectors.h,
  some CRT/multimod functions for nmod
     (using <= 4 primes, reduce requires AVX, CRT is ok without AVX)

  macOS: AVX test for small in CRT_CRT (#include <stdlib.h>), for all in CRT_reduce 

  ​	tests: test_multimod_CRT_CRT stops after i=1 in the main loop 

- nmod_mat_extra   (AVX flags TO DO)
creating random matrices with particular properties
row and column rotations
matrix multiplication with AVX (small modulus)  (TO BE CLEANED + AVX CHECKS)
PLUQ (uses AVX2 but from things already incorporated in FLINT: should compile as such, to be checked)
left nullspace

- nmod_mat_poly_extra
types and basic functions for polynomial matrices stored as a vector of matrix coefficients
mbasis (efficient, but can be made faster if input is m x n with n > m/2)

- nmod_poly_extra
evaluation and interpolation at points in geometric sequence
there used to be sd_fft functions, to be reintroduced (tests are still here)
   (see commit 2ed64e4acd3137acf45495bc73131275a9f87a6d, Sat Aug 17 22:22:24 2024 +0200)

- nmod_poly_mat_extra
basic utilities for univariate polynomial matrices over nmod
tools for handling reduced/weak Popov/Popov/Hermite forms (testing form, pivot degree, leading matrix...)
some tools for printing degree matrix / leading matrix
utils: permutations, rotations, truncate, shift, reverse, random, ...
multiply: some middle-product and multiplication routines using various approaches (top urgent goal is to have here a fast FFT-based multiplication)


- nmod_vec_extra   (AVX flags TO DO)
some basic routines
some dot product routines with unbalanced sizes
dot product with two vectors as input (TODO : see if still interesting, and clean)
some dot product experiments for avx with moduli > 32 bits
dot product multi: basically a vector matrix product, often faster than what is available through FLINT, but not integrated yet

---

# others

- sagemath_extra
utility functions for printing in SageMath-usable formats

---

