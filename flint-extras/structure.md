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
  macOS (with configure): ok  (+ #include <stdlib.h>)  ==> ok 
- fmpz_mat_extra.   ==> produits modulo p 
  multimodular matrix multiplication
  ~~macOS~~: cf nmod_mat_mul_small_modulus
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

- nmod_extra.    reprendre les fichiers 
  some small utility functions,
  some additions to machine_vectors.h,
  some CRT/multimod functions for nmod
     (using <= 4 primes, reduce requires AVX, CRT is ok without AVX)

  ~~macOS~~ (with configure), AVX test for small in CRT_CRT (+ #include <stdlib.h>), for all in CRT_reduce 
  test_multimod_CRT_CRT, nmod_multimod_CRT_CRT blocks for N=1 and 50 bits 

- nmod_mat_extra   (AVX flags TO DO)
  creating random matrices with particular properties
  row and column rotations
  matrix multiplication with AVX (small modulus)  (TO BE CLEANED + AVX CHECKS)
  PLUQ (uses AVX2 but from things already incorporated in FLINT: should compile as such, to be checked)
  left nullspace
  macOS (with configure): ok but AVX tests in mul_newdot.c and small_modulus.c
  check to be improved, test_nmod_mat_mul_newdot.bak test_nmod_mat_mul_small_modulus.bak

- nmod_mat_poly_extra
  types and basic functions for polynomial matrices stored as a vector of matrix coefficients
  mbasis (efficient, but can be made faster if input is m x n with n > m/2)
  macOS (with configure): ok 
  tests: todo (seems to be mostly ok)
  
- nmod_poly_extra
  evaluation and interpolation at points in geometric sequence
  there used to be sd_fft functions, to be reintroduced (tests are still here)
   (see commit 2ed64e4acd3137acf45495bc73131275a9f87a6d, Sat Aug 17 22:22:24 2024 +0200)
  macOS (with configure):  n_fft_evaluate_bench.bak (flint include pb)
  test: .bak for fft and tft, todo/to see for the rest 

- nmod_poly_mat_extra
  basic utilities for univariate polynomial matrices over nmod
  tools for handling reduced/weak Popov/Popov/Hermite forms (testing form, pivot degree, leading matrix...)
  some tools for printing degree matrix / leading matrix
  utils: permutations, rotations, truncate, shift, reverse, random, ...
  multiply: some middle-product and multiplication routines using various approaches (top urgent goal is to have here a fast FFT-based multiplication)
  macOS (with configure): almost ok, several .bak see calls to nmod_mat_mul_pml (+ #include <stdlib.h>)
  tests: todo (seems to be mostly ok)
  warning \# define FLINT_HAVE_NATIVE_mpn_add_n_sub_n  for nmod_poly_mat_mul_geometric.c


- nmod_vec_extra   (AVX flags TO DO)
  some basic routines
  some dot product routines with unbalanced sizes
  dot product with two vectors as input (TODO : see if still interesting, and clean)
  some dot product experiments for avx with moduli > 32 bits
  dot product multi: basically a vector matrix product, often faster than what is available through FLINT, but not integrated yet
  macOS (with configure): global AVX test in small_modulus.c and dot_product.c
  test integer_dot_product ok (others .bak)

---

# others

- sagemath_extra
utility functions for printing in SageMath-usable formats
macOS (with configure): seems to be ok, todo tests (printing works, strange?)  

---

