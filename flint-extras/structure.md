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
some CRT/multimod functions for nmod (using <= 4 primes)

nmod_mat_extra
nmod_mat_poly_extra
nmod_poly_extra
nmod_poly_mat_extra
nmod_vec_extra

---

# others

sagemath_extra
