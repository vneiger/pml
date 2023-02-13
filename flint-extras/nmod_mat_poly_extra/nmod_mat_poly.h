/**
 * \file nmod_mat_poly.h
 * \brief Main header for univariate polynomial with matrix coefficients modulo word-size prime
 * \version 0.0
 * \date 2023-02-13
 *
 * This is the main header for functions for univariate polynomials with
 * coefficients that are matrices over a finite field Z/pZ for "small" p
 * (word-size, so that we use a representation with flint's ``nmod``). Only 
 * operations for which this `nmod_mat_poly` format is best suited are
 * provided. For other operations it is suggested to rely the conversions to and
 * from `nmod_poly_mat` format and use the relevant functions for that format.
 *
 * \todo benchmark performance
 * \todo test for memory leaks
 * \todo Note: all parameters are supposed init
 *
 */

#ifdef NMOD_MAT_POLY_INLINES_C
#define NMOD_MAT_POLY_INLINE FLINT_DLL
#else
#define NMOD_MAT_POLY_INLINE static __inline__
#endif

#ifndef NMOD_MAT_POLY_H
#define NMOD_MAT_POLY_H

#include <flint/perm.h>
#include <flint/nmod_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/* struct, typedefs                                           */
/*------------------------------------------------------------*/

/** Struct for matrix polynomials.
 *
 * Storage is a dynamic array of matrices `nmod_mat`. The maximum number of
 * coefficients is `alloc` and the actual number of coefficients (which is the
 * degree plus 1) is `length`. The number of rows and columns are `r` and `c`.
 * The modulus is stored as an `nmod_t`. In the provided functions, e.g. for
 * modifying a coefficient, it is not checked that the dimensions of the
 * modified coefficient are indeed `r x c`.
 */
typedef struct
{
    nmod_mat_struct * coeffs; /**< array of coefficients */
    slong alloc;              /**< allocated length */
    slong length;             /**< actual length */
    slong r;                  /**< number of rows */
    slong c;                  /**< number of columns */
    nmod_t mod;               /**< modulus */
} nmod_mat_poly_struct;

/* nmod_mat_poly_t allows reference-like semantics for nmod_mat_poly_struct */
typedef nmod_mat_poly_struct nmod_mat_poly_t[1];

/*------------------------------------------------------------*/
/* memory management                                          */
/*------------------------------------------------------------*/

/** @name Memory management for matrix polynomials
 */
//@{

/** Initialises `matp`. It will have dimensions `r x c` and coefficients modulo
 * `n`. */
FLINT_DLL void nmod_mat_poly_init(nmod_mat_poly_t matp,
                                  slong r, slong c,
                                  mp_limb_t n);

/** Initialises `matp`. It will have dimensions `r x c` and coefficients modulo
 * `n`. The caller supplies a precomputed inverse limb generated by
 * `n_preinvert_limb()`. */
FLINT_DLL void nmod_mat_poly_init_preinv(nmod_mat_poly_t matp,
                                         slong r, slong c,
                                         mp_limb_t n, mp_limb_t ninv);

/** Initialises `matp`. It will have dimensions `r x c` and coefficients modulo
 * `n`. Up to alloc coefficients may be stored in `matp`. Implementation note:
 * the `alloc` matrix coefficients are not initialized. */
FLINT_DLL void nmod_mat_poly_init2(nmod_mat_poly_t matp,
                                   slong r, slong c,
                                   mp_limb_t n,
                                   slong alloc);

/** Initialises `matp`. It will have dimensions `r x c` and coefficients
 * modulo~`n`. The caller supplies a precomputed inverse limb generated by
 * n_preinvert_limb(). Up to alloc coefficients may be stored in `matp`.
 * Implementation note: the `alloc` matrix coefficients are not initialized. */
FLINT_DLL void nmod_mat_poly_init2_preinv(nmod_mat_poly_t matp, 
                                          slong r, slong c,
                                          mp_limb_t n, mp_limb_t ninv,
                                          slong alloc);

/** Reallocates `matp` to the given length `alloc`. If the current `length` is
 * more than `alloc`, the polynomial is truncated and normalised. If `alloc` is
 * zero, the polynomial is cleared. */
FLINT_DLL void nmod_mat_poly_realloc(nmod_mat_poly_t matp, slong alloc);

/** Clears `matp` and releases any memory it used. The matrix polynomial cannot be
used again until it is initialised. */
FLINT_DLL void nmod_mat_poly_clear(nmod_mat_poly_t matp);

/** Ensures `matp` has space for at least `alloc` coefficients. This function
 * only ever grows the allocated space, so no data loss can occur. */
FLINT_DLL void nmod_mat_poly_fit_length(nmod_mat_poly_t matp, slong alloc);

/** Initialises `matp` using an already initialised modulus `mod`. */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_init_mod(nmod_mat_poly_t matp,
                       slong r, slong c,
                       const nmod_t mod)
{
    matp->coeffs = NULL;
    matp->alloc = 0;
    matp->length = 0;
    matp->r = r;
    matp->c = c;
    matp->mod = mod;
}

//NMOD_MAT_POLY_INLINE
//void nmod_mat_poly_set_mod(nmod_mat_poly_t matp, const nmod_t mod)
//{
//    matp->mod = mod;
//}

/** Sets the length of `matp` to `length`. If `length < matp->length` then the
 * matrix coefficients of `matp` beyond `length` are cleared; otherwise, the
 * matrix coefficients of `matp` are initialized up to `length`. Note: `matp`
 * may not be normalized after this. The provided `length` must be less than
 * `matp->alloc`. */
NMOD_MAT_POLY_INLINE void
_nmod_mat_poly_set_length(nmod_mat_poly_t matp, slong length)
{
    if (matp->length > length)
        for (slong i = length; i < matp->length; i++)
            nmod_mat_clear(matp->coeffs + i); 
    else
        for (slong i = matp->length; i < length; i++)
            nmod_mat_init(matp->coeffs + i, matp->r, matp->c, matp->mod.n); 
    matp->length = length;
}

/** Normalises a matrix polynomial `matp` so that the top coefficient, if there
 * is one at all, is not zero. */
NMOD_MAT_POLY_INLINE
void _nmod_mat_poly_normalise(nmod_mat_poly_t matp)
{
    while (matp->length && nmod_mat_is_zero(matp->coeffs + matp->length - 1))
    {
        nmod_mat_clear(matp->coeffs + matp->length - 1);
        matp->length--;
    }
}

//@} // doxygen group:  Memory management for matrix polynomials

/*------------------------------------------------------------*/
/* Zero and Identity                                          */
/*------------------------------------------------------------*/

/** @name Zero and Identity
 * \todo TODO doc
 */
//@{

/** Sets `matp` to the zero matrix polynomial */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_zero(nmod_mat_poly_t matp)
{
   _nmod_mat_poly_set_length(matp, 0);
}

/** \def nmod_mat_poly_is_zero(matp)
 * Tests whether `matp` is the zero matrix polynomial
 */
#define nmod_mat_poly_is_zero(matp) \
    ((matp)->length == 0)

/** Sets `matp` to the identity matrix polynomial */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_one(nmod_mat_poly_t matp)
{
    nmod_mat_poly_fit_length(matp, 1);
    nmod_mat_one(matp->coeffs + 0);
    _nmod_mat_poly_set_length(matp, 1);
}

/** Tests whether `matp` is the identity matrix polynomial */
NMOD_MAT_POLY_INLINE int
nmod_mat_poly_is_one(const nmod_mat_poly_t matp)
{
    return (matp->length) == 1 && (nmod_mat_is_one(matp->coeffs + 0));
}

//@} // doxygen group:  Zero and Identity


/*------------------------------------------------------------*/
/* Accessing struct info and coefficients                     */
/*------------------------------------------------------------*/

/** @name Accessing struct info and matrix coefficients
 * \todo TODO doc
 */
//@{

/** Returns the number of rows of `matp`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_nrows(const nmod_mat_poly_t matp)
{
    return matp->r;
}

/** Returns the number of cols of `matp`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_ncols(const nmod_mat_poly_t matp)
{
    return matp->c;
}

/** Returns the length of `matp`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_length(const nmod_mat_poly_t matp)
{
    return matp->length;
}

/** Returns the degree of `matp`. By convention the zero matrix polynomial has
 * degree `-1`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_degree(const nmod_mat_poly_t matp)
{
    return matp->length - 1;
}

/** Leading matrix coefficient of `matp`. */
NMOD_MAT_POLY_INLINE nmod_mat_struct *
nmod_mat_poly_lead(const nmod_mat_poly_t matp)
{
    if (matp->length)
        return matp->coeffs + (matp->length - 1);
    else
        return (nmod_mat_struct *)NULL;
}

/** \def nmod_mat_poly_get_coeff_ptr(matp, n)
 * Returns a reference to the coefficient of `x**n` in the matrix polynomial
 * `matp`. This function is provided so that individual coefficients can be
 * accessed and operated on by functions in the `nmod_mat` module. This
 * function does not make a copy of the data, but returns a reference
 * `nmod_mat_struct *` to the actual coefficient. Returns `NULL` when `n`
 * exceeds the degree of the polynomial.
 */
#define nmod_mat_poly_get_coeff_ptr(matp, n) \
    ((n) < (matp)->length ? (matp)->coeffs + (n) : NULL)

/** \def nmod_mat_poly_lead(const nmod_mat_poly_t poly)
 * Returns a reference to the leading coefficient of the matrix polynomial, as
 * an `nmod_mat_struct *`. This function is provided so that the leading
 * coefficient can be easily accessed and operated on by functions in the
 * `nmod_mat` module. This function does not make a copy of the data, but
 * returns a reference to the actual coefficient.  Returns `NULL` when the
 * polynomial is zero.
 */
#define nmod_mat_poly_lead(matp) \
    ((matp)->length ? (matp)->coeffs + (matp)->length - 1 : NULL)

/** \def nmod_mat_poly_entry(matp,k,i,j)
 * Directly accesses the entry in the coefficient of `matp` of degree `k`, in
 * row `i` and column `j` (indexed from zero). No bounds checking is performed.
 * This macro can be used both for reading and writing coefficients.
 */
#define nmod_mat_poly_entry(matp,k,i,j) \
    (((matp)->coeffs + (k))->rows[(i)][(j)])

/** Get the entry at row `i` and column `j` in the coefficient of
 * degree `k` of the matrix polynomial `matp`. */
NMOD_MAT_POLY_INLINE mp_limb_t
nmod_mat_poly_get_entry(const nmod_mat_poly_t matp,
                        slong k, slong i, slong j)
{
   return (matp->coeffs + k)->rows[i][j];
}

/** Return a pointer to the entry at row `i` and column `j` of the coefficient
 * of degree `k` of the matrix polynomial `matp` */
NMOD_MAT_POLY_INLINE mp_limb_t *
nmod_mat_poly_entry_ptr(const nmod_mat_poly_t matp,
                        slong k, slong i, slong j)
{
   return (matp->coeffs + k)->rows[i] + j;
}

/** Set to `x` the entry at row `i` and column `j` in the coefficient of degree
 * `k` of the matrix polynomial `matp`. */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_set_entry(nmod_mat_poly_t matp,
                        slong k, slong i, slong j,
                        mp_limb_t x)
{
    nmod_mat_poly_entry(matp, k, i, j) = x;
}

//@} // doxygen group:  Accessing struct info and matrix coefficients



/** @name Truncate, shift, reverse.
 * \todo TODO
 */
//@{

/** Truncates `matp` to the given `order` and normalises it. If `order` is
 * greater than the current length of `matp`, then nothing happens. */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_truncate(nmod_mat_poly_t matp, slong order)
{
    if (matp->length > order)
    {
        for (slong i = order; i < matp->length; i++)
            nmod_mat_clear(matp->coeffs + i);
        matp->length = order;
        _nmod_mat_poly_normalise(matp);
    }  
}

//@} // doxygen group:  Truncate, shift, reverse




























































/** void nmod_mat_poly_naive_mul_coef(nmod_mat_t res,
 *				       const nmod_mat_poly_t A,
 *				       const nmod_mat_poly_t B, slong k);
 *
 * A = sum^{deg_A}_{i=0} a_i x^i and B = sum^{deg_B}_{i=0} b_i x^i
 * Compute the coefficient k of the product AB
 * C = AB = sum^{deg_A + deg_B}_{i=0} c_i x^i
 * res = c_k = sum^{k}_{i=0} a_i b_{k-i}
 *
 */
void nmod_mat_poly_naive_mul_coef(nmod_mat_t res,
                                  const nmod_mat_poly_t A,
                                  const nmod_mat_poly_t B,
                                  slong k);

void nmod_mat_poly_print(const nmod_mat_poly_t A);


#ifdef __cplusplus
}
#endif

#endif /* NMOD_MAT_POLY_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
