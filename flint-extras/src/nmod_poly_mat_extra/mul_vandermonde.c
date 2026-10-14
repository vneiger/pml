/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_multiply.h"

/*
    Products and middle products of polynomial matrices by evaluation and
    interpolation at small points, with every step a product by a
    Vandermonde matrix or its inverse (Doliskani, Giorgi, Lebreton, Schost):

    - vandermonde1 uses the points 0, 1, ..., P-1;
    - vandermonde2 uses the pairs of points +-1, +-2, ..., +-P, through the
      even and odd parts f(x) = f_e(x^2) + x f_o(x^2), so that its
      Vandermonde matrices are on the points 1, 4, ..., P^2 and have half
      the size.

    The coefficients of the entries of a polynomial matrix are laid out as
    the columns of a constant matrix ("packing"): column u*c + v holds the
    coefficients of entry (u, v), so that a Vandermonde matrix times it
    gives, in its row t, the evaluation of the whole polynomial matrix at
    the t-th point, as an r x c matrix stored row major. The products at
    the points are then nmod_mat_mul calls on those rows, borrowed in
    place.

    The middle products are the transposes of products (transposition
    principle). With S the operand of smaller length ls, L the other one
    and N = nhi - nlo, the coefficients nlo, ..., nhi-1 of S*L are

        mid_t = sum_j rev(S)_j z_{t+j},    z_s = L_{s + nlo - ls + 1},

    i.e. the transpose of the multiplication y -> rev(S) * y on vectors y
    of length N, applied to a window z of length ls + N - 1 of L. The
    product algorithm for that multiplication is

        y -> iV . diag(V_ls rev(S)) . V_N y       (P = ls + N - 1 points)

    and its transpose

        z -> V_N^T . diag(V_ls rev(S)) . iV^T z,

    with the same products of constant matrices, transposed, and the same
    pointwise products. Compared with the full product followed by a
    truncation, the number of points goes from len1 + len2 - 1 down to
    ls + nhi - nlo - 1: for the balanced middle product of the approximant
    basis algorithms (len1 = nlo + 1, len2 = nhi = 2 nlo), from 3 nlo down
    to 2 nlo. No constraint on len1, len2, nlo, nhi is needed for
    correctness: coefficients of L outside [0, len) are taken as zero.
*/

/*------------------------------------------------------------*/
/* constant matrices                                          */
/*------------------------------------------------------------*/

/* V = (pts[i]^j), 0 <= i < V->r, 0 <= j < V->c */
static void
_vdm_powers(nmod_mat_t V, nn_srcptr pts, nmod_t mod)
{
    slong i, j;
    if (V->c == 0)
        return;
    for (i = 0; i < V->r; i++)
    {
        nn_ptr row = nmod_mat_entry_ptr(V, i, 0);
        row[0] = (mod.n == 1) ? 0 : 1;
        for (j = 1; j < V->c; j++)
            row[j] = nmod_mul(row[j-1], pts[i], mod);
    }
}

/*
    Lag (npts x npts, npts = Lag->r): row i holds the coefficients of
    the Lagrange polynomial L_i = prod_{k != i} (x - pts[k]) / (pts[i] - pts[k]),
    i.e. Lag is the transpose of the inverse of the square Vandermonde
    matrix on pts. With M = prod_k (x - pts[k]), L_i = (M / (x - pts[i])) /
    M'(pts[i]): the master polynomial by a product tree, its derivative at the
    points by multipoint evaluation, then one synthetic division of M per row
    -- the only quadratic part. The points must be distinct with invertible
    differences (e.g. distinct modulo a prime).
*/
static void
_vdm_lagrange(nmod_mat_t Lag, nn_srcptr pts, nmod_t mod)
{
    const slong npts = Lag->r;
    slong i;
    nn_ptr M = _nmod_vec_init(npts + 1);
    nn_ptr w = _nmod_vec_init(2 * npts);
    nn_ptr winv = w + npts;

    _nmod_poly_product_roots_nmod_vec(M, pts, npts, mod);

    /* w[i] = M'(pts[i]) = prod_{k != i} (pts[i] - pts[k]); the derivative
       of M has length npts, written over the first npts entries of winv */
    _nmod_poly_derivative(winv, M, npts + 1, mod);
    _nmod_poly_evaluate_nmod_vec(w, winv, npts, pts, npts, mod);
    _nmod_vec_invert(winv, w, npts, mod);

    for (i = 0; i < npts; i++)
    {
        nn_ptr q = nmod_mat_entry_ptr(Lag, i, 0);
        _nmod_poly_div_root(q, M, npts + 1, pts[i], mod);
        _nmod_vec_scalar_mul_nmod(q, q, npts, winv[i], mod);
    }

    _nmod_vec_clear(M);
    _nmod_vec_clear(w);
}

/*------------------------------------------------------------*/
/* layout of polynomial matrices as constant matrices         */
/*------------------------------------------------------------*/

/*
    M (M->r x r*c): column u*c + v holds the coefficients of indices
    first, first + step, ..., first + (M->r - 1)*step of entry (u, v) of
    pmat, those outside [0, length) being zero. step may be negative
    (reversal).
*/
static void
_vdm_pack(nmod_mat_t M, const nmod_poly_mat_t pmat, slong first, slong step)
{
    const slong r = pmat->r, c = pmat->c;
    slong u, v, k;

    for (u = 0; u < r; u++)
        for (v = 0; v < c; v++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(pmat, u, v);
            const slong e = u * c + v;
            for (k = 0; k < M->r; k++)
            {
                const slong j = first + k * step;
                nmod_mat_entry(M, k, e) = (j >= 0 && j < f->length) ? f->coeffs[j] : 0;
            }
        }
}

/*
    The unpacking functions write the r x c result, which is the leading
    submatrix of pmat (as for the products of FLINT, rows and columns of
    pmat beyond are left untouched).
*/

/* the entries of the r x c leading submatrix of pmat are given length
   len, zero coefficients */
static void
_vdm_unpack_init(nmod_poly_mat_t pmat, slong r, slong c, slong len)
{
    slong u, v;
    for (u = 0; u < r; u++)
        for (v = 0; v < c; v++)
        {
            nmod_poly_struct * f = nmod_poly_mat_entry(pmat, u, v);
            nmod_poly_fit_length(f, len);
            _nmod_vec_zero(f->coeffs, len);
            _nmod_poly_set_length(f, len);
        }
}

/* coefficient first + k*step of entry (u, v) of pmat <- M[k][u*c + v],
   for u < r, v < c and the k whose index is below the length of the entry */
static void
_vdm_unpack(nmod_poly_mat_t pmat, slong r, slong c, const nmod_mat_t M,
            slong first, slong step)
{
    slong u, v, k;

    for (u = 0; u < r; u++)
        for (v = 0; v < c; v++)
        {
            nmod_poly_struct * f = nmod_poly_mat_entry(pmat, u, v);
            const slong e = u * c + v;
            for (k = 0; k < M->r && first + k * step < f->length; k++)
                f->coeffs[first + k * step] = nmod_mat_entry(M, k, e);
        }
}

static void
_vdm_unpack_normalise(nmod_poly_mat_t pmat, slong r, slong c)
{
    slong u, v;
    for (u = 0; u < r; u++)
        for (v = 0; v < c; v++)
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, u, v));
}

/*------------------------------------------------------------*/
/* products of constant matrices                              */
/*------------------------------------------------------------*/

/* C = A[:C->r, :k] * B[:k, :], with C->r <= A->r, C->c == B->c and
   k <= min(A->c, B->r); zero when k == 0 */
static void
_vdm_mul_window(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, slong k)
{
    const slong rows = C->r;

    if (rows == 0 || C->c == 0)
        return;
    if (k == 0)
    {
        nmod_mat_zero(C);
        return;
    }

    if (rows == A->r && k == A->c && k == B->r)
        nmod_mat_mul(C, A, B);
    else
    {
        nmod_mat_t Aw, Bw;
        nmod_mat_window_init(Aw, A, 0, 0, rows, k);
        nmod_mat_window_init(Bw, B, 0, 0, k, B->c);
        nmod_mat_mul(C, Aw, Bw);
        nmod_mat_window_clear(Aw);
        nmod_mat_window_clear(Bw);
    }
}

/*
    The products at the points: for each row t of the three matrices,
    C_t = A_t * B_t, where row t of Aval (resp. Bval, Cval) is an m x n
    (resp. n x p, m x p) matrix stored row major.
*/
static void
_vdm_pointwise(nmod_mat_t Cval, const nmod_mat_t Aval, const nmod_mat_t Bval,
               slong m, slong n, slong p)
{
    const slong npts = Aval->r;
    slong t;

    if (m == 0 || p == 0)
        return;
    if (n == 0)
    {
        nmod_mat_zero(Cval);
        return;
    }

    for (t = 0; t < npts; t++)
    {
        nmod_mat_struct At, Bt, Ct;

        At.entries = nmod_mat_entry_ptr(Aval, t, 0);
        At.r = m; At.c = n; At.stride = n; At.mod = Aval->mod;

        Bt.entries = nmod_mat_entry_ptr(Bval, t, 0);
        Bt.r = n; Bt.c = p; Bt.stride = p; Bt.mod = Aval->mod;

        Ct.entries = nmod_mat_entry_ptr(Cval, t, 0);
        Ct.r = m; Ct.c = p; Ct.stride = p; Ct.mod = Aval->mod;

        nmod_mat_mul(&Ct, &At, &Bt);
    }
}

/* row t of M <- c[t] * (row t of M) */
static void
_vdm_scale_rows(nmod_mat_t M, nn_srcptr c, nmod_t mod)
{
    slong t;
    for (t = 0; t < M->r; t++)
        _nmod_vec_scalar_mul_nmod(nmod_mat_entry_ptr(M, t, 0),
                                  nmod_mat_entry_ptr(M, t, 0), M->c, c[t], mod);
}

/* (X, Y) <- (X + Y, X - Y), rowwise */
static void
_vdm_butterfly(nmod_mat_t X, nmod_mat_t Y, nmod_t mod)
{
    slong t, j;
    for (t = 0; t < X->r; t++)
    {
        nn_ptr x = nmod_mat_entry_ptr(X, t, 0);
        nn_ptr y = nmod_mat_entry_ptr(Y, t, 0);
        for (j = 0; j < X->c; j++)
        {
            const ulong a = x[j], b = y[j];
            x[j] = nmod_add(a, b, mod);
            y[j] = nmod_sub(a, b, mod);
        }
    }
}

/*
    vandermonde2: the evaluations of the polynomial matrix whose packed
    even and odd parts are the coefficients first_e + k*step and
    first_o + k*step of pmat, k < s_e (resp. s_o), at the points +x_t:
    Plus at +x_t, Minus at -x_t, both V->r x (pmat->r * pmat->c), where V
    is the Vandermonde matrix on the squares of the x_t.
*/
static void
_vdm2_eval(nmod_mat_t Plus, nmod_mat_t Minus, const nmod_mat_t V, nn_srcptr x,
           const nmod_poly_mat_t pmat, slong s_e, slong first_e,
           slong s_o, slong first_o, slong step, nmod_t mod)
{
    nmod_mat_t T;

    nmod_mat_init(T, s_e, pmat->r * pmat->c, mod.n);
    _vdm_pack(T, pmat, first_e, step);
    _vdm_mul_window(Plus, V, T, s_e);
    nmod_mat_clear(T);

    nmod_mat_init(T, s_o, pmat->r * pmat->c, mod.n);
    _vdm_pack(T, pmat, first_o, step);
    _vdm_mul_window(Minus, V, T, s_o);
    nmod_mat_clear(T);

    _vdm_scale_rows(Minus, x, mod);     /* x f_o(x^2) */
    _vdm_butterfly(Plus, Minus, mod);   /* f(x), f(-x) */
}

/* inv2x[t] = 1 / (2 x_t), t < npts */
static void
_vdm2_inv_2x(nn_ptr inv2x, nn_srcptr x, slong npts, nmod_t mod)
{
    nn_ptr t2x = _nmod_vec_init(npts);
    _nmod_vec_add(t2x, x, x, npts, mod);
    _nmod_vec_invert(inv2x, t2x, npts, mod);
    _nmod_vec_clear(t2x);
}

/* x_t = t + 1 and its square, for t < npts */
static void
_vdm2_points(nn_ptr x, nn_ptr x2, slong npts, nmod_t mod)
{
    slong t;
    for (t = 0; t < npts; t++)
    {
        x[t] = nmod_set_ui(t + 1, mod);
        x2[t] = nmod_mul(x[t], x[t], mod);
    }
}

/*------------------------------------------------------------*/
/* products                                                   */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void nmod_poly_mat_mul_vandermonde1(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b)
{
    const slong m = a->r, n = a->c, p = b->c;
    const slong ellA = nmod_poly_mat_max_length(a);
    const slong ellB = nmod_poly_mat_max_length(b);
    slong npts, t;
    nmod_t mod;
    nn_ptr pts;
    nmod_mat_t V, Lag, iV, Ta, Tb, Aval, Bval, Cval, Cc;

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }

    nmod_init(&mod, a->modulus);
    npts = ellA + ellB - 1;

    pts = _nmod_vec_init(npts);
    nmod_mat_init(V, npts, FLINT_MAX(ellA, ellB), mod.n);
    nmod_mat_init(Lag, npts, npts, mod.n);
    nmod_mat_init(iV, npts, npts, mod.n);
    nmod_mat_init(Ta, ellA, m * n, mod.n);
    nmod_mat_init(Tb, ellB, n * p, mod.n);
    nmod_mat_init(Aval, npts, m * n, mod.n);
    nmod_mat_init(Bval, npts, n * p, mod.n);
    nmod_mat_init(Cval, npts, m * p, mod.n);
    nmod_mat_init(Cc, npts, m * p, mod.n);

    for (t = 0; t < npts; t++)
        pts[t] = nmod_set_ui(t, mod);
    _vdm_powers(V, pts, mod);

    /* evaluations of a and b; both are read before c is written,
       so that c may alias them */
    _vdm_pack(Ta, a, 0, 1);
    _vdm_mul_window(Aval, V, Ta, ellA);
    _vdm_pack(Tb, b, 0, 1);
    _vdm_mul_window(Bval, V, Tb, ellB);

    _vdm_pointwise(Cval, Aval, Bval, m, n, p);

    /* interpolation */
    _vdm_lagrange(Lag, pts, mod);
    nmod_mat_transpose(iV, Lag);
    _vdm_mul_window(Cc, iV, Cval, npts);

    _vdm_unpack_init(c, m, p, npts);
    _vdm_unpack(c, m, p, Cc, 0, 1);
    _vdm_unpack_normalise(c, m, p);

    nmod_mat_clear(V);
    nmod_mat_clear(Lag);
    nmod_mat_clear(iV);
    nmod_mat_clear(Ta);
    nmod_mat_clear(Tb);
    nmod_mat_clear(Aval);
    nmod_mat_clear(Bval);
    nmod_mat_clear(Cval);
    nmod_mat_clear(Cc);
    _nmod_vec_clear(pts);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* same as above, with the points +-1, +-2, ..., +-npts       */
/*------------------------------------------------------------*/
void nmod_poly_mat_mul_vandermonde2(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b)
{
    const slong m = a->r, n = a->c, p = b->c;
    const slong ellA = nmod_poly_mat_max_length(a);
    const slong ellB = nmod_poly_mat_max_length(b);
    slong npts;
    ulong inv2;
    nmod_t mod;
    nn_ptr x, x2, cst;
    nmod_mat_t V, Lag, iV, Aplus, Aminus, Bplus, Bminus, Cplus, Cminus, Ce, Co;

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }

    nmod_init(&mod, a->modulus);
    /* 2*npts >= ellA + ellB - 1 */
    npts = (ellA + ellB - 2) / 2 + 1;

    x = _nmod_vec_init(npts);
    x2 = _nmod_vec_init(npts);
    cst = _nmod_vec_init(npts);
    nmod_mat_init(V, npts, (FLINT_MAX(ellA, ellB) + 1) / 2, mod.n);
    nmod_mat_init(Lag, npts, npts, mod.n);
    nmod_mat_init(iV, npts, npts, mod.n);
    nmod_mat_init(Aplus, npts, m * n, mod.n);
    nmod_mat_init(Aminus, npts, m * n, mod.n);
    nmod_mat_init(Bplus, npts, n * p, mod.n);
    nmod_mat_init(Bminus, npts, n * p, mod.n);
    nmod_mat_init(Cplus, npts, m * p, mod.n);
    nmod_mat_init(Cminus, npts, m * p, mod.n);
    nmod_mat_init(Ce, npts, m * p, mod.n);
    nmod_mat_init(Co, npts, m * p, mod.n);

    _vdm2_points(x, x2, npts, mod);
    _vdm_powers(V, x2, mod);

    /* evaluations of a and b at +-x_t; both are read before c is
       written, so that c may alias them */
    _vdm2_eval(Aplus, Aminus, V, x, a, (ellA + 1) / 2, 0, ellA / 2, 1, 2, mod);
    _vdm2_eval(Bplus, Bminus, V, x, b, (ellB + 1) / 2, 0, ellB / 2, 1, 2, mod);

    _vdm_pointwise(Cplus, Aplus, Bplus, m, n, p);
    _vdm_pointwise(Cminus, Aminus, Bminus, m, n, p);

    /* values of the even part (C(x) + C(-x)) / 2 and of the odd part
       (C(x) - C(-x)) / (2x) at x_t^2 */
    _vdm_butterfly(Cplus, Cminus, mod);
    inv2 = nmod_inv(2, mod);
    nmod_mat_scalar_mul(Cplus, Cplus, inv2);
    _vdm2_inv_2x(cst, x, npts, mod);
    _vdm_scale_rows(Cminus, cst, mod);

    /* interpolation */
    _vdm_lagrange(Lag, x2, mod);
    nmod_mat_transpose(iV, Lag);
    _vdm_mul_window(Ce, iV, Cplus, npts);
    _vdm_mul_window(Co, iV, Cminus, npts);

    _vdm_unpack_init(c, m, p, 2 * npts);
    _vdm_unpack(c, m, p, Ce, 0, 2);
    _vdm_unpack(c, m, p, Co, 1, 2);
    _vdm_unpack_normalise(c, m, p);

    nmod_mat_clear(V);
    nmod_mat_clear(Lag);
    nmod_mat_clear(iV);
    nmod_mat_clear(Aplus);
    nmod_mat_clear(Aminus);
    nmod_mat_clear(Bplus);
    nmod_mat_clear(Bminus);
    nmod_mat_clear(Cplus);
    nmod_mat_clear(Cminus);
    nmod_mat_clear(Ce);
    nmod_mat_clear(Co);
    _nmod_vec_clear(x);
    _nmod_vec_clear(x2);
    _nmod_vec_clear(cst);
}

/*------------------------------------------------------------*/
/* middle products                                            */
/*------------------------------------------------------------*/

/*
    Common preamble of the middle products: clamps nhi, handles the empty
    cases (returns 0 after zeroing res), and chooses the operand S which
    is reversed (the shorter one) and the one L which is windowed. Sets
    *ls (length of S), *off (z_s = L_{s + off}), *s_left (S is pmat1).
*/
static int
_vdm_mulmid_setup(nmod_poly_mat_t res,
                  const nmod_poly_mat_t pmat1, slong len1,
                  const nmod_poly_mat_t pmat2, slong len2,
                  slong nlo, slong * nhi, slong * ls, slong * off, int * s_left)
{
    if (pmat1->r == 0 || pmat2->c == 0)
        return 0;

    if (len1 <= 0 || len2 <= 0 || pmat1->c == 0)
    {
        nmod_poly_mat_zero(res);
        return 0;
    }

    *nhi = FLINT_MIN(*nhi, len1 + len2 - 1);
    if (nlo >= *nhi)
    {
        nmod_poly_mat_zero(res);
        return 0;
    }

    *s_left = (len1 <= len2);
    *ls = *s_left ? len1 : len2;
    *off = nlo - *ls + 1;
    return 1;
}

void _nmod_poly_mat_mulmid_vandermonde1(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1, slong len1,
                                        const nmod_poly_mat_t pmat2, slong len2,
                                        slong nlo, slong nhi)
{
    const slong m = pmat1->r, n = pmat1->c, p = pmat2->c;
    slong ls, off, N, npts, t;
    int s_left;
    nmod_t mod;
    nn_ptr pts;
    nmod_mat_t V, VT, Lag, Ts, Tl, Sval, Zhat, Cval, Cc;

    if (!_vdm_mulmid_setup(res, pmat1, len1, pmat2, len2, nlo, &nhi, &ls, &off, &s_left))
        return;

    const nmod_poly_mat_struct * S = s_left ? pmat1 : pmat2;
    const nmod_poly_mat_struct * L = s_left ? pmat2 : pmat1;

    nmod_init(&mod, pmat1->modulus);
    N = nhi - nlo;
    npts = ls + N - 1;

    pts = _nmod_vec_init(npts);
    nmod_mat_init(V, npts, FLINT_MAX(ls, N), mod.n);
    nmod_mat_init(VT, V->c, npts, mod.n);
    nmod_mat_init(Lag, npts, npts, mod.n);
    nmod_mat_init(Ts, ls, S->r * S->c, mod.n);
    nmod_mat_init(Tl, npts, L->r * L->c, mod.n);
    nmod_mat_init(Sval, npts, S->r * S->c, mod.n);
    nmod_mat_init(Zhat, npts, L->r * L->c, mod.n);
    nmod_mat_init(Cval, npts, m * p, mod.n);
    nmod_mat_init(Cc, N, m * p, mod.n);

    for (t = 0; t < npts; t++)
        pts[t] = nmod_set_ui(t, mod);
    _vdm_powers(V, pts, mod);
    _vdm_lagrange(Lag, pts, mod);

    /* evaluations of rev(S) */
    _vdm_pack(Ts, S, ls - 1, -1);
    _vdm_mul_window(Sval, V, Ts, ls);

    /* iV^T applied to the window of L; both operands are read before
       res is written, so that res may alias them */
    _vdm_pack(Tl, L, off, 1);
    _vdm_mul_window(Zhat, Lag, Tl, npts);

    if (s_left)
        _vdm_pointwise(Cval, Sval, Zhat, m, n, p);
    else
        _vdm_pointwise(Cval, Zhat, Sval, m, n, p);

    /* V_N^T: the first N rows of V^T */
    nmod_mat_transpose(VT, V);
    _vdm_mul_window(Cc, VT, Cval, npts);

    _vdm_unpack_init(res, m, p, N);
    _vdm_unpack(res, m, p, Cc, 0, 1);
    _vdm_unpack_normalise(res, m, p);

    nmod_mat_clear(V);
    nmod_mat_clear(VT);
    nmod_mat_clear(Lag);
    nmod_mat_clear(Ts);
    nmod_mat_clear(Tl);
    nmod_mat_clear(Sval);
    nmod_mat_clear(Zhat);
    nmod_mat_clear(Cval);
    nmod_mat_clear(Cc);
    _nmod_vec_clear(pts);
}

void _nmod_poly_mat_mulmid_vandermonde2(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1, slong len1,
                                        const nmod_poly_mat_t pmat2, slong len2,
                                        slong nlo, slong nhi)
{
    const slong m = pmat1->r, n = pmat1->c, p = pmat2->c;
    slong ls, off, N, npts;
    int s_left;
    ulong inv2;
    nmod_t mod;
    nn_ptr x, x2, cst;
    nmod_mat_t V, VT, Lag, Te, To, Splus, Sminus, Ze, Zo, Dplus, Dminus, Ye, Yo;

    if (!_vdm_mulmid_setup(res, pmat1, len1, pmat2, len2, nlo, &nhi, &ls, &off, &s_left))
        return;

    const nmod_poly_mat_struct * S = s_left ? pmat1 : pmat2;
    const nmod_poly_mat_struct * L = s_left ? pmat2 : pmat1;

    nmod_init(&mod, pmat1->modulus);
    N = nhi - nlo;
    /* the transposed product has length ls + N - 1 <= 2*npts */
    npts = (ls + N - 2) / 2 + 1;

    x = _nmod_vec_init(npts);
    x2 = _nmod_vec_init(npts);
    cst = _nmod_vec_init(npts);
    nmod_mat_init(V, npts, FLINT_MAX((ls + 1) / 2, (N + 1) / 2), mod.n);
    nmod_mat_init(VT, V->c, npts, mod.n);
    nmod_mat_init(Lag, npts, npts, mod.n);
    nmod_mat_init(Te, npts, L->r * L->c, mod.n);
    nmod_mat_init(To, npts, L->r * L->c, mod.n);
    nmod_mat_init(Splus, npts, S->r * S->c, mod.n);
    nmod_mat_init(Sminus, npts, S->r * S->c, mod.n);
    nmod_mat_init(Ze, npts, L->r * L->c, mod.n);
    nmod_mat_init(Zo, npts, L->r * L->c, mod.n);
    nmod_mat_init(Dplus, npts, m * p, mod.n);
    nmod_mat_init(Dminus, npts, m * p, mod.n);
    nmod_mat_init(Ye, (N + 1) / 2, m * p, mod.n);
    nmod_mat_init(Yo, N / 2, m * p, mod.n);

    _vdm2_points(x, x2, npts, mod);
    _vdm_powers(V, x2, mod);
    _vdm_lagrange(Lag, x2, mod);

    /* rev(S) at +-x_t: its even part is the coefficients ls-1, ls-3, ...
       of S, its odd part ls-2, ls-4, ... */
    _vdm2_eval(Splus, Sminus, V, x, S, (ls + 1) / 2, ls - 1, ls / 2, ls - 2, -2, mod);

    /* the transpose of the interpolation, on the even and odd parts of the
       window z of length 2*npts; both operands are read before res is
       written, so that res may alias them */
    _vdm_pack(Te, L, off, 2);
    _vdm_mul_window(Ze, Lag, Te, npts);
    _vdm_pack(To, L, off + 1, 2);
    _vdm_mul_window(Zo, Lag, To, npts);

    /* transpose of (C(x), C(-x)) -> ((C(x)+C(-x))/2, (C(x)-C(-x))/(2x)):
       (Ze, Zo) -> (Ze/2 + Zo/(2x), Ze/2 - Zo/(2x)) */
    inv2 = nmod_inv(2, mod);
    nmod_mat_scalar_mul(Ze, Ze, inv2);
    _vdm2_inv_2x(cst, x, npts, mod);
    _vdm_scale_rows(Zo, cst, mod);
    _vdm_butterfly(Ze, Zo, mod);

    /* pointwise products by rev(S)(x) and rev(S)(-x) */
    if (s_left)
    {
        _vdm_pointwise(Dplus, Splus, Ze, m, n, p);
        _vdm_pointwise(Dminus, Sminus, Zo, m, n, p);
    }
    else
    {
        _vdm_pointwise(Dplus, Ze, Splus, m, n, p);
        _vdm_pointwise(Dminus, Zo, Sminus, m, n, p);
    }

    /* transpose of the evaluation of y at +-x_t from its even and odd
       parts: (Dplus, Dminus) -> (Dplus + Dminus, x (Dplus - Dminus)),
       then the first (N+1)/2, resp. N/2, rows of V^T on each */
    _vdm_butterfly(Dplus, Dminus, mod);
    _vdm_scale_rows(Dminus, x, mod);

    nmod_mat_transpose(VT, V);
    _vdm_mul_window(Ye, VT, Dplus, npts);
    _vdm_mul_window(Yo, VT, Dminus, npts);

    _vdm_unpack_init(res, m, p, N);
    _vdm_unpack(res, m, p, Ye, 0, 2);
    _vdm_unpack(res, m, p, Yo, 1, 2);
    _vdm_unpack_normalise(res, m, p);

    nmod_mat_clear(V);
    nmod_mat_clear(VT);
    nmod_mat_clear(Lag);
    nmod_mat_clear(Te);
    nmod_mat_clear(To);
    nmod_mat_clear(Splus);
    nmod_mat_clear(Sminus);
    nmod_mat_clear(Ze);
    nmod_mat_clear(Zo);
    nmod_mat_clear(Dplus);
    nmod_mat_clear(Dminus);
    nmod_mat_clear(Ye);
    nmod_mat_clear(Yo);
    _nmod_vec_clear(x);
    _nmod_vec_clear(x2);
    _nmod_vec_clear(cst);
}
