#include <string.h>  // for memcpy, memset

#include "machine_vectors.h"
#include "nmod_mat_poly.h"
#include "nmod_mat_poly_extra/impl.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

void nmod_mat_poly_init_set_from_nmod_mat(nmod_mat_poly_t matp,
                                          const nmod_mat_t cmat)
{
    nmod_mat_poly_init2_preinv(matp, cmat->r, cmat->c, cmat->mod.n, cmat->mod.ninv, 1);
    // if cmat is zero: do nothing more
    // otherwise set length to 1 and copy data
    if (! nmod_mat_is_zero(cmat))
    {
        nmod_mat_t coeff;
        _nmod_mat_poly_set_length(matp, 1);
        nmod_mat_poly_coeff_attach(coeff, matp, 0);
        nmod_mat_set(coeff, cmat);
    }
}

void nmod_mat_poly_set_from_nmod_mat(nmod_mat_poly_t matp,
                                     const nmod_mat_t cmat)
{
    if (nmod_mat_is_zero(cmat))
        nmod_mat_poly_zero(matp);
    else
    {
        nmod_mat_t coeff;
        nmod_mat_poly_fit_length(matp, 1);
        _nmod_mat_poly_set_length(matp, 1);
        nmod_mat_poly_coeff_attach(coeff, matp, 0);
        nmod_mat_set(coeff, cmat);
    }
}

/*------------------------------------------------------------------------*/
/* Transposition between two collections of separately allocated rows      */
/*                                                                         */
/* Both conversions between nmod_poly_mat and nmod_mat_poly are the same    */
/* transposition, with the two sides exchanged.  Writing `n = r*c` for the  */
/* number of matrix entries and `len` for the number of coefficients, an    */
/* nmod_poly_mat is an `n x len` array whose rows (the coefficient arrays   */
/* of the polynomial entries) are contiguous, and an nmod_mat_poly is a     */
/* `len x n` array whose rows (the entry arrays of the matrix coefficients) */
/* are contiguous.  ::_pml_transpose does one such transposition and is     */
/* used by both directions; only which table is passed as the source and    */
/* which as the destination differs.                                       */
/*                                                                         */
/* The naive element-by-element loop uses 8 bytes of each cache line it     */
/* touches on one of the two sides before moving on, so as soon as that     */
/* side stops fitting in cache it pays a miss per *word* instead of per     */
/* cache line.  The routines below move a `WxW` block at a time, `W = 4` or */
/* `8`, so that a whole cache line, or at least half of one, is filled at   */
/* every visit.  Measured against the naive loop, that is worth 5x to 9x on */
/* Zen 4 and 2x to 5x on Apple M4 over matrices 4x4 to 128x128 and lengths  */
/* 32 to 2048.                                                             */
/*                                                                         */
/* Alignment decides how wide a kernel is worth using, and the two          */
/* directions differ there: the coefficients of an nmod_mat_poly are        */
/* 64-byte aligned (see ::_nmod_mat_poly_coeff_alloc) so a 64-byte access   */
/* lands on one cache line, while the coefficient array of an nmod_poly is  */
/* whatever `flint_realloc` returns, so the same access straddles two.      */
/* Hence the two conversions do not pick the same default kernel.          */
/*------------------------------------------------------------------------*/

/* Load W words of `s` starting at index `k`, zero-padding past `l`.
   `l` is the (already truncated) length of the source polynomial: reading
   at or beyond it is not allowed, the coefficients there are not stored. */
#define _PML_LOADPAD(buf, s, l, k, W)                                        \
    do {                                                                     \
        const slong _av = (l) - (k);                                         \
        if (_av >= (W))                                                      \
            memcpy((buf), (s) + (k), (W) * sizeof(ulong));                   \
        else                                                                 \
        {                                                                    \
            memset((buf), 0, (W) * sizeof(ulong));                           \
            if (_av > 0)                                                     \
                memcpy((buf), (s) + (k), (size_t) _av * sizeof(ulong));      \
        }                                                                    \
    } while (0)

/* --- 8x8 scalar kernel (portable fallback) --- */
#define _PML_KERNEL_SCALAR(e, k)                                             \
    do {                                                                     \
        ulong _b[8][8];                                                      \
        for (slong _t = 0; _t < 8; _t++)                                     \
            _PML_LOADPAD(_b[_t], src[(e) + _t], slen[(e) + _t], (k), 8);     \
        for (slong _u = 0; _u < 8; _u++)                                     \
        {                                                                    \
            nn_ptr _d = dst[(k) + _u] + (e);                                 \
            for (slong _t = 0; _t < 8; _t++)                                 \
                _d[_t] = _b[_t][_u];                                         \
        }                                                                    \
    } while (0)

/* --- 4x4 kernel on ulong machine vectors: AVX2 and NEON --- */
#if PML_HAVE_MACHINE_VECTORS
# define _PML_HAVE_CONV_VEC4 1

FLINT_FORCE_INLINE vec4n _pml_conv_load4(nn_srcptr s, slong l, slong k)
{
    if (k + 4 <= l)
        return vec4n_load_unaligned(s + k);
    ulong b[4];
    _PML_LOADPAD(b, s, l, k, 4);
    return vec4n_load_unaligned(b);
}

#define _PML_KERNEL_VEC4(e, k)                                               \
    do {                                                                     \
        vec4n _a0 = _pml_conv_load4(src[(e) + 0], slen[(e) + 0], (k));       \
        vec4n _a1 = _pml_conv_load4(src[(e) + 1], slen[(e) + 1], (k));       \
        vec4n _a2 = _pml_conv_load4(src[(e) + 2], slen[(e) + 2], (k));       \
        vec4n _a3 = _pml_conv_load4(src[(e) + 3], slen[(e) + 3], (k));       \
        vec4n _z0, _z1, _z2, _z3;                                            \
        VEC4N_TRANSPOSE(_z0, _z1, _z2, _z3, _a0, _a1, _a2, _a3);             \
        vec4n_store_unaligned(dst[(k) + 0] + (e), _z0);                      \
        vec4n_store_unaligned(dst[(k) + 1] + (e), _z1);                      \
        vec4n_store_unaligned(dst[(k) + 2] + (e), _z2);                      \
        vec4n_store_unaligned(dst[(k) + 3] + (e), _z3);                      \
    } while (0)
#endif  /* _PML_HAVE_CONV_VEC4 */

/* --- 8x8 kernel, AVX-512 --- */
#if PML_HAVE_AVX512
# define _PML_HAVE_CONV_VEC8 1

FLINT_FORCE_INLINE __m512i _pml_conv_load8(nn_srcptr s, slong l, slong k)
{
    if (k + 8 <= l)
        return _mm512_loadu_si512((const void *) (s + k));
    ulong b[8];
    _PML_LOADPAD(b, s, l, k, 8);
    return _mm512_loadu_si512((const void *) b);
}

#define _PML_KERNEL_VEC8(e, k)                                               \
    do {                                                                     \
        __m512i _a0 = _pml_conv_load8(src[(e) + 0], slen[(e) + 0], (k));     \
        __m512i _a1 = _pml_conv_load8(src[(e) + 1], slen[(e) + 1], (k));     \
        __m512i _a2 = _pml_conv_load8(src[(e) + 2], slen[(e) + 2], (k));     \
        __m512i _a3 = _pml_conv_load8(src[(e) + 3], slen[(e) + 3], (k));     \
        __m512i _a4 = _pml_conv_load8(src[(e) + 4], slen[(e) + 4], (k));     \
        __m512i _a5 = _pml_conv_load8(src[(e) + 5], slen[(e) + 5], (k));     \
        __m512i _a6 = _pml_conv_load8(src[(e) + 6], slen[(e) + 6], (k));     \
        __m512i _a7 = _pml_conv_load8(src[(e) + 7], slen[(e) + 7], (k));     \
        __m512i _z0, _z1, _z2, _z3, _z4, _z5, _z6, _z7;                      \
        VEC8N_TRANSPOSE(_z0, _z1, _z2, _z3, _z4, _z5, _z6, _z7,              \
                        _a0, _a1, _a2, _a3, _a4, _a5, _a6, _a7);             \
        _mm512_storeu_si512((void *) (dst[(k) + 0] + (e)), _z0);             \
        _mm512_storeu_si512((void *) (dst[(k) + 1] + (e)), _z1);             \
        _mm512_storeu_si512((void *) (dst[(k) + 2] + (e)), _z2);             \
        _mm512_storeu_si512((void *) (dst[(k) + 3] + (e)), _z3);             \
        _mm512_storeu_si512((void *) (dst[(k) + 4] + (e)), _z4);             \
        _mm512_storeu_si512((void *) (dst[(k) + 5] + (e)), _z5);             \
        _mm512_storeu_si512((void *) (dst[(k) + 6] + (e)), _z6);             \
        _mm512_storeu_si512((void *) (dst[(k) + 7] + (e)), _z7);             \
    } while (0)
#endif  /* _PML_HAVE_CONV_VEC8 */

/* One transposition routine, for a given block width, block kernel and
   schedule.

   DST_MAJOR == 1: the destination row index is the outer loop, so the stores
   walk each destination row from left to right (long sequential write
   streams, scattered reads).  DST_MAJOR == 0 is the transposed schedule
   (long sequential read streams, scattered writes).

   Either way the two residual bands, at most W-1 wide each, are finished with
   plain scalar loops. */
#define _PML_MK_TRANSPOSE(NAME, W, KERNEL, DST_MAJOR)                           \
static void NAME(nn_ptr * dst, slong len,                                    \
                 nn_srcptr * src, const slong * slen, slong n)               \
{                                                                            \
    const slong nW = n - (n % (W));                                          \
    const slong lW = len - (len % (W));                                      \
                                                                             \
    if (DST_MAJOR)                                                         \
        for (slong k = 0; k < lW; k += (W))                                  \
            for (slong e = 0; e < nW; e += (W))                              \
                KERNEL(e, k);                                                \
    else                                                                     \
        for (slong e = 0; e < nW; e += (W))                                  \
            for (slong k = 0; k < lW; k += (W))                              \
                KERNEL(e, k);                                                \
                                                                             \
    for (slong k = lW; k < len; k++)                                         \
    {                                                                        \
        nn_ptr d = dst[k];                                                   \
        for (slong e = 0; e < n; e++)                                        \
            d[e] = (k < slen[e]) ? src[e][k] : UWORD(0);                     \
    }                                                                        \
    for (slong e = nW; e < n; e++)                                           \
    {                                                                        \
        nn_srcptr s = src[e];                                                \
        const slong l = FLINT_MIN(slen[e], lW);                              \
        for (slong k = 0; k < l; k++)                                        \
            dst[k][e] = s[k];                                                \
        for (slong k = l; k < lW; k++)                                       \
            dst[k][e] = UWORD(0);                                            \
    }                                                                        \
}

_PML_MK_TRANSPOSE(_pml_tr_sca_dm, 8, _PML_KERNEL_SCALAR, 1)
_PML_MK_TRANSPOSE(_pml_tr_sca_sm, 8, _PML_KERNEL_SCALAR, 0)
#if _PML_HAVE_CONV_VEC4
_PML_MK_TRANSPOSE(_pml_tr_v4_dm, 4, _PML_KERNEL_VEC4, 1)
_PML_MK_TRANSPOSE(_pml_tr_v4_sm, 4, _PML_KERNEL_VEC4, 0)
#endif
#if _PML_HAVE_CONV_VEC8
_PML_MK_TRANSPOSE(_pml_tr_v8_dm, 8, _PML_KERNEL_VEC8, 1)
_PML_MK_TRANSPOSE(_pml_tr_v8_sm, 8, _PML_KERNEL_VEC8, 0)
#endif

/* Dispatch on (kernel, schedule).  Documented in nmod_mat_poly_extra/impl.h. */
void _pml_transpose(nn_ptr * dst, slong ndst, nn_srcptr * src,
                    const slong * slen, slong nsrc, int kern, int dmaj)
{
#if _PML_HAVE_CONV_VEC8
    if (kern == NMOD_MAT_POLY_CONV_VEC8)
    {
        if (dmaj) _pml_tr_v8_dm(dst, ndst, src, slen, nsrc);
        else      _pml_tr_v8_sm(dst, ndst, src, slen, nsrc);
        return;
    }
#endif
#if _PML_HAVE_CONV_VEC4
    if (kern == NMOD_MAT_POLY_CONV_VEC4)
    {
        if (dmaj) _pml_tr_v4_dm(dst, ndst, src, slen, nsrc);
        else      _pml_tr_v4_sm(dst, ndst, src, slen, nsrc);
        return;
    }
#endif
    if (dmaj) _pml_tr_sca_dm(dst, ndst, src, slen, nsrc);
    else      _pml_tr_sca_sm(dst, ndst, src, slen, nsrc);
}

/* Narrow `kern` to a kernel whose block fits inside a `ndst x nsrc`
   transposition: a wider one would leave all the work to the scalar residual
   bands.  This matters for small matrices -- an 8x8 kernel does nothing at
   all on a 2x2 matrix (4 entries), where the 4x4 one is 4x faster. */
int _pml_transpose_narrow_kernel(int kern, slong ndst, slong nsrc)
{
    if (kern == NMOD_MAT_POLY_CONV_VEC8 && (nsrc < 8 || ndst < 8))
        kern = NMOD_MAT_POLY_CONV_VEC4;
    if (kern == NMOD_MAT_POLY_CONV_VEC4 && (nsrc < 4 || ndst < 4))
        kern = NMOD_MAT_POLY_CONV_SCALAR;
    return kern;
}

/* The widest kernel this build provides. */
int _pml_transpose_widest_kernel(void)
{
#if _PML_HAVE_CONV_VEC8
    return NMOD_MAT_POLY_CONV_VEC8;
#elif _PML_HAVE_CONV_VEC4
    return NMOD_MAT_POLY_CONV_VEC4;
#else
    return NMOD_MAT_POLY_CONV_SCALAR;
#endif
}

/*------------------------------------------------------------------------*/
/* blocked implementation                                                  */
/*------------------------------------------------------------------------*/

/* Below this many words the pointer tables and their allocation dominate;
   just run the naive loop. */
#define _PML_CONV_TINY 512

void _nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                            const nmod_poly_mat_t pmat,
                                            slong order,
                                            int kern,
                                            int dmaj)
{
    const slong len = nmod_poly_mat_max_length(pmat);
    if (order > len)
        order = len;

    // allocate and init coefficients
    nmod_mat_poly_fit_length(matp, order);
    _nmod_mat_poly_set_length(matp, order);

    const slong r = matp->r;
    const slong c = matp->c;

    if (order == 0 || r == 0 || c == 0)
    {
        if (order < len)
            _nmod_mat_poly_normalise(matp);
        return;
    }

    if ((double) r * (double) c * (double) order < (double) _PML_CONV_TINY)
    {
        for (slong k = 0; k < order; k++)
            for (slong i = 0; i < r; i++)
                for (slong j = 0; j < c; j++)
                    nmod_mat_poly_entry(matp, k, i, j) = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k);
        if (order < len)
            _nmod_mat_poly_normalise(matp);
        return;
    }

    /* Effective number of "rows" of the transposition: the whole matrix when
       the coefficients are contiguous, one matrix row otherwise. */
    const slong nrows = (matp->stride == c) ? r * c : c;

    /* Default kernel: the widest one available.  On Zen 4 the 8-wide AVX-512
       kernel is 20% to 40% faster than the 4-wide one as long as the data
       fits in cache (0.20 vs 0.28 ns per word at 32x32 and length 128), and
       10% to 20% slower once it does not (0.95 vs 0.86 at 64x64 and length
       2048); close enough overall to keep the wider one for its advantage in
       the common, cache-resident range.  This relies on the destination rows
       here being 64-byte aligned, as ::_nmod_mat_poly_coeff_alloc makes them;
       the conversion in the other direction, whose destination rows are not,
       makes a different choice. */
    if (kern < 0 || kern > NMOD_MAT_POLY_CONV_VEC8)
        kern = _pml_transpose_widest_kernel();

    kern = _pml_transpose_narrow_kernel(kern, order, nrows);

    /* Schedule; here destination-major is the coefficient-major schedule (the
       stores walk each output matrix coefficient from left to right) and
       source-major is the entry-major one (the loads walk each input
       polynomial).  See nmod_mat_poly_extra/impl.h for why the default
       depends on the target. */
    if (dmaj < 0 || dmaj > 1)
        dmaj = PML_CONV_DST_MAJOR;

    nn_srcptr * src = (nn_srcptr *) flint_malloc(r * c * sizeof(nn_srcptr));
    slong * slen = (slong *) flint_malloc(r * c * sizeof(slong));

    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        {
            const nmod_poly_struct * p = nmod_poly_mat_entry(pmat, i, j);
            src[i * c + j] = p->coeffs;
            slen[i * c + j] = FLINT_MIN(p->length, order);
        }

    if (matp->stride == c)
    {
        /* the usual case: the entries of a coefficient are one contiguous run
           of r*c words, 64-byte aligned, so the whole matrix is a single row
           of the transposition */
        _pml_transpose(matp->coeffs, order, src, slen, r * c, kern, dmaj);
    }
    else
    {
        /* padded or otherwise nonstandard stride: one transposition per row */
        nn_ptr * dst = (nn_ptr *) flint_malloc(order * sizeof(nn_ptr));
        for (slong i = 0; i < r; i++)
        {
            for (slong k = 0; k < order; k++)
                dst[k] = matp->coeffs[k] + i * matp->stride;
            _pml_transpose(dst, order, src + i * c, slen + i * c, c, kern, dmaj);
        }
        flint_free(dst);
    }

    flint_free(src);
    flint_free(slen);

    // normalize (useless if order==len)
    if (order < len)
        _nmod_mat_poly_normalise(matp);
}

void nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                           const nmod_poly_mat_t pmat,
                                           slong order)
{
    _nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, order, -1, -1);
}
