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
/* Conversion nmod_poly_mat -> nmod_mat_poly                               */
/*                                                                         */
/* This is a transposition.  Writing `n = r*c` for the number of matrix     */
/* entries and `len` for the number of coefficients, the input is an        */
/* `n x len` array whose rows (the coefficient arrays of the polynomial     */
/* entries) are contiguous, and the output is a `len x n` array whose rows  */
/* (the entry arrays of the matrix coefficients) are contiguous.  Both      */
/* sides are collections of separately allocated rows.                      */
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
/* The output rows being 64-byte aligned (see ::_nmod_mat_poly_coeff_alloc) */
/* is what makes the 8-wide kernel worth having: a 64-byte store lands on   */
/* one cache line rather than straddling two, which it would do at every    */
/* single store if the rows were only 16-byte aligned, as a plain           */
/* `flint_calloc` leaves them.                                             */
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

/* One conversion routine, for a given block width, block kernel and schedule.

   COEFF_MAJOR == 1: the coefficient index is the outer loop, so the stores
   walk each output matrix from left to right (long sequential write streams,
   scattered reads).  COEFF_MAJOR == 0 is the transposed schedule (long
   sequential read streams, scattered writes).

   Either way the two residual bands, at most W-1 wide each, are finished with
   plain scalar loops. */
#define _PML_MK_CONV(NAME, W, KERNEL, COEFF_MAJOR)                           \
static void NAME(nn_ptr * dst, slong len,                                    \
                 nn_srcptr * src, const slong * slen, slong n)               \
{                                                                            \
    const slong nW = n - (n % (W));                                          \
    const slong lW = len - (len % (W));                                      \
                                                                             \
    if (COEFF_MAJOR)                                                         \
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

_PML_MK_CONV(_pml_conv_sca_cm, 8, _PML_KERNEL_SCALAR, 1)
_PML_MK_CONV(_pml_conv_sca_em, 8, _PML_KERNEL_SCALAR, 0)
#if _PML_HAVE_CONV_VEC4
_PML_MK_CONV(_pml_conv_v4_cm, 4, _PML_KERNEL_VEC4, 1)
_PML_MK_CONV(_pml_conv_v4_em, 4, _PML_KERNEL_VEC4, 0)
#endif
#if _PML_HAVE_CONV_VEC8
_PML_MK_CONV(_pml_conv_v8_cm, 8, _PML_KERNEL_VEC8, 1)
_PML_MK_CONV(_pml_conv_v8_em, 8, _PML_KERNEL_VEC8, 0)
#endif

/* dispatch on (kernel, schedule) */
static void _pml_conv(nn_ptr * dst, slong len, nn_srcptr * src,
                      const slong * slen, slong n, int kern, int cmaj)
{
#if _PML_HAVE_CONV_VEC8
    if (kern == NMOD_MAT_POLY_CONV_VEC8)
    {
        if (cmaj) _pml_conv_v8_cm(dst, len, src, slen, n);
        else      _pml_conv_v8_em(dst, len, src, slen, n);
        return;
    }
#endif
#if _PML_HAVE_CONV_VEC4
    if (kern == NMOD_MAT_POLY_CONV_VEC4)
    {
        if (cmaj) _pml_conv_v4_cm(dst, len, src, slen, n);
        else      _pml_conv_v4_em(dst, len, src, slen, n);
        return;
    }
#endif
    if (cmaj) _pml_conv_sca_cm(dst, len, src, slen, n);
    else      _pml_conv_sca_em(dst, len, src, slen, n);
}

/* Default kernel: the widest one the build provides.
 *
 * On Zen 4 the 8-wide AVX-512 kernel is 20% to 40% faster than the 4-wide one
 * as long as the data fits in cache (0.20 vs 0.28 ns per word at 32x32 and
 * length 128), and 10% to 20% slower once it does not (0.95 vs 0.86 at 64x64
 * and length 2048).  The two are close enough overall that the wider one is
 * kept for its advantage in the common, cache-resident range.
 *
 * The ranking depends on the output rows being 64-byte aligned, as
 * ::_nmod_mat_poly_coeff_alloc makes them: with the 16-byte alignment a plain
 * `flint_calloc` leaves, every 64-byte store straddles two cache lines and
 * the 4-wide kernel wins instead. */
static int _pml_conv_default_kernel(void)
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
                                            int cmaj)
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

    if (kern < 0 || kern > NMOD_MAT_POLY_CONV_VEC8)
        kern = _pml_conv_default_kernel();

    /* A kernel whose block does not fit inside the problem would leave all
       the work to the scalar residual bands, so step down to a narrower one.
       This matters for small matrices: an 8x8 kernel does nothing at all on a
       2x2 matrix (4 entries), where the 4x4 one is 4x faster. */
    if (kern == NMOD_MAT_POLY_CONV_VEC8 && (nrows < 8 || order < 8))
        kern = NMOD_MAT_POLY_CONV_VEC4;
    if (kern == NMOD_MAT_POLY_CONV_VEC4 && (nrows < 4 || order < 4))
        kern = NMOD_MAT_POLY_CONV_SCALAR;

    /* Schedule.  Whichever loop is outermost, one side of the transposition
       is visited in long sequential runs and the other one block at a time,
       from rows that are separately allocated.  Entry-major puts the long
       runs on the loads from the input polynomials and scatters the stores
       into the output matrices; coefficient-major does the opposite.

       Which one wins depends on how much of a cache line one scattered store
       covers.  On x86, where a line is 64 bytes, a block is a whole line
       (8-wide kernel) or half of one (4-wide), the store needs no line to be
       kept around, and entry-major wins: on Zen 4 it is 1.6x to 2.4x faster
       than coefficient-major from 16 MB of data upwards, and within noise
       below.

       On Apple silicon a line is 128 bytes, so a scattered store covers a
       quarter of one (the 4-wide NEON kernel is the widest available there):
       the line has to be fetched for ownership and then kept until the
       remaining quarters are written, which only happens after a full sweep
       of the inner loop.  Entry-major then loses badly -- on an M4 it is 3x
       to 6x slower than coefficient-major as soon as the matrix reaches
       16x16, and slower than the naive loop itself from 32x32 on.  With the
       stores sequential instead, consecutive blocks fill each line back to
       back and the line size stops mattering. */
    if (cmaj < 0 || cmaj > 1)
        cmaj = PML_CONV_COEFF_MAJOR;

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
        _pml_conv(matp->coeffs, order, src, slen, r * c, kern, cmaj);
    }
    else
    {
        /* padded or otherwise nonstandard stride: one transposition per row */
        nn_ptr * dst = (nn_ptr *) flint_malloc(order * sizeof(nn_ptr));
        for (slong i = 0; i < r; i++)
        {
            for (slong k = 0; k < order; k++)
                dst[k] = matp->coeffs[k] + i * matp->stride;
            _pml_conv(dst, order, src + i * c, slen + i * c, c, kern, cmaj);
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
