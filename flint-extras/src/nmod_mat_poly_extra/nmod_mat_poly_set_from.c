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
/* The naive element-by-element loop uses 8 bytes of each 64-byte line it   */
/* touches on one of the two sides before moving on, so as soon as that     */
/* side stops fitting in cache it pays a miss per *word* instead of per     */
/* cache line: an 8x read amplification.  The routines below move a `WxW`   */
/* block at a time, `W = 4` or `8`, so that every cache line touched is     */
/* read, or written, in full.                                              */
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

/* Software prefetching.
 *
 * Whichever schedule is used, one of the two sides is visited one cache line
 * per row, from rows that are separately allocated: an access pattern the
 * hardware prefetchers cannot follow, so every one of those accesses is a
 * demand miss (plus, beyond a few thousand rows, a page walk).  The addresses
 * are all known in advance -- they are in the pointer tables -- so they can be
 * prefetched a few blocks ahead.  Measured on a Skylake-SP this is worth 8% to
 * 30% from about 1 MiB of data upwards, and is a small loss below, hence the
 * size test in _nmod_mat_poly_set_trunc_from_poly_mat. */
#if defined(__GNUC__) || defined(__clang__)
# define _PML_PREFETCH_R(p) __builtin_prefetch((const void *) (p), 0, 3)
# define _PML_PREFETCH_W(p) __builtin_prefetch((const void *) (p), 1, 3)
#else
# define _PML_PREFETCH_R(p) ((void) 0)
# define _PML_PREFETCH_W(p) ((void) 0)
#endif

/* how far ahead, counted in rows of the scattered side */
#define _PML_CONV_PFD 32

/* Body shared by all instantiations.

   COEFF_MAJOR == 1: the coefficient index is the outer loop, so the stores
   walk each output matrix from left to right (long sequential write streams,
   scattered reads).  COEFF_MAJOR == 0 is the transposed schedule (long
   sequential read streams, scattered writes).

   Whichever schedule is used, the residual bands are finished with plain
   scalar loops; they are at most W-1 wide.  */
#define _PML_CONV_BODY(W, KERNEL, COEFF_MAJOR, PF)                           \
    const slong nW = n - (n % (W));                                          \
    const slong lW = len - (len % (W));                                      \
                                                                             \
    if (COEFF_MAJOR)                                                         \
    {                                                                        \
        for (slong k = 0; k < lW; k += (W))                                  \
            for (slong e = 0; e < nW; e += (W))                              \
            {                                                                \
                if (PF)                                                      \
                {                                                            \
                    const slong ep = (e + _PML_CONV_PFD < nW)                \
                                   ? e + _PML_CONV_PFD : e;                  \
                    for (slong _t = 0; _t < (W); _t++)                       \
                        _PML_PREFETCH_R(src[ep + _t] + k);                   \
                }                                                            \
                KERNEL(e, k);                                                \
            }                                                                \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        for (slong e = 0; e < nW; e += (W))                                  \
            for (slong k = 0; k < lW; k += (W))                              \
            {                                                                \
                if (PF)                                                      \
                {                                                            \
                    const slong kp = (k + _PML_CONV_PFD < lW)                \
                                   ? k + _PML_CONV_PFD : k;                  \
                    for (slong _u = 0; _u < (W); _u++)                       \
                        _PML_PREFETCH_W(dst[kp + _u] + e);                   \
                }                                                            \
                KERNEL(e, k);                                                \
            }                                                                \
    }                                                                        \
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
    }

#define _PML_MK_CONV(NAME, W, KERNEL, COEFF_MAJOR, PF)                       \
static void NAME(nn_ptr * dst, slong len,                                    \
                 nn_srcptr * src, const slong * slen, slong n)               \
{                                                                            \
    _PML_CONV_BODY(W, KERNEL, COEFF_MAJOR, PF)                               \
}

_PML_MK_CONV(_pml_conv_sca_cm, 8, _PML_KERNEL_SCALAR, 1, 0)
_PML_MK_CONV(_pml_conv_sca_em, 8, _PML_KERNEL_SCALAR, 0, 0)
_PML_MK_CONV(_pml_conv_sca_cm_pf, 8, _PML_KERNEL_SCALAR, 1, 1)
_PML_MK_CONV(_pml_conv_sca_em_pf, 8, _PML_KERNEL_SCALAR, 0, 1)
#if _PML_HAVE_CONV_VEC4
_PML_MK_CONV(_pml_conv_v4_cm, 4, _PML_KERNEL_VEC4, 1, 0)
_PML_MK_CONV(_pml_conv_v4_em, 4, _PML_KERNEL_VEC4, 0, 0)
_PML_MK_CONV(_pml_conv_v4_cm_pf, 4, _PML_KERNEL_VEC4, 1, 1)
_PML_MK_CONV(_pml_conv_v4_em_pf, 4, _PML_KERNEL_VEC4, 0, 1)
#endif
#if _PML_HAVE_CONV_VEC8
_PML_MK_CONV(_pml_conv_v8_cm, 8, _PML_KERNEL_VEC8, 1, 0)
_PML_MK_CONV(_pml_conv_v8_em, 8, _PML_KERNEL_VEC8, 0, 0)
_PML_MK_CONV(_pml_conv_v8_cm_pf, 8, _PML_KERNEL_VEC8, 1, 1)
_PML_MK_CONV(_pml_conv_v8_em_pf, 8, _PML_KERNEL_VEC8, 0, 1)
#endif

/* dispatch on (kernel, schedule, prefetch) */
static void _pml_conv(nn_ptr * dst, slong len, nn_srcptr * src,
                      const slong * slen, slong n, int kern, int cmaj, int pf)
{
#if _PML_HAVE_CONV_VEC8
    if (kern == NMOD_MAT_POLY_CONV_VEC8)
    {
        if (cmaj) { if (pf) _pml_conv_v8_cm_pf(dst, len, src, slen, n);
                    else    _pml_conv_v8_cm(dst, len, src, slen, n); }
        else      { if (pf) _pml_conv_v8_em_pf(dst, len, src, slen, n);
                    else    _pml_conv_v8_em(dst, len, src, slen, n); }
        return;
    }
#endif
#if _PML_HAVE_CONV_VEC4
    if (kern == NMOD_MAT_POLY_CONV_VEC4)
    {
        if (cmaj) { if (pf) _pml_conv_v4_cm_pf(dst, len, src, slen, n);
                    else    _pml_conv_v4_cm(dst, len, src, slen, n); }
        else      { if (pf) _pml_conv_v4_em_pf(dst, len, src, slen, n);
                    else    _pml_conv_v4_em(dst, len, src, slen, n); }
        return;
    }
#endif
    if (cmaj) { if (pf) _pml_conv_sca_cm_pf(dst, len, src, slen, n);
                else    _pml_conv_sca_cm(dst, len, src, slen, n); }
    else      { if (pf) _pml_conv_sca_em_pf(dst, len, src, slen, n);
                else    _pml_conv_sca_em(dst, len, src, slen, n); }
}

/* Default kernel: the widest one the build provides.
 *
 * Measured on a Skylake-SP, with the output rows 64-byte aligned as
 * nmod_mat_poly now allocates them, the 8-wide AVX-512 kernel is 10% to 40%
 * faster than the 4-wide one on every shape tested -- and that is a lower
 * bound for the machines PML_HAVE_AVX512 selects, since those have IFMA, hence
 * are Ice Lake or later or Zen 4 or later, and do not pay the 512-bit
 * frequency licence that Skylake-SP pays here.  (Before the coefficients were
 * 64-byte aligned the ranking was the other way round, by a similar margin:
 * every 64-byte store then straddled two cache lines.) */
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

/* Above this many bytes of input, prefetch the scattered side (see above).
   Measured crossover on a Skylake-SP is between 256 KiB and 1 MiB. */
#ifndef PML_CONV_PREFETCH_BYTES
# define PML_CONV_PREFETCH_BYTES (1024 * 1024)
#endif

void _nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                            const nmod_poly_mat_t pmat,
                                            slong order,
                                            int kern,
                                            int cmaj,
                                            int pf)
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
       is visited in long sequential runs and the other one cache line at a
       time, from rows that are separately allocated.

       Entry-major is the default: it reads the input polynomials in long
       sequential streams, which the hardware prefetchers follow for free, and
       its scattered side is *stores* of one whole, 64-byte aligned cache line
       each -- which is only cheap because the coefficients of an
       nmod_mat_poly are aligned (see ::_nmod_mat_poly_coeff_alloc). With
       16-byte aligned coefficients, as a plain flint_calloc leaves them,
       every one of those stores straddled two cache lines and the ranking was
       the other way round.

       Measured on a Skylake-SP, entry-major is the faster schedule on every
       shape tested from 2x2 to 128x128 and lengths 8 to 4096 but one isolated
       pocket (1024 entries, length 2048, where the two power-of-two strides
       conflict); its immediate neighbours in both directions prefer
       entry-major again, so no special case is made for it. */
    if (cmaj < 0 || cmaj > 1)
        cmaj = 0;

    /* Prefetch the scattered side -- the output matrices under the
       entry-major schedule, the input polynomials under the other one -- as
       soon as it no longer fits in the first cache levels. Pointless, and a
       small loss, when that side has fewer rows than the prefetch distance,
       since the prefetches then all land on rows already in flight. */
    if (pf < 0 || pf > 1)
        pf = ((double) r * (double) c * (double) order * sizeof(ulong)
                  >= (double) PML_CONV_PREFETCH_BYTES)
             && ((cmaj ? nrows : order) > _PML_CONV_PFD);

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
        _pml_conv(matp->coeffs, order, src, slen, r * c, kern, cmaj, pf);
    }
    else
    {
        /* padded or otherwise nonstandard stride: one transposition per row */
        nn_ptr * dst = (nn_ptr *) flint_malloc(order * sizeof(nn_ptr));
        for (slong i = 0; i < r; i++)
        {
            for (slong k = 0; k < order; k++)
                dst[k] = matp->coeffs[k] + i * matp->stride;
            _pml_conv(dst, order, src + i * c, slen + i * c, c, kern, cmaj, pf);
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
    _nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, order, -1, -1, -1);
}
