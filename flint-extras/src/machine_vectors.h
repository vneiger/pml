#ifndef __MACHINE_VECTORS__H
#define __MACHINE_VECTORS__H

#include <flint/flint.h>

#include "pml.h"

#if PML_HAVE_MACHINE_VECTORS
# include <flint/machine_vectors.h>
#endif

#if PML_HAVE_AVX2
# if defined(__GNUC__)
#  include <immintrin.h>
# elif defined(_MSC_VER)
#  include <intrin.h>
# endif
#endif

/*-------------------------------------------------*/
/* some basic additions to FLINT's machine vectors */
/*-------------------------------------------------*/

#if PML_HAVE_MACHINE_VECTORS


#if PML_HAVE_AVX2  /* already in flint's machine_vectors for NEON/ARM64 */
/* returns a + b mod n, assuming a,b reduced mod n            */
FLINT_FORCE_INLINE vec1n vec1n_addmod(vec1n a, vec1n b, vec1n n)
{
    return n - b > a ? a + b : a + b - n;
}
#endif  /* PML_HAVE_AVX2 */

/* returns a + b mod n, assuming a,b reduced mod n            */
FLINT_FORCE_INLINE vec1d vec1d_addmod(vec1d a, vec1d b, vec1d n)
{
    return a + b - n >= 0 ? a + b - n : a + b;
}

/* returns a + b mod n, assuming a,b reduced mod n            */
FLINT_FORCE_INLINE vec4d vec4d_addmod(vec4d a, vec4d b, vec4d n)
{
    return vec4d_reduce_2n_to_n(vec4d_add(a, b), n);
}

/* loads a vec4n from a and converts it to double             */
#if PML_HAVE_AVX512
FLINT_FORCE_INLINE vec4d vec4d_load_unaligned_nn_ptr(nn_ptr a)
{
    return  _mm256_setr_m128d( _mm_cvtepi64_pd(_mm_loadu_si128((vec2n *) a)),
                               _mm_cvtepi64_pd(_mm_loadu_si128((vec2n *) (a + 2))) );
}
#else
FLINT_FORCE_INLINE vec4d vec4d_load_unaligned_nn_ptr(nn_ptr a)
{
    return vec4n_convert_limited_vec4d(vec4n_load_unaligned(a));
}
#endif

/* converts a vec4d to vec4n and stores it                    */
FLINT_FORCE_INLINE void vec4d_store_unaligned_nn_ptr(nn_ptr dest, vec4d a)
{
    vec4n_store_unaligned(dest, vec4d_convert_limited_vec4n(a));
}

#endif



/*----------------------------------------------------*/
/* some AVX-only additions to FLINT's machine vectors */
/*----------------------------------------------------*/

#if PML_HAVE_AVX2

FLINT_FORCE_INLINE void vec4n_store_aligned(ulong* z, vec4n a)
{
    _mm256_store_si256((__m256i*) z, a);
}

/* reduce_to_pm1n(a, n, ninv): return a mod n in (-n,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_to_pm1no(vec2d a, vec2d n, vec2d ninv)
{
    return _mm_fnmadd_pd(_mm_round_pd(_mm_mul_pd(a, ninv), 4), n, a);
}

/* reduce_pm1no_to_0n(a, n): return a mod n in [0,n) assuming a in (-n,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_pm1no_to_0n(vec2d a, vec2d n)
{
    return _mm_blendv_pd(a, _mm_add_pd(a, n), a);
}

/* reduce_to_0n(a, n, ninv): return a mod n in [0,n) */
FLINT_FORCE_INLINE vec2d vec2d_reduce_to_0n(vec2d a, vec2d n, vec2d ninv)
{
    return vec2d_reduce_pm1no_to_0n(vec2d_reduce_to_pm1no(a, n, ninv), n);
}

FLINT_FORCE_INLINE vec2d vec2d_set_d2(double a1, double a0)
{
    return _mm_set_pd(a0, a1);
}

#define vec4n_bit_shift_right_45(a) vec4n_bit_shift_right((a), 45)


// avx2 horizontal sum  (already in FLINT)
FLINT_FORCE_INLINE ulong _mm256_hsum(__m256i a) {
    __m256i a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    __m256i sum_lo = _mm256_add_epi64(a, a_hi);
    __m128i sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (ulong) _mm_cvtsi128_si64(sum);
}

#endif /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512
// avx512 horizontal sum
FLINT_FORCE_INLINE ulong _mm512_hsum(__m512i a) {
    return _mm512_reduce_add_epi64(a);
}
#endif /* PML_HAVE_AVX512 */


#endif /* ifndef __MACHINE_VECTORS__H */
