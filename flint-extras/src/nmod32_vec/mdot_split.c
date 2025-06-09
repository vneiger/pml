#include "flint/nmod_vec.h"

#include "nmod32_vec.h"

void _nmod32_vec_mdot_split(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                            slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    for (slong i = 0; i < nrows; i++)
        mv[i] = _nmod32_vec_dot_split(mat + i*stride, vec, len, mod, pow2_precomp);
}

#if PML_HAVE_AVX2

void _nmod32_vec_mdot_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                 slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    for (slong i = 0; i < nrows; i++)
        mv[i] = _nmod32_vec_dot_split_avx2(mat + i*stride, vec, len, mod, pow2_precomp);
}

#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512

void _nmod32_vec_mdot_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                   slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    for (slong i = 0; i < nrows; i++)
        mv[i] = _nmod32_vec_dot_split_avx512(mat + i*stride, vec, len, mod, pow2_precomp);
}

#endif  /* PML_HAVE_AVX512 */
