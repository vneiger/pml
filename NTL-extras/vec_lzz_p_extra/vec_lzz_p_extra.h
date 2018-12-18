#ifndef __VEC_LZZ_P_EXTRA_H
#define __VEC_LZZ_P_EXTRA_H

/** \brief Additional functions for vectors over `zz_p`
 *
 * \file vec_lzz_p_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-18
 *
 */

#include <NTL/mat_lzz_p.h>
#include <NTL/vec_lzz_p.h>
#include <NTL/version.h>

NTL_CLIENT

#if ( (NTL_MAJOR_VERSION < 10) || ((NTL_MAJOR_VERSION == 10) && (NTL_MINOR_VERSION < 4)) )
/** Builds and returns a random vector over `zz_p` of length `n` (defined only
 * for versions of NTL that do not include this feature) */
inline Vec<zz_p> random_vec_zz_p(long n)
{ Vec<zz_p> x; random(x, n); return x; }
#endif

#if ( (NTL_MAJOR_VERSION < 10) || ((NTL_MAJOR_VERSION == 10) && (NTL_MINOR_VERSION < 4)) )
/** Builds and returns a random matrix vector over `zz_p`'s of dimensions `d x
 * e` (defined only for versions of NTL that do not include this feature) */
inline Mat<zz_p> random_mat_zz_p(long d, long e)
{ Mat<zz_p> x; random(x, d, e); return x; }
#endif

/** Computes the vector `invA`, which is `A` with each entry inverted. The OUT
 * parameter `invA` may alias the IN parameter `A`. */
void inv_naive(Vec<zz_p> & invA, const Vec<zz_p> & A);

/** Computes the vector `invA`, which is `A` with each entry inverted. The OUT
 * parameter `invA` may alias the IN parameter `A`. */
void inv(Vec<zz_p> & invA, const Vec<zz_p> & A);

/** Computes and returns the vector which is `A` with each entry inverted. */
inline Vec<zz_p> inv(const Vec<zz_p> & A)
{ Vec<zz_p> x; inv(x, A); return x; }

/** Computes a vector `precon` containing the `mulmod_precon_t` object
 * for each entry of the input vector `A`  */
void precomp(Vec<mulmod_precon_t> & precon, const Vec<zz_p> & A);

/** Computes and returns a vector containing the `mulmod_precon_t` object for
 * each entry of the input vector `A`  */
inline Vec<mulmod_precon_t> precomp(const Vec<zz_p> & A)
{ Vec<mulmod_precon_t> precon; precomp(precon, A); return precon; }


#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
