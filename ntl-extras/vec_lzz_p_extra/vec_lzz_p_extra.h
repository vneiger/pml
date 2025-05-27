#ifndef __VEC_LZZ_P_EXTRA_H
#define __VEC_LZZ_P_EXTRA_H

#include <NTL/mat_lzz_p.h>
#include <NTL/vec_lzz_p.h>
#include "util.h"

PML_OPEN_NNS
NTL_USE_NNS

/** \brief Additional functions for vectors over `zz_p`
 *
 * \file vec_lzz_p_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-18
 *
 */


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


PML_CLOSE_NNS

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
