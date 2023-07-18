#ifndef __NMOD_VEC_EXTRA__H
#define __NMOD_VEC_EXTRA__H

/** \brief Extra functions for vectors over Z/nZ
 *
 * \file nmod_vec_extra.h
 * \version 0.0
 * \date 2023-01-23
 *
 * Some functions to deal with vectors over `nmod`.
 *
 */

#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#ifdef __cplusplus
extern "C" {
#endif


// TODO augment/add documentation

/** Random */
// to be completed: random sparse? randtest small? randtest nonzero ? ...

/** Fills the entries `0`, .., `len-1` of vector with uniformly random entries.
 * Vector must already be allocated with length at least `len`. */
void _nmod_vec_rand(mp_ptr vec,
            		flint_rand_t state,
            		slong len,
            		nmod_t mod);




/** Input/Output */
// to be completed: print to file, print to sagemath..

/** Prints the entries `0`, .., `len-1` of vector `vec`. Vector must already be
 * initialized, with length at least `len`. */
void _nmod_vec_print_pretty(mp_ptr vec,
                            slong len,
                            nmod_t mod);




/*------------------------------------------------------------*/
/* y[i] = a*x[i] mod mod, i=0..len-1                          */
/*------------------------------------------------------------*/
static inline void _nmod_vec_scalar_mul(mp_ptr y, mp_srcptr x, slong len, mp_limb_t a, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
    {
        y[i] = nmod_mul(x[i], a, mod);
    }
}

#ifdef __cplusplus
}
#endif

#endif // ifndef  __NMOD_VEC_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
