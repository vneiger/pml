#ifndef __LZZ_P_EXTRA__H
#define __LZZ_P_EXTRA__H

/** \brief Additional functions for `lzz_p`
 *
 * \file lzz_p_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-15
 *
 */

/** Determines the maximum number of bits for "small" primes */
#define SMALL_PRIME_SIZE 25
/** Prime of type _FFT prime_ (2 other types: small prime, larges prime) */
#define TYPE_FFT_PRIME 0
/** Prime of type _small prime_ (2 other types: FFT prime, larges prime) */
#define TYPE_SMALL_PRIME 1
/** Prime of type _large prime_ (2 other types: FFT prime, small prime) */
#define TYPE_LARGE_PRIME 2

#include <NTL/lzz_p.h>

NTL_CLIENT

/** Computes and returns the multiplicative order of `a` if it is a unit;
 * returns `-1` otherwise */
long order(const zz_p & a);

/** Computes an element `a` of order at least `ord` using `nb_trials` trials
 * (default: 100 trials); if no element is found, set `a` to `0`. If `strong`
 * is `true` (the default), then it is further required that `a^i - 1` be a
 * unit for all `i` from `1` to `ord-1`. */
void element_of_order(
                      zz_p & a,
                      long ord,
                      long nb_trials = 100,
                      bool strong = true
                     );

/** Computes and returns an element of order at least `ord` using `nb_trials`
 * trials (default: 100 trials); if no element is found, set `a` to `0`. If
 * `strong` is `true` (the default), then it is further required that the
 * output element `a` be such that `a^i - 1` is a unit for all `i` from `1` to
 * `ord-1`. */
inline zz_p element_of_order(
                             long ord,
                             long nb_trials = 100,
                             bool strong = true
                            )
{ zz_p a; element_of_order(a, ord, nb_trials, strong); return a; }

/** Returns `true` if the current modulus can be used as an FFT prime; returns
 * `false` otherwise */
inline bool is_FFT_prime()
{ return (zz_pInfo->p_info != NULL); }

/** Returns `true` if the current modulus can be used as an FFT prime with
 * transforms of size up to `2^k`; returns `false` otherwise */
inline bool is_FFT_ready(long k)
{ return (zz_pInfo->p_info != NULL && k <= zz_pInfo->MaxRoot); }

/** Returns the type of prime we are using (see macros #TYPE_FFT_PRIME,
 * #TYPE_SMALL_PRIME, #TYPE_LARGE_PRIME) */
inline long type_of_prime()
{
    if (zz_pInfo->p_info != NULL) // i.e. if is_FFT_prime()
        return TYPE_FFT_PRIME;
    if (NumBits(zz_p::modulus()) <= SMALL_PRIME_SIZE)
        return TYPE_SMALL_PRIME;
    return TYPE_LARGE_PRIME;
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
