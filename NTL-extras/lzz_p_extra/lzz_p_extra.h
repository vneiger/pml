#ifndef __LZZ_P_EXTRA__H
#define __LZZ_P_EXTRA__H

NTL_CLIENT

/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a is a non-unit                                      */
/*------------------------------------------------------------*/
long order(const zz_p& a);

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* a = 0 if no element was found                              */
/* does (by default) 100 trials                               */
/* by default, asks that all (a^i-1) are units, i=1..ord-1    */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord, long nb_trials = 100, long strong = 1);

inline zz_p element_of_order(long ord, long nb_trials = 100, long strong = 1)
{
    zz_p a;
    element_of_order(a, ord, nb_trials, strong);
    return a;
}

/*------------------------------------------------------------*/
/* 1 if the current prime can be used as an FFT prime         */
/*------------------------------------------------------------*/
long is_FFT_prime();

/*------------------------------------------------------------*/
/* 1 if the current prime can be used as an FFT prime         */
/* with transforms of size up to 2^k                          */
/*------------------------------------------------------------*/
long is_FFT_ready(long k);

/*------------------------------------------------------------*/
/* number of bits for "small primes"                          */
/*------------------------------------------------------------*/
#define SMALL_PRIME_SIZE 25

/*------------------------------------------------------------*/
/* 3 types of primes                                          */
/*------------------------------------------------------------*/
#define TYPE_FFT_PRIME 0
#define TYPE_SMALL_PRIME 1
#define TYPE_LARGE_PRIME 2

/*------------------------------------------------------------*/
/* finds what kind of prime we are using                      */
/*------------------------------------------------------------*/
inline long type_of_prime()
{
    if (is_FFT_prime())
        return TYPE_FFT_PRIME;
    long p = zz_p::modulus();
    if (NumBits(p) <= SMALL_PRIME_SIZE)
        return TYPE_SMALL_PRIME;
    else
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
