#include "flint/flint.h"
#include "flint/longlong.h"
#include "flint/longlong_div_gnu.h"
#include "flint/profiler.h"
#include "flint/nmod_vec.h"
#include "flint/ulong_extras.h"

// must be a multiple of 6
#define LEN 1002

typedef struct
{
    flint_bitcnt_t bits;
} info_t;

/*-----------------------------*/
/* modular reduction functions */
/*-----------------------------*/

// modular reduction with precomputation,
// using n_pr = floor(2**FLINT_BITS / n),
// e.g. computed via n_mulmod_precomp_shoup(1L, n)
static inline
ulong n_modred_precomp(ulong a, ulong n_pr, ulong n)
{
    // the next call is equivalent to the code below,
    // and any decent compiler should compile it right
    return n_mulmod_shoup(1L, a, n_pr, n);

    //ulong p_hi, p_lo;

    //umul_ppmm(p_hi, p_lo, n_pr, a);

    //a -= p_hi * n;

    //if (a >= n)
    //    a -= n;

    //return a;
}

// modular reduction with precomputation,
// using n_pr =
// Explanation:
// Suppose n has nn bits, 2**(nn-1) <= n < 2**nn
//         a has <= aa bits, a < 2**aa, and we assume nn <= aa
// Take some parameter tt >= 0,
// and compute n_pr = floor(2**tt / n); this can be precomputed
// for fixed n and bound aa.
// One has n_pr <= 2**tt / n <= 2**tt / 2**(nn-1), with equalities iff n is a
// power of 2. Typically aa = FLINT_BITS and one wants n_pr < 2**aa:
// this holds as soon as tt <= nn+aa-1 and n is not a power of 2.
//
// Then, to reduce any a with aa bits,
//    q <- floor(a * n_pr / 2**tt)
//          [if aa = FLINT_BITS and n_pr < 2**aa, this is obtained by shifting
//          by tt - aa the high part of the double-word product a * n_pr;
//          for the case where n is a power of 2, i.e. n = 2**(nn-1), one can
//          check that taking tt = aa, one can compute q without modifying the
//          algorithm and it will be correct q = floor(a/n) in that case]
//    r <- a - q * n
//          [in the above context, this is a single-word multiplication and a
//          subtraction]
//
// Write a = Qn + R with 0 <= R <= n-1
// How far is q from Q = floor(a/n), and r from R = a mod n?
//
//    . By construction, q <= floor(a * (2**tt / n) / 2**tt) = Q.
//
//    . a * floor(2**tt / n) == a ((2**tt + k) / n  - 1),
//      where k >= 1 is the smallest integer such that n divides 2**tt + k;
//      thus k <= n, with equality iff n is a power of 2.
//
//      We get q = floor(a/n + k * a / (n * 2**tt) - a / 2**tt)     [1]
//               = Q + floor(R/n +  k * R / (n * 2**tt) + (k*Q - a) / 2**tt)
//               = Q + floor( 2**(-tt) * (R/n * (2**tt + k) + k*Q - a) )  =: F
//
//      Looking at [1], it seems a good idea to require tt >= aa, because it
//      ensures a / 2**tt < 1; then, in all cases, F >= Q - 1, meaning that
//      either q = Q and r = R, or q = Q-1 and r = r + n.
//
//      We can say more about the floor in the definition of F:
//          it is >= 0  iff  R/n * (2**tt + k) + k*Q - a >= 0
//                      iff  R * (2**tt + k) >= (a - k*Q) * n = a*(n-k) + k*R
//                      iff  R * 2**tt  >=  a * (n-k).
//
//      We seek tt so as to maximize the chance of this happening.
//
//      A particularly bad case is when n divides a (Q = a/n and R = 0). Then
//      whatever choice for tt, q = Q + floor((k*Q - a) / 2**tt). In most cases
//      (as soon as a > 0 and n is not a power of 2), one has k*Q - a < n*Q - a
//      = 0, and thus q = Q-1 and r = R+n = n.
//
//      On the other hand, for all a such that R is >= 2, we just need to
//      ensure 2**(tt+1) >= a * (n-k), which holds for tt = nn+aa-1 (this tt
//      also guarantees that the precomputed n_pr fits in aa bits).
//
//      Say we set tt = nn+aa-1. We know that a's with R >= 2 are fine (r is
//      correct), and a's with R = 0 are understood (we get r = 0 or r = n).
//      There remains the case R == 1. The inequality 2**tt >= a * (n-1) may
//      not hold for some a, typically those close to 2**aa, and we may get r =
//      n+1 instead of r = 1.
//       Example: n = 1095128386781456219, which has nn = 60 bits and divides 2**123 + 1
//                a = 9*n + 1 which has aa = 64 bits, and (Q,R) = (9,1)
//                here nn+aa-1 = 123, and if we take tt = 123,
//                2**123 < a * (n-1) = 2**(123.0215...)
//                One can check that indeed the formulas yield q == 8 and r = n+1
//
//      A further interesting remark is that, if n-k happens to have < nn bits,
//      then R * 2**tt  >=  a * (n-k) holds for tt = nn+aa-1 even when R = 1.
//      So, in that case, only a divisible by n will not lead to a reduced
//      result, and the result is always in [0,n+1). Note that the property "n-k
//      has < nn bits" is known at precomputation stage: it depends on n and
//      aa (provided tt = nn+aa-1), but it does not depend on a. In contexts where
//      one can choose the modulus, one could hopefully ensure this property,
//      so as to get results in [0,n+1).
//
//      All in all, this may not give the correct reduced remainder, but still
//      gives the remainder possibly +n, and always in the range [0,n+2) (i.e.
//      at most 2 away from the "ideal" range [0,n)). For lazy techniques / delayed reductions
//      where one accumulates computation results that might be outside the ideal range,
//      it is much better to have this result in [0,n+2) than the only guarantee
//      [0,2*n). In other words, in contexts where this reduced a mod n is not
//      the end result, and is to be used for further computations, this very
//      small excess may well be manageable without requiring a test such as
//         if (r >= n) r -= n;

// below we take aa = 64, tt = nbits(n) + 64 - 1
// sets shift to tt - aa = nbits(n) - 1,
// sets n_pr to floor(2**tt / n)
// TODO handle special cases: n == 1, and power of 2
static inline
void n_precomp_modred_precomp2(flint_bitcnt_t * shift, ulong * n_pr, ulong n)
{
    ulong tmp;

    *shift = FLINT_BIT_COUNT(n) - 1;
    udiv_qrnnd(*n_pr, tmp, UWORD(1) << *shift, UWORD(0), n);
}

static inline
ulong n_modred_precomp2(ulong a, ulong n_pr, flint_bitcnt_t shift, ulong n)
{
    // q <- floor(a * n_pr / 2**tt)
    ulong q, tmp;

    umul_ppmm(q, tmp, a, n_pr);
    q >>= shift;

    // r <- a - q * n
    a -= q * n;
    if (a >= n)
        a -= n;

    return a;
}

static inline
ulong n_modred_precomp2_lazy(ulong a, ulong n_pr, flint_bitcnt_t shift, ulong n)
{
    // q <- floor(a * n_pr / 2**tt)
    ulong q, tmp;

    umul_ppmm(q, tmp, a, n_pr);
    q >>= shift;

    // r <- a - q * n
    a -= q * n;

    return a;
}

#define SAMPLE(fun)                                       \
void sample_##fun(void * arg, ulong count)                \
{                                                         \
    ulong n;                                              \
    nmod_t mod;                                           \
    info_t * info = (info_t *) arg;                       \
    flint_bitcnt_t bits = info->bits;                     \
    nn_ptr vec = _nmod_vec_init(LEN);                     \
    nn_ptr vec2 = _nmod_vec_init(LEN);                    \
    FLINT_TEST_INIT(state);                               \
                                                          \
    for (ulong j = 0; j+2 < LEN; j+=3)                    \
    {                                                     \
        /* already reduced, for NMOD_RED3 */              \
        vec[j+0] = n_randbits(state, bits-1);             \
        vec[j+1] = n_randlimb(state);                     \
        vec[j+2] = n_randlimb(state);                     \
    }                                                     \
                                                          \
    prof_start();                                         \
    for (ulong i = 0; i < count; i++)                     \
    {                                                     \
        n = n_randbits(state, bits);                      \
        nmod_init(&mod, n);                               \
                                                          \
        prof_nmod_vec_reduce_##fun(vec2, vec, LEN, mod);  \
    }                                                     \
    prof_stop();                                          \
                                                          \
    _nmod_vec_clear(vec);                                 \
    _nmod_vec_clear(vec2);                                \
    FLINT_TEST_CLEAR(state);                              \
}

/*------------------------------*/
/* reduce res[i] = vec[i] % n   */
/*------------------------------*/

void prof_nmod_vec_reduce_nmod_red(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    for (slong i = 0 ; i < len; i++)
        NMOD_RED(res[i], vec[i], mod);
}

void prof_nmod_vec_reduce_nmod_red_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    slong i;
    for (i = 0; i+3 < len; i+=4)
    {
        NMOD_RED(res[i+0], vec[i+0], mod);
        NMOD_RED(res[i+1], vec[i+1], mod);
        NMOD_RED(res[i+2], vec[i+2], mod);
        NMOD_RED(res[i+3], vec[i+3], mod);
    }
    for ( ; i < len; i++)
        NMOD_RED(res[i], vec[i], mod);
}

// note: in this particular instance of n_mulmod_shoup, no restriction on n
void prof_nmod_vec_reduce_precomp(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    const ulong one_precomp = n_mulmod_precomp_shoup(1L, mod.n);
    for (slong i = 0 ; i < len; i++)
        res[i] = n_modred_precomp(vec[i], one_precomp, mod.n);
}

// note: in this particular instance of n_mulmod_shoup, no restriction on n
void prof_nmod_vec_reduce_precomp_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    const ulong one_precomp = n_mulmod_precomp_shoup(1L, mod.n);
    slong i;
    for (i = 0 ; i+3 < len; i+=4)
    {
        res[i+0] = n_modred_precomp(vec[i+0], one_precomp, mod.n);
        res[i+1] = n_modred_precomp(vec[i+1], one_precomp, mod.n);
        res[i+2] = n_modred_precomp(vec[i+2], one_precomp, mod.n);
        res[i+3] = n_modred_precomp(vec[i+3], one_precomp, mod.n);
    }
    for ( ; i < len; i++)
        res[i] = n_modred_precomp(vec[i], one_precomp, mod.n);
}

void prof_nmod_vec_reduce_precomp2(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    flint_bitcnt_t shift;
    ulong n_pr;
    n_precomp_modred_precomp2(&shift, &n_pr, mod.n);
    for (slong i = 0 ; i < len; i++)
        res[i] = n_modred_precomp2(vec[i], n_pr, shift, mod.n);
}

void prof_nmod_vec_reduce_precomp2_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    flint_bitcnt_t shift;
    ulong n_pr;
    n_precomp_modred_precomp2(&shift, &n_pr, mod.n);
    slong i;
    for (i = 0; i+3 < len; i+=4)
    {
        res[i+0] = n_modred_precomp2(vec[i+0], n_pr, shift, mod.n);
        res[i+1] = n_modred_precomp2(vec[i+1], n_pr, shift, mod.n);
        res[i+2] = n_modred_precomp2(vec[i+2], n_pr, shift, mod.n);
        res[i+3] = n_modred_precomp2(vec[i+3], n_pr, shift, mod.n);
    }
    for ( ; i < len; i++)
        res[i] = n_modred_precomp2(vec[i], n_pr, shift, mod.n);
}


/*----------------------------------------*/
/*        res[i] = reduce double word     */
/* (hi, lo) = (vec[i], vec[i+1]) modulo n */
/*----------------------------------------*/
// warning: using mulmod_shoup here restricts n to 63 bits

void prof_nmod_vec_reduce_nmod2_red2(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    for (slong i = 0 ; i+1 < len; i+=2)
        NMOD2_RED2(res[i], vec[i], vec[i+1], mod);
}

void prof_nmod_vec_reduce_nmod2_red2_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    slong i;
    for (i = 0 ; i+7 < len; i+=8)
    {
        NMOD2_RED2(res[i+0], vec[i+0], vec[i+1], mod);
        NMOD2_RED2(res[i+2], vec[i+2], vec[i+3], mod);
        NMOD2_RED2(res[i+4], vec[i+4], vec[i+5], mod);
        NMOD2_RED2(res[i+6], vec[i+6], vec[i+7], mod);
    }
    for ( ; i+1 < len; i+=2)
        NMOD2_RED2(res[i], vec[i], vec[i+1], mod);
}

void prof_nmod_vec_reduce_precomp22(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n
    W_pr = n_mulmod_precomp_shoup(W, mod.n);

    flint_bitcnt_t shift;
    ulong n_pr;
    n_precomp_modred_precomp2(&shift, &n_pr, mod.n);

    // naive: res[i] = (vec[i] * W mod n + vec[i+1] mod n) mod n
    for (slong i = 0 ; i+1 < len; i+=2)
    {
        const ulong a = n_mulmod_shoup(W, vec[i], W_pr, mod.n);
        //const ulong b = n_modred_precomp(vec[i+1], one_precomp, mod.n);
        const ulong b = n_modred_precomp2_lazy(vec[i+1], n_pr, shift, mod.n);
        res[i] = _nmod_add(a, b, mod);
    }
}

void prof_nmod_vec_reduce_precomp22_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n
    W_pr = n_mulmod_precomp_shoup(W, mod.n);

    flint_bitcnt_t shift;
    ulong n_pr;
    n_precomp_modred_precomp2(&shift, &n_pr, mod.n);

    // naive: res[i] = (vec[i] * W mod n + vec[i+1] mod n) mod n
    slong i;
    for (i = 0; i+7 < len; i+=8)
    {
        res[i+0] = _nmod_add(n_mulmod_shoup(W, vec[i+0], W_pr, mod.n),
                             //n_modred_precomp(vec[i+1], one_precomp, mod.n),
                             n_modred_precomp2(vec[i+1], n_pr, shift, mod.n),
                             mod);
        res[i+2] = _nmod_add(n_mulmod_shoup(W, vec[i+2], W_pr, mod.n),
                             //n_modred_precomp(vec[i+3], one_precomp, mod.n),
                             n_modred_precomp2(vec[i+3], n_pr, shift, mod.n),
                             mod);
        res[i+4] = _nmod_add(n_mulmod_shoup(W, vec[i+4], W_pr, mod.n),
                             //n_modred_precomp(vec[i+5], one_precomp, mod.n),
                             n_modred_precomp2(vec[i+5], n_pr, shift, mod.n),
                             mod);
        res[i+6] = _nmod_add(n_mulmod_shoup(W, vec[i+6], W_pr, mod.n),
                             //n_modred_precomp(vec[i+7], one_precomp, mod.n),
                             n_modred_precomp2(vec[i+7], n_pr, shift, mod.n),
                             mod);
    }

    for ( ; i+1 < len; i+=2)
    {
        res[i] = _nmod_add(n_mulmod_shoup(W, vec[i], W_pr, mod.n),
                           //n_modred_precomp(vec[i+1], one_precomp, mod.n),
                            n_modred_precomp2(vec[i+1], n_pr, shift, mod.n),
                           mod);
    }
}

// version without unrolling is slower
void prof_nmod_vec_reduce_precomp22_bis(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n
    W_pr = n_mulmod_precomp_shoup(W, mod.n);

    // naive: res[i] = (vec[i] * W mod n + vec[i+1] mod n) mod n
    slong i;
    for (i = 0; i+7 < len; i+=8)
    {
        ulong a0, b0, a1, b1, a2, b2, a3, b3;
        a0 = n_mulmod_shoup(W, vec[i+0], W_pr, mod.n);
        a1 = n_mulmod_shoup(W, vec[i+2], W_pr, mod.n);
        a2 = n_mulmod_shoup(W, vec[i+4], W_pr, mod.n);
        a3 = n_mulmod_shoup(W, vec[i+6], W_pr, mod.n);

        b0 = vec[i+1];
        b1 = vec[i+3];
        b2 = vec[i+5];
        b3 = vec[i+7];
        if (b0 >= (mod.n << (mod.norm-1)))
            b0 -= (mod.n << (mod.norm-1));
        if (b1 >= (mod.n << (mod.norm-1)))
            b1 -= (mod.n << (mod.norm-1));
        if (b2 >= (mod.n << (mod.norm-1)))
            b2 -= (mod.n << (mod.norm-1));
        if (b3 >= (mod.n << (mod.norm-1)))
            b3 -= (mod.n << (mod.norm-1));

        res[i+0] = n_modred_precomp(a0+b0, one_precomp, mod.n);
        res[i+2] = n_modred_precomp(a1+b3, one_precomp, mod.n);
        res[i+4] = n_modred_precomp(a2+b3, one_precomp, mod.n);
        res[i+6] = n_modred_precomp(a3+b3, one_precomp, mod.n);
    }
    for ( ; i+1 < len; i+=2)
    {
        ulong a, b;
        a = n_mulmod_shoup(W, vec[i], W_pr, mod.n);
        b = vec[i+1];
        if (b >= (mod.n << (mod.norm-1)))
            b -= (mod.n << (mod.norm-1));
        res[i] = n_modred_precomp(a+b, one_precomp, mod.n);
    }
}

// version without unrolling is slower
void prof_nmod_vec_reduce_precomp22_ter(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n
    W_pr = n_mulmod_precomp_shoup(W, mod.n);

    flint_bitcnt_t shift;
    ulong n_pr;
    n_precomp_modred_precomp2(&shift, &n_pr, mod.n);

    // naive: res[i] = (vec[i] * W mod n + vec[i+1] mod n) mod n
    slong i;
    for (i = 0; i+7 < len; i+=8)
    {
        ulong a0 = n_mulmod_shoup(W, vec[i+0], W_pr, mod.n) + n_modred_precomp2_lazy(vec[i+1], n_pr, shift, mod.n);
        ulong a1 = n_mulmod_shoup(W, vec[i+2], W_pr, mod.n) + n_modred_precomp2_lazy(vec[i+3], n_pr, shift, mod.n);
        ulong a2 = n_mulmod_shoup(W, vec[i+4], W_pr, mod.n) + n_modred_precomp2_lazy(vec[i+5], n_pr, shift, mod.n);
        ulong a3 = n_mulmod_shoup(W, vec[i+6], W_pr, mod.n) + n_modred_precomp2_lazy(vec[i+7], n_pr, shift, mod.n);

        if (a0 >= mod.n) a0 -= mod.n;
        if (a1 >= mod.n) a1 -= mod.n;
        if (a2 >= mod.n) a2 -= mod.n;
        if (a3 >= mod.n) a3 -= mod.n;

        res[i+0] = a0;
        res[i+2] = a1;
        res[i+4] = a2;
        res[i+6] = a3;
    }
    for ( ; i+1 < len; i+=2)
    {
        ulong a, b;
        a = n_mulmod_shoup(W, vec[i], W_pr, mod.n);
        b = vec[i+1];
        if (b >= (mod.n << (mod.norm-1)))
            b -= (mod.n << (mod.norm-1));
        res[i] = n_modred_precomp(a+b, one_precomp, mod.n);
    }
}

/*------------------------------------------------------*/
/*        res[i] = reduce 3-limb word                   */
/* (hi, me, lo) = (vec[i], vec[i+1], vec[i+2]) modulo n */
/* the nmod_red3 variant requires vec[i] < n            */
/*------------------------------------------------------*/

void prof_nmod_vec_reduce_nmod_red3(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    for (slong i = 0 ; i+2 < len; i+=3)
        NMOD_RED3(res[i], vec[i], vec[i+1], vec[i+2], mod);
}

void prof_nmod_vec_reduce_nmod_red3_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    slong i;
    for (i = 0; i+5 < len; i+=6)
    {
        NMOD_RED3(res[i], vec[i], vec[i+1], vec[i+2], mod);
        NMOD_RED3(res[i+3], vec[i+3], vec[i+4], vec[i+5], mod);
    }
    for ( ; i+2 < len; i+=3)
    {
        NMOD_RED3(res[i], vec[i], vec[i+1], vec[i+2], mod);
    }
}

// warning: using mulmod_shoup here restricts n to 63 bits
void prof_nmod_vec_reduce_precomp3(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr_quo, W_pr_rem, W2, W2_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n

    n_mulmod_precomp_shoup_quo_rem(&W_pr_quo, &W_pr_rem, W, mod.n);
    n_mulmod_and_precomp_shoup(&W2, &W2_pr, W, W, W_pr_quo, W_pr_rem, W_pr_quo, mod.n);
    // --> W2 = remainder of division of 2**(2*FLINT_BITS) by n

    // naive: res[i] = (vec[i] * W2 mod n + vec[i+1] * W mod n + vec[i+2] mod n) mod n
    for (slong i = 0 ; i+2 < len; i+=3)
    {
        res[i] = _nmod_add(n_mulmod_shoup(W2, vec[i], W2_pr, mod.n),
                           _nmod_add(n_mulmod_shoup(W, vec[i+1], W_pr_quo, mod.n),
                                     n_modred_precomp(vec[i+2], one_precomp, mod.n),
                                     mod),
                           mod);
    }
}

void prof_nmod_vec_reduce_precomp3_bis(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr_quo, W_pr_rem, W2, W2_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n

    n_mulmod_precomp_shoup_quo_rem(&W_pr_quo, &W_pr_rem, W, mod.n);
    n_mulmod_and_precomp_shoup(&W2, &W2_pr, W, W, W_pr_quo, W_pr_rem, W_pr_quo, mod.n);
    // --> W2 = remainder of division of 2**(2*FLINT_BITS) by n

    // naive: res[i] = (vec[i] * W2 mod n + vec[i+1] * W mod n + vec[i+2] mod n) mod n
    for (slong i = 0 ; i+2 < len; i+=3)
    {
        ulong buf = 0;
        ulong c0 = n_mulmod_shoup(W2, vec[i+0], W2_pr, mod.n);
        ulong c1 = n_mulmod_shoup(W, vec[i+1], W_pr_quo, mod.n);
        c0 = _nmod_add(c0, c1, mod);
        if (vec[i+2] >= (mod.n << (mod.norm-1)))
            buf = (mod.n << (mod.norm-1));
        c0 += vec[i+2] - buf;
        umul_ppmm(c1, buf, one_precomp, c0);
        c0 -= c1 * mod.n;
        if (c0 >= mod.n)
            c0 -= mod.n;
        res[i+0] = c0;
    }
}

void prof_nmod_vec_reduce_precomp3_bis_unroll(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    ulong one_precomp, W, W_pr_quo, W_pr_rem, W2, W2_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    // --> W = remainder of division of 2**FLINT_BITS by n

    n_mulmod_precomp_shoup_quo_rem(&W_pr_quo, &W_pr_rem, W, mod.n);
    n_mulmod_and_precomp_shoup(&W2, &W2_pr, W, W, W_pr_quo, W_pr_rem, W_pr_quo, mod.n);
    // --> W2 = remainder of division of 2**(2*FLINT_BITS) by n

    // naive: res[i] = (vec[i] * W2 mod n + vec[i+1] * W mod n + vec[i+2] mod n) mod n
    for (slong i = 0 ; i+5 < len; i+=6)
    {
        ulong c0 = n_mulmod_shoup(W2, vec[i+0], W2_pr, mod.n);
        ulong c1 = n_mulmod_shoup(W, vec[i+1], W_pr_quo, mod.n);
        c0 = _nmod_add(c0, c1, mod);

        ulong c3 = n_mulmod_shoup(W2, vec[i+3], W2_pr, mod.n);
        ulong c4 = n_mulmod_shoup(W, vec[i+4], W_pr_quo, mod.n);
        c3 = _nmod_add(c0, c1, mod);

        ulong buf = 0;
        if (vec[i+2] >= (mod.n << (mod.norm-1)))
            buf = (mod.n << (mod.norm-1));
        c0 += vec[i+2] - buf;
        buf = 0;
        if (vec[i+5] >= (mod.n << (mod.norm-1)))
            buf = (mod.n << (mod.norm-1));
        c3 += vec[i+5] - buf;

        umul_ppmm(c1, buf, one_precomp, c0);
        c0 -= c1 * mod.n;
        if (c0 >= mod.n)
            c0 -= mod.n;
        res[i+0] = c0;

        umul_ppmm(c4, buf, one_precomp, c3);
        c3 -= c4 * mod.n;
        if (c3 >= mod.n)
            c3 -= mod.n;
        res[i+3] = c3;
    }
}

// Horner style
void prof_nmod_vec_reduce_precomp3_ter(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
{
    // --> W = remainder of division of 2**FLINT_BITS by n
    ulong one_precomp, W, W_pr;
    n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
    W_pr = n_mulmod_precomp_shoup(W, mod.n);

    // Horner for res[i] = (vec[i] * W + vec[i+1]) * W + vec[i+2]
    for (slong i = 0 ; i+2 < len; i+=3)
    {
        ulong v0 = vec[i+0];
        ulong v1 = vec[i+1];
        ulong v2 = vec[i+2];
        if (v1 >= (mod.n << (mod.norm-1)))
            v1 -= (mod.n << (mod.norm-1));
        if (v2 >= (mod.n << (mod.norm-1)))
            v2 -= (mod.n << (mod.norm-1));
        v1 += n_mulmod_shoup(W, v0, W_pr, mod.n);
        v2 += n_mulmod_shoup(W, v1, W_pr, mod.n);
        res[i+0] = n_modred_precomp(v2, one_precomp, mod.n);
    }
}


int check(flint_bitcnt_t bits)
{
    ulong n;
    nmod_t mod;

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < LEN; i++)
    {
        n = n_randbits(state, bits);
        nmod_init(&mod, n);

        ulong a_lo = n_randlimb(state);
        ulong a_me = n_randlimb(state);
        ulong a_hi = n_randlimb(state);
        ulong a_hi_red = n_randint(state, n);

        ulong witness;
        { // 1 limb
            NMOD_RED(witness, a_lo, mod);

            const ulong one_precomp = n_mulmod_precomp_shoup(1L, mod.n);
            ulong candidate = n_modred_precomp(a_lo, one_precomp, mod.n);

            if (witness != candidate)
            { printf("\n\n\nFAILURE 1 limb modred_precomp!!\n\n\n"); return 0; }

            flint_bitcnt_t shift;
            ulong n_pr;
            n_precomp_modred_precomp2(&shift, &n_pr, n);
            candidate = n_modred_precomp2(a_lo, n_pr, shift, n);

            if (witness != candidate)
            { printf("\n\n\nFAILURE 1 limb modred_precomp2!!\n\n\n"); return 0; }
        }

        if (bits < FLINT_BITS)
        { // 2 limbs, direct
            NMOD2_RED2(witness, a_hi, a_lo, mod);

            ulong one_precomp, W, W_pr;
            n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
            W_pr = n_mulmod_precomp_shoup(W, mod.n);

            ulong candidate =
                _nmod_add(n_mulmod_shoup(W, a_hi, W_pr, mod.n),
                           n_modred_precomp(a_lo, one_precomp, mod.n),
                           mod);

            if (witness != candidate)
            { printf("\n\n\nFAILURE 2 limbs direct!!\n\n\n"); return 0; }
        }

        if (bits < FLINT_BITS)
        { // 2 limbs, less direct
            NMOD2_RED2(witness, a_hi, a_lo, mod);

            ulong one_precomp, W, W_pr;
            n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
            W_pr = n_mulmod_precomp_shoup(W, mod.n);

            ulong candidate = n_mulmod_shoup(W, a_hi, W_pr, mod.n);
            ulong correct = 0;
            if (a_lo >= (mod.n << (mod.norm-1)))
                correct = (mod.n << (mod.norm-1));
            candidate = n_modred_precomp(candidate+a_lo-correct, one_precomp, mod.n);

            if (witness != candidate)
            //{ printf("\n\n\nFAILURE 2 limbs!!\n\n\n"); printf("%lu, %lu, %lu, %lu, %lu\n", mod.n, a_lo, a_hi, candidate, witness); return 0; }
            { printf("\n\n\nFAILURE 2 limbs!!\n\n\n"); return 0; }
        }

        if (bits < FLINT_BITS-1)
        { // 2 limbs, less direct, bis
          // <= 62 bits allows to remove one conditional (see below), but
          // this does not seem to gain anything in the vector reduction
          // context at least (tried within precomp22)
            NMOD2_RED2(witness, a_hi, a_lo, mod);

            ulong one_precomp, W, W_pr;
            n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
            W_pr = n_mulmod_precomp_shoup(W, mod.n);

            ulong candidate, tmp, p_hi;
            umul_ppmm(candidate, tmp, W_pr, a_hi);
            candidate = W * a_hi - candidate * mod.n;
            //if (candidate >= mod.n)
                //candidate -= mod.n;
            ulong correct = 0;
            if (a_lo >= (mod.n << (mod.norm-1)))
                correct = (mod.n << (mod.norm-1));
            candidate += a_lo - correct;
            umul_ppmm(p_hi, tmp, one_precomp, candidate);
            candidate -= p_hi * mod.n;
            if (candidate >= mod.n)
                candidate -= mod.n;

            if (witness != candidate)
            //{ printf("\n\n\nFAILURE 2 limbs bis!!\n\n\n"); printf("%lu, %lu, %lu, %lu, %lu\n", mod.n, a_lo, a_hi, candidate, witness); return 0; }
            { printf("\n\n\nFAILURE 2 limbs!!\n\n\n"); return 0; }
        }

        if (bits < FLINT_BITS)
        { // 3 limbs, direct
            NMOD_RED3(witness, a_hi_red, a_me, a_lo, mod);

            ulong one_precomp, W, W_pr_quo, W_pr_rem, W2, W2_pr;
            n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
            // --> W = remainder of division of 2**FLINT_BITS by n

            n_mulmod_precomp_shoup_quo_rem(&W_pr_quo, &W_pr_rem, W, mod.n);
            n_mulmod_and_precomp_shoup(&W2, &W2_pr, W, W, W_pr_quo, W_pr_rem, W_pr_quo, mod.n);
            // --> W2 = remainder of division of 2**(2*FLINT_BITS) by n

            // naive: res[i] = (vec[i] * W2 mod n + vec[i+1] * W mod n + vec[i+2] mod n) mod n
            ulong candidate = _nmod_add(n_mulmod_shoup(W2, a_hi_red, W2_pr, mod.n),
                               _nmod_add(n_mulmod_shoup(W, a_me, W_pr_quo, mod.n),
                                         n_modred_precomp(a_lo, one_precomp, mod.n),
                                         mod),
                               mod);

            if (witness != candidate)
            { printf("\n\n\nFAILURE 3 limbs!!\n\n\n"); return 0; }
        }

        if (bits < FLINT_BITS)
        { // 3 limbs, less direct
            NMOD_RED3(witness, a_hi_red, a_me, a_lo, mod);

            ulong one_precomp, W, W_pr_quo, W_pr_rem, W2, W2_pr;
            n_mulmod_precomp_shoup_quo_rem(&one_precomp, &W, 1L, mod.n);
            // --> W = remainder of division of 2**FLINT_BITS by n

            n_mulmod_precomp_shoup_quo_rem(&W_pr_quo, &W_pr_rem, W, mod.n);
            n_mulmod_and_precomp_shoup(&W2, &W2_pr, W, W, W_pr_quo, W_pr_rem, W_pr_quo, mod.n);
            // --> W2 = remainder of division of 2**(2*FLINT_BITS) by n

            // naive: res[i] = (vec[i] * W2 mod n + vec[i+1] * W mod n + vec[i+2] mod n) mod n
            ulong c0 = n_mulmod_shoup(W2, a_hi_red, W2_pr, mod.n);
            ulong c1 = n_mulmod_shoup(W, a_me, W_pr_quo, mod.n);
            ulong candidate = _nmod_add(c0, c1, mod);
            ulong correct = 0;
            if (a_lo >= (mod.n << (mod.norm-1)))
                correct = (mod.n << (mod.norm-1));
            candidate += a_lo - correct;
            umul_ppmm(c1, c0, one_precomp, candidate);
            candidate -= c1 * mod.n;
            if (candidate >= mod.n)
                candidate -= mod.n;

            if (witness != candidate)
            { printf("\n\n\nFAILURE 3 limbs less direct!!\n\n\n"); return 0; }
        }
    }

    FLINT_TEST_CLEAR(state);
    return 1;
}

SAMPLE(nmod_red)
SAMPLE(nmod_red_unroll)
SAMPLE(precomp)
SAMPLE(precomp_unroll)
SAMPLE(precomp2)
SAMPLE(precomp2_unroll)

SAMPLE(nmod2_red2)
SAMPLE(nmod2_red2_unroll)
SAMPLE(precomp22)
SAMPLE(precomp22_unroll)
SAMPLE(precomp22_bis)
SAMPLE(precomp22_ter)

SAMPLE(nmod_red3)
SAMPLE(nmod_red3_unroll)
SAMPLE(precomp3)
SAMPLE(precomp3_bis)
SAMPLE(precomp3_bis_unroll)
SAMPLE(precomp3_ter)

int main(void)
{
    double min[16], max[16];
    info_t info;
    flint_bitcnt_t i;

    //flint_printf("bits\tcheck\tred\txx_ur\tprcp\txx_ur\t||\t2red2\txx_ur\tprcp22\txx_ur\n");
    flint_printf("bits\tcheck\tred\txx_ur\tprcp\txx_ur\t||\t2red2\tprcp22\tp22_bis\t||\tred3\txx_ur\tprcp3\txx_ur\t||\tnew1\txx_ur\thorner3\n");
    for (i = 58; i <= FLINT_BITS; i+=1)
    {
        info.bits = i;

        int correct = check(i);

        prof_repeat(min+0, max+0, sample_nmod_red, (void *) &info);
        prof_repeat(min+1, max+1, sample_nmod_red_unroll, (void *) &info);
        prof_repeat(min+2, max+2, sample_precomp, (void *) &info);
        prof_repeat(min+3, max+3, sample_precomp_unroll, (void *) &info);
        prof_repeat(min+4, max+4, sample_nmod2_red2, (void *) &info);
        //prof_repeat(min+5, max+5, sample_nmod2_red2_unroll, (void *) &info);
        if (i == FLINT_BITS) min[6] = 0.; else prof_repeat(min+6, max+6, sample_precomp22, (void *) &info);
        if (i == FLINT_BITS) min[7] = 0.; else prof_repeat(min+7, max+7, sample_precomp22_bis, (void *) &info);
        //if (i == FLINT_BITS) min[7] = 0.; else prof_repeat(min+7, max+7, sample_precomp22_ter, (void *) &info);
        //if (i == FLINT_BITS) min[8] = 0.; else prof_repeat(min+8, max+8, sample_precomp22_unroll, (void *) &info);
        prof_repeat(min+9, max+9, sample_nmod_red3, (void *) &info);
        prof_repeat(min+10, max+10, sample_nmod_red3_unroll, (void *) &info);
        //if (i == FLINT_BITS) min[11] = 0.; else prof_repeat(min+11, max+11, sample_precomp3, (void *) &info);
        if (i == FLINT_BITS) min[11] = 0.; else prof_repeat(min+11, max+11, sample_precomp3_bis, (void *) &info);
        if (i == FLINT_BITS) min[12] = 0.; else prof_repeat(min+12, max+12, sample_precomp3_bis_unroll, (void *) &info);
        prof_repeat(min+13, max+13, sample_precomp2, (void *) &info);
        prof_repeat(min+14, max+14, sample_precomp2_unroll, (void *) &info);
        prof_repeat(min+15, max+15, sample_precomp3_ter, (void *) &info);

        flint_printf("%wd\t%s\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t||\t%.2lf\t%.2lf\t%.2lf\t||\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t||\t%.2lf\t%.2lf\t%.2lf\n",
                i,
                correct ? "pass" : "fail",
                (min[0]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[1]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[2]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[3]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[4]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                //(min[5]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[6]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[7]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                //(min[8]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[9]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[10]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[11]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[12]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[13]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[14]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN,
                (min[15]/(double)FLINT_CLOCK_SCALE_FACTOR)/LEN
                );
    }

    return 0;
}
