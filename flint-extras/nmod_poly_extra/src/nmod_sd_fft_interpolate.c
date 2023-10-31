#include <flint/fft_small.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

#define TRY_AVX2

FLINT_FORCE_INLINE vec1d vec1d_addmod_limited(vec1d a, vec1d b, vec1d p)
{
    vec1d s = a + b;
    if (s > p)
        return s - p;
    else
        return s;
}

/*------------------------------------------------------------*/
/* base cases are functions from sd_ifft.c                    */
/* note: they are not static in the flint source tree         */
/*       but the corresponding ones in sd_fft.c are           */
/* k=4 is a wrapper for a forced inline function              */
/*------------------------------------------------------------*/
void do_sd_ifft_basecase_4_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_ifft_basecase_5_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_ifft_basecase_6_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_ifft_basecase_7_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_ifft_basecase_8_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_ifft(sd_fft_lctx_t Q, double *X, long k); /* use for k >= 8 */


/*------------------------------------------------------------*/
/*  in-place size 1 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void sd_ifft_1(vec1d * x)
{
    vec1d u0, u1;
    
    u0 = x[0];
    u1 = x[1];
    x[0] = u0 + u1;
    x[1] = u0 - u1;
}

/*------------------------------------------------------------*/
/*  in-place size 2 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void sd_ifft_2(vec1d *x, const vec1d w, const vec1d p, const vec1d pinv)
{
    vec1d u0, u1, u2, u3, v0, v1, v2, v3;
    
    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];

    v0 = u0 + u1;
    v1 = u0 - u1;
    v2 = u2 + u3;
    v3 = u2 - u3;
    
    x[0] = v0 + v2;
    x[2] = v0 - v2;
    v3 = vec1d_mulmod(v3, w, p, pinv);      
    x[1] = v1 - v3;
    x[3] = v1 + v3;
}

/*------------------------------------------------------------*/
/* in-place inverse FFT in size 2^3                           */
/*------------------------------------------------------------*/
static void sd_ifft_3(vec1d *x, const vec1d t1, const vec1d t2, const vec1d t3, const vec1d p, const vec1d pinv)
{
    vec1d u0, u1, u2, u3, v0, v1, v2, v3;
    
    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];
    
    v0 = u0 + u1;
    v2 = u2 + u3;
    v1 = u0 - u1;
    v3 = vec1d_mulmod(u3 - u2, t1, p, pinv);
    
    x[0] = v0 + v2; 
    x[1] = v1 + v3; 
    x[2] = v0 - v2; 
    x[3] = v1 - v3; 

    u0 = x[0 + 4];
    u1 = x[1 + 4];
    u2 = x[2 + 4];
    u3 = x[3 + 4];
    
    v0 = u0 + u1;
    v2 = u2 + u3;
    v1 = u0 - u1;
    v3 = vec1d_mulmod(u3 - u2, t1, p, pinv);
    
    x[0 + 4] = v0 + v2; 
    x[1 + 4] = v1 + v3; 
    x[2 + 4] = v0 - v2; 
    x[3 + 4] = v1 - v3; 
    
    u0 = x[0];
    u1 = x[4];
    x[0] = u0 + u1;
    x[4] = u0 - u1;
    u0 = x[1];
    u1 = vec1d_mulmod(x[1 + 4], t2, p, pinv);
    x[1] = u0 - u1;
    x[1 + 4] = u0 + u1;
    u0 = x[2];
    u1 = vec1d_mulmod(x[2 + 4], t1, p, pinv);
    x[2] = u0 - u1;
    x[2 + 4] = u0 + u1;
    u0 = x[3];
    u1 = vec1d_mulmod(x[3 + 4], t3, p, pinv);
    x[3] = u0 - u1;
    x[3 + 4] = u0 + u1;
}

/*---------------------------------------------------------*/
/* reduces s mod (X^sz-1), assuming s has length N         */
/* assumes sz divides N                                    */
/* input in [0..p), output in [0..p)                       */
/*---------------------------------------------------------*/
static inline
void fold_minus(vec1d *x, vec1d *s, ulong sz, ulong N, vec1d p)
{
    ulong i, j;
    // we don't need sz = 1
    switch(sz)
    {
        case 2:
            x[0] = s[0];
            x[1] = s[1];
            i = 2;
            while (i < N)
            {
                x[0] = vec1d_reduce_2n_to_n(x[0] + s[i], p);
                x[1] = vec1d_reduce_2n_to_n(x[1] + s[i+1], p);
                i+=2;
            }
            break;
        default:
            for (i = 0; i < sz; i++)
            {
                x[i] = s[i];
            }
            while (i < N)
            {
                for (j = 0; j < sz; j++, i++)
                {
                    x[j] = vec1d_reduce_2n_to_n(x[j] + s[i], p);
                }
            }
            break;
    }
}

/*-----------------------------------------------------------*/
/* applies the CRT map to x in place                         */
/* x has length n                                            */
/* input is in [0,p), output is in [0,p)                     */
/* tmp is a temporary workspace, of length at least 2n       */
/*-----------------------------------------------------------*/
static inline
void CRT(vec1d *x, vec1d *tmp, ulong n, vec1d p, vec1d pinv)
{
    ulong a, b, b2, n2, i, nn;
    vec1d half;
    vec4d p4, pinv4, half4;
    
    nn = n;
    if (nn <= 1)
    {
        return;
    }

    half = (vec1d) (1 + ((ulong)(p) >> 1));
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);
    half4 = vec4d_set_d(half);

    a = 1L;
    n2 = n >> 1;
    while (a <= n2)
    {
        a <<= 1;
    }
    
    while (n != a)
    {
        n = n - a;
        b = 1L;
        n2 = (n >> 1);
        while (b <= n2)
        {
            b <<= 1;
        }
        b2 = b << 1;

        fold_minus(tmp, x, b2, a, p);

#ifdef TRY_AVX2
        i = 0;
        if (b >= 4)
            for(; i < b; i += 4)
            {
                vec4d u;
                u = vec4d_load_aligned(tmp + b + i);
                vec4d_store_aligned(x + a + i, vec4d_reduce_2n_to_n(vec4d_add(vec4d_reduce_2n_to_n(vec4d_add(u, u), p4), vec4d_load_aligned(x + a + i)), p4));
            }
        for(; i < b; i++)
        {
            x[a+i] = vec1d_reduce_2n_to_n(vec1d_reduce_2n_to_n(tmp[b+i] + tmp[b+i], p) + x[a+i], p);
        }
#else
        for(i = 0; i < b; i++)
        {
            x[a+i] = vec1d_reduce_2n_to_n(vec1d_reduce_2n_to_n(tmp[b+i] + tmp[b+i], p) + x[a+i], p);
        }
#endif
        tmp += b2;
        x += a;
        a = b;
    }

    while (n != nn)
    {
        b = a;
        a <<= 1;
        b2 = a;
        while (!(a & nn))
        {
            a <<= 1;
        }
        tmp -= b2;
        x -= a;

#ifdef TRY_AVX2
        i = 0;
        if (b >= 4)
            for (; i < b; i += 4)
            {
                vec4d u, x4;
                u = vec4d_add(vec4d_load(tmp + i), vec4d_load(tmp + b + i));
                u = vec4d_sub(vec4d_load(x + a + i), u);
                u = vec4d_mulmod(half4, u, p4, pinv4);
                x4 = vec4d_reduce_2n_to_n(vec4d_add(u, p4), p4);
                vec4d_store_aligned(x + a + i, x4);
                vec4d_store_aligned(x + i, vec4d_reduce_2n_to_n(vec4d_add(vec4d_load(x + i), x4), p4));
            }
        for (; i < b; i++)
        {
            vec1d u;
            u = tmp[i] + tmp[b+i];                 // [0..2p)
            u = x[a+i] - u;                        // (-2p..p)
            u = vec1d_mulmod(half, u, p, pinv);    // (-p..p)
            x[a+i] = vec1d_reduce_2n_to_n(u + p, p); // [0..p)
            x[i] = vec1d_reduce_2n_to_n(x[i] + x[a+i], p); // [0..p)
        }
#else
        for (i = 0; i < b; i++)
        {
            vec1d u;
            u = tmp[i] + tmp[b+i];                 // [0..2p)
            u = x[a+i] - u;                        // (-2p..p)
            u = vec1d_mulmod(half, u, p, pinv);    // (-p..p)
            x[a+i] = vec1d_reduce_2n_to_n(u + p, p); // [0..p)
            x[i] = vec1d_reduce_2n_to_n(x[i] + x[a+i], p); // [0..p)
        }
#endif

#ifdef TRY_AVX2
        i = b;
        if (b >= 4)
            for (; i < n; i += 4)
            {
                vec4d u, x4;
                u = vec4d_mulmod(half4, vec4d_load_aligned(x + a + i), p4, pinv4);  // (-p..p)
                x4 = vec4d_reduce_2n_to_n(vec4d_add(u, p4), p4);
                vec4d_store_aligned(x + a + i, x4);
                vec4d_store_aligned(x + i, vec4d_reduce_2n_to_n(vec4d_add(vec4d_load_aligned(x + i), x4), p4));
            }
        for (; i < n; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[a+i], p, pinv);  // (-p..p)
            x[a+i] = vec1d_reduce_2n_to_n(u + p, p);  // [0..p)
            x[i] = vec1d_reduce_2n_to_n(x[i] + x[a+i], p); // [0..p)
        }
#else
        for (i = b; i < n; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[a+i], p, pinv);  // (-p..p)
            x[a+i] = vec1d_reduce_2n_to_n(u + p, p);  // [0..p)
            x[i] = vec1d_reduce_2n_to_n(x[i] + x[a+i], p); // [0..p)
        }
#endif   
        n = n + a;
    }
}

/*------------------------------------------------------------*/
/* main inverse fft routine, in place                         */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE void _ifft(vec1d *wk, sd_fft_lctx_t Q, vec1d p, vec1d pinv, ulong k)
{
    switch(k)
    {
        case 0:
            break;
        case 1:
            sd_ifft_1(wk);
            break;
        case 2:
            sd_ifft_2(wk, Q->w2tab[1][0], p, pinv);
            break;
        case 3:
            sd_ifft_3(wk, Q->w2tab[1][0], Q->w2tab[2][1], Q->w2tab[2][0], p, pinv);
            break;
        case 4:
            do_sd_ifft_basecase_4_1(Q, wk, 0, 0);
            break;
        case 5:
            sd_ifft_basecase_5_1(Q, wk, 0, 0);
            break;
        case 6:
            sd_ifft_basecase_6_1(Q, wk, 0, 0);
            break;
        case 7:
            sd_ifft_basecase_7_1(Q, wk, 0, 0);
            break;
        case 8:
            sd_ifft_basecase_8_1(Q, wk, 0, 0);
            break;
        default:
            sd_ifft(Q, wk, k);
            break;
    }
}                       

/*------------------------------------------------------------*/
/* fft interpolation                                          */
/*------------------------------------------------------------*/
void nmod_sd_fft_interpolate(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, const ulong k)
{
    ulong i, N;
    vec1d p, pinv;
    vec1d *dx;
    
    N = 1L << k;
    p = Q->p;
    pinv = Q->pinv;

    dx = (vec1d *) aligned_alloc(32, FLINT_MAX(N, 4) * sizeof(vec1d));
    
    for (i = 0; i < N; i++)
    {
        dx[i] = (vec1d) x[i];
    }
 
    _ifft(dx, Q, p, pinv, k);

    nmod_poly_fit_length(poly, 1L << k);
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = (mp_limb_t) vec1d_reduce_to_0n(dx[i], p, pinv);
    }
    _nmod_poly_normalise(poly);
    flint_free(dx);
}

    
/*------------------------------------------------------------*/
/* transposed fft evaluation                                  */
/* obtained by using the inverse root and ifft routines       */
/*------------------------------------------------------------*/
void nmod_sd_fft_evaluate_t(mp_ptr x, const nmod_poly_t poly, sd_fft_lctx_t Qt, const ulong k)
{
    ulong i, N, d;
    vec1d *dx;
    vec1d p, pinv;
    
    N = 1L << k;
    d = poly->length - 1;
    p = Qt->p;
    pinv = Qt->pinv;

    dx = (vec1d *) aligned_alloc(32, FLINT_MAX(N, 4) * sizeof(vec1d));
    for (i = 0; i <= d; i++)
    {
        dx[i] = (vec1d) poly->coeffs[i];
    }
    for (; i < N; i++)
    {
        dx[i] = 0;
    }

    _ifft(dx, Qt, p, pinv, k);
    for (i = 0; i < N; i++)
    {
        x[i] = (mp_limb_t) vec1d_reduce_to_0n(dx[i], p, pinv);
    }
    flint_free(dx);
}

/*------------------------------------------------------------*/
/* tft interpolation                                          */
/*------------------------------------------------------------*/
void nmod_sd_tft_interpolate(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, nmod_sd_fft_t F, const ulong N)
{
    ulong i, nn, k, aa, Nup;
    vec1d  *wk, *wk2, *powers_w;
    vec1d p, pinv;
    vec4d p4, pinv4;
    
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    Nup = FLINT_MAX(4, n_next_pow2m1(N-1) + 1);
    wk = (vec1d *) aligned_alloc(32, 3 * Nup * sizeof(vec1d));
    
    for (i = 0; i < N; i++)
    {
        wk[i] = x[i];
    }
    for (; i < 3 * Nup; i++)
    {
        wk[i] = 0;
    }

    wk2 = wk;
    p = Q->p;
    pinv = Q->pinv;
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);

    nn = N;
    k = 0;
    aa = 1;
    
    while (aa <= nn/2)
    {
        k++;
        aa <<= 1;
    }

    do
    {
        _ifft(wk, Q, p, pinv, k); // (-3p..3p), expectedly
        powers_w = F->powers_inv_w_over_2[k]; // (-p/2..p/2)

#ifdef TRY_AVX2
        i = 0;
        if (aa >= 4)
            for (; i < aa; i += 4)
            {
                vec4d z;
                z = vec4d_mulmod(vec4d_load_aligned(wk + i),
                                 vec4d_load_aligned(powers_w + i), p4, pinv4);                  // (-p..p)
                vec4d_store_aligned(wk + i, vec4d_reduce_2n_to_n(vec4d_add(z, p4), p4));  // [0..p)
            }
        for (; i < aa; i++)
        {
            vec1d z = vec1d_mulmod(wk[i], powers_w[i], p, pinv);  // (-p..p)
            wk[i] = vec1d_reduce_2n_to_n(z + p, p);               // [0..p)
        }
#else
        for (i = 0; i < aa; i++)
        {
            vec1d z = vec1d_mulmod(wk[i], powers_w[i], p, pinv);  // (-p..p)
            wk[i] = vec1d_reduce_2n_to_n(z + p, p);               // [0..p)
        }
#endif
        wk += aa;
        nn = nn-aa;
        
        aa = 1;
        k = 0;
        while (aa <= nn/2)
        {
            k++;
            aa <<= 1;
        }
    }
    while (nn != 0);

    wk = wk2;
    CRT(wk, wk + Nup, N, p, pinv); // [0..p)

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = wk[i];
    }
    
    flint_free(wk);
    _nmod_poly_normalise(poly);
}

/*-----------------------------------------------------------*/
/* transpose tft in length N                                 */
/*-----------------------------------------------------------*/
void nmod_sd_tft_evaluate_t(mp_ptr x, mp_srcptr A, sd_fft_lctx_t Qt, nmod_sd_fft_t F, ulong N)
{
    ulong Nup, n, a, b, b2, t, i, j, k, ell, lambda;
    vec1d p, pinv;
    vec4d p4, pinv4;
    vec1d *wk, *wk2, *powers_w;
    
    if (N == 0)
    {
        return;
    }
    
    if (N == 1)
    {
        x[0] = A[0];
        return;
    }

    Nup = 2 * FLINT_MAX(4, n_next_pow2m1(N-1) + 1);
    wk = (vec1d *) aligned_alloc(32, Nup * sizeof(vec1d));
    
    for (i = 0; i < N; i++)
    {
        wk[i] = A[i];
    }
    for (; i < Nup; i++)
    {
        wk[i] = 0;
    }

    p = Qt->p;
    pinv = Qt->pinv;
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);
    
    k = 0;
    a = 1;
    n = N;
    
    while (a <= n/2)
    {
        k++;
        a <<= 1;
    }

    wk2 = wk;
    do
    {

        _ifft(wk, Qt, p, pinv, k);

        if (k > 0)
        {
            powers_w = F->powers_w[k+1];
#ifdef TRY_AVX2
            i = 0;
            if (a >= 4) // should remove this test
                for (; i < a; i += 4)
                {
                    vec4d z;
                    z = vec4d_mulmod(vec4d_load_aligned(wk + i),
                                     vec4d_load_aligned(powers_w + i), p4, pinv4);                  // (-p..p)
                    vec4d_store_aligned(wk + i, vec4d_reduce_2n_to_n(vec4d_add(z, p4), p4));  // [0..p)
                }
            for (; i < a; i++)
            {
                vec1d z;
                z = vec1d_mulmod(wk[i], powers_w[i], p, pinv);                  // (-p..p)
                wk[i] = vec1d_reduce_2n_to_n(z + p, p);  // [0..p)
            }
#else
            for (i = 0; i < a; i++)
            {
                vec1d z;
                z = vec1d_mulmod(wk[i], powers_w[i], p, pinv);                  // (-p..p)
                wk[i] = vec1d_reduce_2n_to_n(z + p, p);  // [0..p)
            }
#endif

        }
        wk += a;
        n = n-a;
       
        a = 1;
        k = 0;
        while (a <= n/2)
        {
            k++;
            a <<= 1;
        }

    }
    while (n != 0);
    wk = wk2;



    a = 1;
    k = 0;
    while (! (a & N))
    {
        k++;
        a <<= 1;
    }
    
    for (i = 0; i < a; i++)
    {
        if (wk[N-a+i] != 0)
        {
            wk[N+i] = p - wk[N-a+i];
        }
    }
    
    n = a;
    while (n != N)
    {
        b = a;
        ell = k;
        a <<= 1;
        k++;
        while (!(a & N))
        {
            a <<= 1;
            k++;
        }
        
        b2 = b << 1;
        lambda = 1L << (k-ell-1);

        n = n + a;
        t = b2;
        // we have to deal with i=0 separately
        // (since otherwise we would erase the source)
        switch(b)
        {
            case 1:
// not using AVX2 here
#ifdef TRY_AVX2
                i = 1;
                for (; i < lambda; i++)
                {
                    vec1d u, v, w;
                    u = wk[N-n+t];
                    v = wk[N-n+a];
                    w = v + u;
                    v = v - u;
                    w = vec1d_reduce_2n_to_n(w, p);
                    v = vec1d_reduce_2n_to_n(v + p, p);
                    wk[N-n+t] = w;
                    wk[N-n+a+t] = v;
                    t++;
                    u = wk[N-n+t];
                    v = wk[N-n+a+1];
                    w = v + u;
                    v = v - u;
                    w = vec1d_reduce_2n_to_n(w, p);
                    v = vec1d_reduce_2n_to_n(v + p, p);
                    wk[N-n+t] = w;
                    wk[N-n+a+t] = v;
                    t++;
                }
#else
                for (i = 1; i < lambda; i++)
                {
                    vec1d u, v, w;
                    u = wk[N-n+t];
                    v = wk[N-n+a];
                    w = v + u;
                    v = v - u;
                    w = vec1d_reduce_2n_to_n(w, p);
                    v = vec1d_reduce_2n_to_n(v + p, p);
                    wk[N-n+t] = w;
                    wk[N-n+a+t] = v;
                    t++;
                    u = wk[N-n+t];
                    v = wk[N-n+a+1];
                    w = v + u;
                    v = v - u;
                    w = vec1d_reduce_2n_to_n(w, p);
                    v = vec1d_reduce_2n_to_n(v + p, p);
                    wk[N-n+t] = w;
                    wk[N-n+a+t] = v;
                    t++;
                }
#endif
                 
                vec1d u, v, w;
                u = wk[N-n];
                v = wk[N-n+a];
                w = v + u;
                v = v - u;
                w = vec1d_reduce_2n_to_n(w, p);
                v = vec1d_reduce_2n_to_n(v + p, p);
                wk[N-n] = w;
                wk[N-n+a] = v;
                u = wk[N-n+1];
                v = wk[N-n+a+1];
                w = v + u;
                v = v - u;
                w = vec1d_reduce_2n_to_n(w, p);
                v = vec1d_reduce_2n_to_n(v + p, p);
                wk[N-n+1] = w;
                wk[N-n+a+1] = v;
                break;
                
            case 2:
            default:
                for (i = 1; i < lambda; i++)
                {
#ifdef TRY_AVX2
                    j = 0;
                    if (b2 >= 4)
                        for (; j < b2; j += 4)
                        {
                            vec4d u, v, w;
                            u = vec4d_load_aligned(wk + N - n + t);
                            v = vec4d_load_aligned(wk + N - n + a + j);
                            w = vec4d_add(v, u);
                            v = vec4d_sub(v, u);
                            w = vec4d_reduce_2n_to_n(w, p4);
                            v = vec4d_reduce_2n_to_n(vec4d_add(v, p4), p4);
                            vec4d_store_aligned(wk + N - n + t, w);
                            vec4d_store_aligned(wk + N - n +  a + t, v);
                            t += 4;
                    }
                    for (; j < b2; j++)
                    {
                        vec1d u, v, w;
                        u = wk[N-n+t];
                        v = wk[N-n+a+j];
                        w = v + u;
                        v = v - u;
                        w = vec1d_reduce_2n_to_n(w, p);
                        v = vec1d_reduce_2n_to_n(v + p, p);
                        wk[N-n+t] = w;
                        wk[N-n+a+t] = v;
                        t++;
                    }
#else
                    for (j = 0; j < b2; j++)
                    {
                        vec1d u, v, w;
                        u = wk[N-n+t];
                        v = wk[N-n+a+j];
                        w = v + u;
                        v = v - u;
                        w = vec1d_reduce_2n_to_n(w, p);
                        v = vec1d_reduce_2n_to_n(v + p, p);
                        wk[N-n+t] = w;
                        wk[N-n+a+t] = v;
                        t++;
                    }
#endif
                }
                

#ifdef TRY_AVX2
                    j = 0;
                    if (b2 >= 4)
                        for (; j < b2; j += 4)
                        {
                            vec4d u, v, w;
                            u = vec4d_load_aligned(wk + N - n + j);
                            v = vec4d_load_aligned(wk + N - n + a + j);
                            w = vec4d_add(v, u);
                            v = vec4d_sub(v, u);
                            w = vec4d_reduce_2n_to_n(w, p4);
                            v = vec4d_reduce_2n_to_n(vec4d_add(v, p4), p4);
                            vec4d_store_aligned(wk + N - n + j, w);
                            vec4d_store_aligned(wk + N - n +  a + j, v);
                    }
                    for (; j < b2; j++)
                    {
                        vec1d u, v, w;
                        u = wk[N-n+j];
                        v = wk[N-n+a+j];
                        w = v + u;
                        v = v - u;
                        w = vec1d_reduce_2n_to_n(w, p);
                        v = vec1d_reduce_2n_to_n(v + p, p);
                        wk[N-n+j] = w;
                        wk[N-n+a+j] = v;
                    }
#else
                for (j = 0; j < b2; j++)
                {
                    vec1d u, v, w;
                    u = wk[N-n+j];
                    v = wk[N-n+a+j];
                    w = v + u;
                    v = v - u;
                    w = vec1d_reduce_2n_to_n(w, p);
                    v = vec1d_reduce_2n_to_n(v + p, p);
                    wk[N-n+j] = w;
                    wk[N-n+a+j] = v;
                }
#endif
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = wk[i];
    }

    flint_free(wk);
}


