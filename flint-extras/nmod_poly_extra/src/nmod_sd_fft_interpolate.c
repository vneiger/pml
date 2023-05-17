#include <flint/fft_small.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

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

    nn = n;
    if (nn <= 1)
    {
        return;
    }

    half = (vec1d) (1 + ((ulong)(p) >> 1));

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
        
        for(i = 0; i < b; i++)
        {
            x[a+i] = vec1d_reduce_2n_to_n(vec1d_reduce_2n_to_n(tmp[b+i] + tmp[b+i], p) + x[a+i], p);
        }

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
        
        for (i = 0; i < b; i++)
        {
            vec1d u;
            u = tmp[i] + tmp[b+i];                 // [0..2p)
            u = x[a+i] - u;                        // (-2p..p)
            u = vec1d_mulmod(half, u, p, pinv);    // (-p..p)
            x[a+i] = (u < 0) ? u + p : u;
            x[i] = vec1d_reduce_2n_to_n(x[i] + x[a+i], p);
        }
        
        for (i = b; i < n; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[a+i], p, pinv);
            x[a+i] = (u < 0) ? u + p : u;
            x[i] = vec1d_reduce_2n_to_n(x[i] + x[a+i], p);
        }
        
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
void nmod_sd_tft_interpolate(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, const ulong N)
{
    ulong i, nn, k, aa, Nup;
    vec1d  *wk, *wk2;
    vec1d p, pinv, half, powerh, w;
    
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    Nup = FLINT_MAX(N, 4);
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

    half = (vec1d) (1 + ((ulong)(p) >> 1));
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
        _ifft(wk, Q, p, pinv, k);
        
        if (k >= 1)
        {
            w = -Q->w2tab[k][(1 << (k-1)) - 1];
            powerh = half;
            for (i = 1; i < k; i++)
            {
                powerh = vec1d_mulmod(powerh, half, p, pinv);
            }
        }
        else
        {
            w = 1;
            powerh = 1;
        }
        
        for (i = 0; i < aa; i++)
        {
            wk[i] = vec1d_mulmod(wk[i], powerh, p, pinv);
            powerh = vec1d_mulmod(w, powerh, p, pinv);                        
        }

        
        for (i = 0; i < aa; i++)
        {
            wk[i] = vec1d_reduce_2n_to_n(wk[i] + p, p);
        }

        
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
    CRT(wk, wk + Nup, N, p, pinv);

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
void nmod_sd_tft_evaluate_t(mp_ptr x, mp_srcptr A, sd_fft_lctx_t Qt, ulong N)
{
    ulong Nup, n, a, b, b2, t, i, j, k, ell, lambda;
    vec1d p, pinv, w, w2;
    vec1d *wk, *wk2;
    
    if (N == 0)
    {
        return;
    }
    
    if (N == 1)
    {
        x[0] = A[0];
        return;
    }

    Nup = FLINT_MAX(N, 4);
    wk = (vec1d *) aligned_alloc(32, 2 * Nup * sizeof(vec1d));
    
    for (i = 0; i < N; i++)
    {
        wk[i] = A[i];
    }
    for (; i < 2*Nup; i++)
    {
        wk[i] = 0;
    }

    p = Qt->p;
    pinv = Qt->pinv;
    
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
            w2 = -Qt->w2tab[k][(1 << (k-1)) - 1];
            w = 1;
            for (i = 0; i < a; i++)
            {
                vec1d z;
                z = vec1d_mulmod(wk[i], w, p, pinv);                  // (-p..p)
                wk[i] = (z < 0) ? z + p : z;
                w = vec1d_mulmod(w, w2, p, pinv);                     // w2 still in (-p..p)
            }
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
                for (i = 1; i < lambda; i++)
                {
                    vec1d u, v, w;
                    u = wk[N-n+t];
                    v = wk[N-n+a];
                    w = v + u;
                    v = v - u;
                    w = (w > p) ? w - p : w;
                    v = (v < 0) ? v + p : v;
                    wk[N-n+t] = w;
                    wk[N-n+a+t] = v;
                    t++;
                    u = wk[N-n+t];
                    v = wk[N-n+a+1];
                    w = v + u;
                    v = v - u;
                    w = (w > p) ? w - p : w;
                    v = (v < 0) ? v + p : v;
                    wk[N-n+t] = w;
                    wk[N-n+a+t] = v;
                    t++;
                }
                
                vec1d u, v, w;
                u = wk[N-n];
                v = wk[N-n+a];
                w = v + u;
                v = v - u;
                w = (w > p) ? w - p : w;
                v = (v < 0) ? v + p : v;
                wk[N-n] = w;
                wk[N-n+a] = v;
                u = wk[N-n+1];
                v = wk[N-n+a+1];
                w = v + u;
                v = v - u;
                w = (w > p) ? w - p : w;
                v = (v < 0) ? v + p : v;
                wk[N-n+1] = w;
                wk[N-n+a+1] = v;
                break;
                
            case 2:
            default:
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j++)
                    {
                        vec1d u, v, w;
                        u = wk[N-n+t];
                        v = wk[N-n+a+j];
                        w = v + u;
                        v = v - u;
                        w = (w > p) ? w - p : w;
                        v = (v < 0) ? v + p : v;
                        wk[N-n+t] = w;
                        wk[N-n+a+t] = v;
                        t++;
                    }
                }
                for (j = 0; j < b2; j++)
                {
                    vec1d u, v, w;
                    u = wk[N-n+j];
                    v = wk[N-n+a+j];
                    w = v + u;
                    v = v - u;
                    w = (w > p) ? w - p : w;
                    v = (v < 0) ? v + p : v;
                    wk[N-n+j] = w;
                    wk[N-n+a+j] = v;
                }
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = wk[i];
    }

    flint_free(wk);
}


