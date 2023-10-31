#include <flint/fft_small.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

#define TRY_AVX2

/*------------------------------------------------------------*/
/* base cases are functions from sd_ifft.c                    */
/* k=4 is a wrapper for a forced inline function              */
/* inputs in (-3p..3p), outputs in (-3p..3p)                  */
/*------------------------------------------------------------*/
void do_sd_fft_basecase_4_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_fft_basecase_5_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_fft_basecase_6_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_fft_basecase_7_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_fft_basecase_8_1(const sd_fft_lctx_t Q, double* X, ulong j_r, ulong j_bits);
void sd_fft(sd_fft_lctx_t Q, double *X, long k); /* use for k >= 8 */

/*------------------------------------------------------------*/
/*  in-place size 2^1 butterfly                               */
/*  o1, o2 = i1+i2, i1-i2                                     */
/*  input in (-p..p), output in (-2p..2p)                     */
/*------------------------------------------------------------*/
static inline void _fft_1(vec1d *x)
{
    vec1d u0, u1;

    u0 = x[0];
    u1 = x[1];
    
    x[0] = u0 + u1;
    x[1] = u0 - u1;
}

/*------------------------------------------------------------*/
/*  in-place size 2^2 butterfly                               */
/*  t1, t2, t3, t4 = i1+i3, i2+i4, i1-i3, w(i2-i4)            */
/*  o1, o2, o3, o4 = t1+t2, t1-t2, t3+t4, t3-t4               */
/*  input in (-p..p), output in (-4p..4p)                     */
/*------------------------------------------------------------*/
static inline void _fft_2(vec1d *x, const vec1d w, const vec1d p, const vec1d pinv)
{
    vec1d u0, u1, u2, u3, v0, v1, v2, v3;

    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];

    v0 = u0 + u2;  // (-2p..2p)
    v2 = u0 - u2;  // (-2p..2p)
    v1 = u1 + u3;  // (-2p..2p)
    v3 = vec1d_mulmod(u1 - u3, w, p, pinv);      // (-p..p)

    x[0] = v0 + v1;  // (-4p..4p) 
    x[1] = v0 - v1;  // (-4p..4p)
    x[2] = v2 + v3;  // (-3p..3p)
    x[3] = v2 - v3;  // (-3p..3p)
}


/*------------------------------------------------------------*/
/* in-place FFT in size 2^3                                   */
/* input in (-p..p), output in (-4p..4p)                      */
/*------------------------------------------------------------*/
static void _fft_3(vec1d *x, const vec1d t1, const vec1d t2, const vec1d t3, const vec1d p, const vec1d pinv)
{
    vec4d w4, p4, pinv4, u0, u1, u, v;
    vec1d z0, z1, z2, z3, v0, v1, v2, v3;

    w4 = vec4d_set_d4(1, t2, t1, t3);
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);
    
    u0 = vec4d_load(x);
    u1 = vec4d_load(x + 4);
    
    u = vec4d_add(u0, u1); // (-2p..2p)
    z0 = u[0];
    z1 = u[1];
    z2 = u[2];
    z3 = u[3];
    v0 = z0 + z2;  // (-4p..4p)
    v2 = z0 - z2;  // (-4p..4p)
    v1 = z1 + z3;  // (-4p..4p)
    v3 = vec1d_mulmod(z1 - z3, t1, p, pinv);      // (-p..p)
    x[0] = vec1d_reduce_to_pm1no(v0 + v1, p, pinv);  // (-p..p)
    x[1] = vec1d_reduce_to_pm1no(v0 - v1, p, pinv);      // (-p..p)
    x[2] = vec1d_reduce_to_pm1no(v2 + v3, p, pinv);  // (-p..p)
    x[3] = vec1d_reduce_to_pm1no(v2 - v3, p, pinv);      // (-p..p)
    
    u = vec4d_sub(u0, u1); // (-2p..2p)
    v = vec4d_mulmod(u, w4, p4, pinv4); // (-p..p)
    z0 = v[0];
    z1 = v[1];
    z2 = v[2];
    z3 = v[3];
    v0 = z0 + z2;  // (-2p..2p)
    v2 = z0 - z2;  // (-2p..2p)
    v1 = z1 + z3;  // (-2p..2p)
    v3 = vec1d_mulmod(z1 - z3, t1, p, pinv);      // (-p..p)
    x[4] = v0 + v1;  // (-4p..4p)
    x[5] = v0 - v1;  // (-4p..4p)
    x[6] = v2 + v3;  // (-3p..3p)
    x[7] = v2 - v3;  // (-3p..3p)
}

/*------------------------------------------------------------*/
/* main fft routine, in place                                 */
/* input in (-p..p), output in (-3p..3p)                      */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE void _fft(vec1d *wk, sd_fft_lctx_t Q, vec1d p, vec1d pinv, ulong k)
{
    switch(k)
    {
        case 0:
            break;
        case 1:
            _fft_1(wk);
            break;
        case 2:
            _fft_2(wk, Q->w2tab[1][0], p, pinv);
            break;
        case 3:
            _fft_3(wk, Q->w2tab[1][0], Q->w2tab[2][0], Q->w2tab[2][1], p, pinv);
            break;
        case 4:
            do_sd_fft_basecase_4_1(Q, wk, 0, 0);
            break;
        case 5:
            sd_fft_basecase_5_1(Q, wk, 0, 0);
            break;
        case 6:
            sd_fft_basecase_6_1(Q, wk, 0, 0);
            break;
        case 7:
            sd_fft_basecase_7_1(Q, wk, 0, 0);
            break;
        case 8:
            sd_fft_basecase_8_1(Q, wk, 0, 0);
            break;
        default:
            sd_fft(Q, wk, k);
            break;
    }
}                       

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^i), n=2^k, in some order             */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_sd_fft_evaluate(mp_ptr x, const nmod_poly_t poly, sd_fft_lctx_t Q, const ulong k)
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
        dx[i] = (vec1d) nmod_poly_get_coeff_ui(poly, i);
    }

    _fft(dx, Q, p, pinv, k);

    for (i = 0; i < N; i++)
    {
        x[i] = (mp_limb_t) vec1d_reduce_to_0n(dx[i], p, pinv);
    }
    flint_free(dx);
}

/*------------------------------------------------------------*/
/* transposed fft interpolation                               */
/* obtained by using the inverse root and fft routines        */
/*------------------------------------------------------------*/
void nmod_sd_fft_interpolate_t(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, const ulong k)
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

    _fft(dx, Q, p, pinv, k);

    nmod_poly_fit_length(poly, 1L << k);
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = (mp_limb_t) vec1d_reduce_to_0n(dx[i], p, pinv);
    }
    _nmod_poly_normalise(poly);
    flint_free(dx);
}

/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_sd_tft_evaluate(mp_ptr x, const nmod_poly_t poly, sd_fft_lctx_t Q, nmod_sd_fft_t F, const ulong N)
{
    ulong i, j, k, a, ell, b, t, b2, lambda, N2;
    vec1d *x2, *x2_bak, *powers_rho;
    vec1d p, pinv;
    vec4d p4, pinv4;

    
    // N <= 1 -> easy case
    if (N <= 1)
    {
        for (i = 0; i < N; i++)
        {
            x[i] = nmod_poly_get_coeff_ui(poly, i);
        }
        return;
    }

    // else, loop through the reduction process
    p = Q->p;
    pinv = Q->pinv;
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);
    
    N2 = FLINT_MAX(4, n_next_pow2m1(2*N-1) + 1);
    x2 = (vec1d *) aligned_alloc(32, N2 * sizeof(vec1d));
    x2_bak = x2;

    // x2 in [0..p)
    for (i = 0; i < N; i++)
    {
        x2[i] = nmod_poly_get_coeff_ui(poly, i);
    }
    for (; i < N2; i++)
    {
        x2[i] = 0;
    }
    
    ell = 0;
    b = 1;
    while (b <= N/2)
    {
        ell++;
        b <<= 1;
    }

    do
    {
        a = b;
        k = ell;
        
        b = a >> 1;
        ell = k-1;
        while (b > 0 && (! (N & b)))
        {
            b >>= 1;
            ell--;
        }

        t = 0;
        b2 = b << 1;

        // assume x2 in [0..p)
        powers_rho = F->powers_w[k+1];

        switch(b)
        {
            case 0:
#ifdef TRY_AVX2
                i = 0;
                if (a >= 4)
                    for (; i < a; i += 4)
                    {
                        vec4d_store_aligned(x2 + i, vec4d_mulmod(
                                                vec4d_sub(vec4d_load_aligned(x2 + i),
                                                          vec4d_load_aligned(x2 + i + a)),
                                                vec4d_load_aligned(powers_rho + i), p4, pinv4));
                    }
                for (; i < a; i++)
                {
                    x2[i] = vec1d_mulmod(vec1d_sub(x2[i], x2[i+a]), powers_rho[i], p, pinv); // (-p..p)
                }
#else
                for (i = 0; i < a; i++)
                {
                    x2[i] = vec1d_mulmod(vec1d_sub(x2[i], x2[i+a]), powers_rho[i], p, pinv); // (-p..p)
                }
#endif
                break;

            default:
                // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
                // i = 0 : special case
#ifdef TRY_AVX2
                if (b2 >= 4)
                    for (; t < b2; t += 4)
                    {
                        vec4d u, v;
                        u = vec4d_load_aligned(x2 + t);
                        v = vec4d_load_aligned(x2 + a + t);
                        vec4d_store_aligned(x2 + t, vec4d_mulmod(vec4d_sub(u, v),
                                                                 vec4d_load_aligned(powers_rho + t), p4, pinv4));
                        vec4d_store_aligned(x2 + a + t, vec4d_reduce_2n_to_n(vec4d_add(u, v), p4));
                    }
                for (; t < b2; t++)
                {
                    vec1d u, v;
                    u = x2[t]; 
                    v = x2[a+t]; 
                    x2[t] = vec1d_mulmod(vec1d_sub(u, v), powers_rho[t], p, pinv); // (-p..p)
                    x2[a+t] = vec1d_reduce_2n_to_n(u+v, p); // [0..p)
                }
#else
                for (; t < b2; t++)
                {
                    vec1d u, v;
                    u = x2[t]; 
                    v = x2[a+t]; 
                    x2[t] = vec1d_mulmod(vec1d_sub(u, v), powers_rho[t], p, pinv); // (-p..p)
                    x2[a+t] = vec1d_reduce_2n_to_n(u+v, p); // [0..p)
                }
#endif
                lambda = (1L << k) >> (ell+1);
                for (i = 1; i < lambda; i++)
                {
#ifdef TRY_AVX2
                    j = 0;
                    if (b2 >= 4)
                        for (j = 0; j < b2; j += 4)
                        {
                            vec4d u, v, tmp;
                            u = vec4d_load_aligned(x2 + t);
                            v = vec4d_load_aligned(x2 + a + t);
                            vec4d_store_aligned(x2 + t, vec4d_mulmod(vec4d_sub(u, v), vec4d_load_aligned(powers_rho + t),
                                                                     p4, pinv4));
                            tmp = vec4d_reduce_2n_to_n(vec4d_add(u, v), p4); // [0..p)
                            vec4d_store_aligned(x2 + a + j, vec4d_reduce_2n_to_n(
                                                    vec4d_add(tmp, vec4d_load_aligned(x2 + a + j)), p4)); // [0..p)
                            t += 4;
                        }
                    for (; j < b2; j++)
                    {
                        vec1d u, v, tmp;
                        u = x2[t];
                        v = x2[a+t];
                        x2[t] = vec1d_mulmod(vec1d_sub(u, v), powers_rho[t], p, pinv); // (-p..p)
                        tmp = vec1d_reduce_2n_to_n(u+v, p); // [0..p)
                        x2[a+j] = vec1d_reduce_2n_to_n(tmp + x2[a+j], p); // [0..p)
                        t++;
                    }
#else
                    for (j = 0; j < b2; j++)
                    {
                        vec1d u, v, tmp;
                        u = x2[t];
                        v = x2[a+t];
                        x2[t] = vec1d_mulmod(vec1d_sub(u, v), powers_rho[t], p, pinv); // (-p..p)
                        tmp = vec1d_reduce_2n_to_n(u+v, p); // [0..p)
                        x2[a+j] = vec1d_reduce_2n_to_n(tmp + x2[a+j], p); // [0..p)
                        t++;
                    }
#endif
                }
        }
        // remaining part of x2 still in [0..p)
        _fft(x2, Q, p, pinv, k);  // (-3p..3p)
        x2 += a;
    }
    while (b != 0);
    
    x2 = x2_bak;

    for (i = 0; i < N; i += 4)
    {
        vec4d u, v;
        u = vec4d_load_aligned(x2 + i);
        v = vec4d_reduce_to_0n(u, p4, pinv4);
        vec4d_store_aligned(x2 + i, v);
    }
    for (; i < N; i++)
    {
        x[i] = vec1d_reduce_to_0n(x2[i], p, pinv); // [0..p)
    }
    for (i = 0; i < N; i++) // we redo a few here, does not matter
    {
        x[i] = x2[i];
    }

    free(x2);
}

/*-----------------------------------------------------------*/
/* transpose CRT (= decomposition of sequences)              */
/* in place                                                  */
/* input/output has length N                                 */
/* tmp has length 2N, all zero                               */
/* inputs and outputs in [0..p)                              */
/*-----------------------------------------------------------*/
static inline
void CRT_t(vec1d *x, vec1d *tmp, ulong n, vec1d half, vec1d p, vec1d pinv)
{
    ulong a, b, b2, n2, i, t, nn;
    vec4d half4, p4, pinv4;
    
    if (n <= 1)
    {
        return;
    }
        
    nn = n;
    a = 1;
    n2 = (nn >> 1);
    while (a <= n2)
    {
        a <<= 1;
    }

    half4 = vec4d_set_d(half);
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);

    while (nn != a)
    {
        
        nn = nn-a;
        b = 1;
        n2 = (nn >> 1);
        while (b <= n2)
        {
            b <<= 1;
        }
        b2 = b << 1;

#ifdef TRY_AVX2
        i = 0;
        if (b >= 4)
        {
            for (; i < b; i += 4)
            {
                vec4d u;
                u = vec4d_mulmod(half4, vec4d_add(vec4d_load_aligned(x + i), vec4d_load_aligned(x + a + i)),
                                 p4, pinv4);
                u = vec4d_reduce_2n_to_n(vec4d_add(u, p4), p4);
                vec4d_store_aligned(x + a + i, u);
                vec4d_store_aligned(tmp + i, u);
            }
        }
        for (; i < b; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[i] + x[a + i], p, pinv);
            u = vec1d_reduce_2n_to_n(u + p, p);
            x[a + i] = u;
            tmp[i] = u;
        }

        i = b;
        if (b >= 4)
            for (; i < nn; i += 4)
            {
                vec4d u;
                u = vec4d_mulmod(half4, vec4d_add(vec4d_load_aligned(x + i), vec4d_load_aligned(x + a + i)),
                                 p4, pinv4);
                u = vec4d_reduce_2n_to_n(vec4d_add(u, p4), p4);
                vec4d_store_aligned(x + a + i, u);
            }
        for (; i < nn; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[i] + x[a + i], p, pinv);
            x[a + i] = vec1d_reduce_2n_to_n(u + p, p);
        }
#else
        for (i = 0; i < b; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[i] + x[a + i], p, pinv);
            u = vec1d_reduce_2n_to_n(u + p, p);
            x[a + i] = u;
            tmp[i] = u;
        }
        for (i = b; i < nn; i++)
        {
            vec1d u;
            u = vec1d_mulmod(half, x[i] + x[a + i], p, pinv);
            x[a + i] = vec1d_reduce_2n_to_n(u + p, p);
        }
#endif
        tmp += b2;
        x += a;
        a = b;
    }

    while (nn != n)
    {
        b = a;
        a <<= 1;
        b2 = a;
        while (!(a & n))
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
                vec4d u, v;
                v = vec4d_load_aligned(x + a + i);
                u =vec4d_sub(vec4d_load_aligned(tmp + i), vec4d_add(v, v));
                u = vec4d_reduce_2n_to_n(u + p4, p4);
                vec4d_store_aligned(tmp + b + i, vec4d_reduce_2n_to_n(u + p4, p4));
            }
        for (; i < b; i++)
        {
            vec1d u;
            u = tmp[i] - (x[a+i] + x[a+i]);
            u = vec1d_reduce_2n_to_n(u + p, p);
            tmp[b + i] = vec1d_reduce_2n_to_n(u + p, p);
        }
#else
        for (i = 0; i < b; i++)
        {
            vec1d u;
            u = tmp[i] - (x[a+i] + x[a+i]);
            u = vec1d_reduce_2n_to_n(u + p, p);
            tmp[b + i] = vec1d_reduce_2n_to_n(u + p, p);
        }
#endif
        
        t = 0;
        while (t < a)
        {
#ifdef TRY_AVX2
            i = 0;
            if (b >= 2)
                for (; i < b2; i += 4, t += 4)
                {
                    vec4d u;
                    u = vec4d_sub(vec4d_load_aligned(x + t), vec4d_load_aligned(tmp + i));
                    vec4d_store_aligned(x + t, vec4d_reduce_2n_to_n(vec4d_add(u, p4), p4));
                }
            for (; i < b2; i++, t++)
            {
                vec1d u;
                u = x[t] - tmp[i];
                x[t] = vec1d_reduce_2n_to_n(u + p, p);
            }
#else
            for (i = 0; i < b2; i++, t++)
            {
                vec1d u;
                u = x[t] - tmp[i];
                x[t] = vec1d_reduce_2n_to_n(u + p, p);
            }
#endif
        }
        nn = nn+a;
    }
}

/*-----------------------------------------------------------*/
/* transpose inverse tft in length N                         */
/*-----------------------------------------------------------*/
void nmod_sd_tft_interpolate_t(mp_ptr x, mp_srcptr A, sd_fft_lctx_t Q, nmod_sd_fft_t F, const ulong N)
{
    ulong i, n, a, k, Nup;
    vec1d half, p, pinv;
    vec4d p4, pinv4;
    vec1d *wk, *wk2, *wk_bak;
    
    if (N == 0)
    {
        return;
    }
    
    if (N == 1)
    {
        x[0] = A[0];
        return;
    }

    Nup = FLINT_MAX(4, n_next_pow2m1(N-1) + 1);
    wk = (vec1d *) aligned_alloc(32, 3 * Nup * sizeof(vec1d));

    wk2 = wk + Nup;
    wk_bak = wk;
    
    for (i = 0; i < N; i++)
    {
        wk[i] = A[i];
    }
    for (; i < 3 * Nup; i++)
    {
        wk[i] = 0;
    }

    p = Q->p;
    pinv = Q->pinv;
    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);

    half = (vec1d) (1 + ((ulong)(p) >> 1));

    CRT_t(wk, wk2, N, half, p, pinv);
    
    n = N;
    a = 1;
    k = 0;
    while (a <= n/2)
    {
        k++;
        a <<= 1;
    }
    
    do
    {
        vec1d *powers_w = F->powers_inv_w_over_2[k];
#ifdef TRY_AVX2        
        i = 0;
        if (a >= 4)
            for (; i < a; i += 4)
            {
                vec4d_store_aligned(wk + i, vec4d_mulmod(vec4d_load_aligned(wk + i),
                                                         vec4d_load_aligned(powers_w + i), p4, pinv4));
            }
        for (; i < a; i++)
        {
            wk[i] = vec1d_mulmod(wk[i], powers_w[i], p, pinv); // (-p..p)
        }
#else
        for (i = 0; i < a; i++)
        {
            wk[i] = vec1d_mulmod(wk[i], powers_w[i], p, pinv);
        }
#endif
        _fft(wk, Q, p, pinv, k);

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


#ifdef TRY_AVX2
    for (i = 0; i < N; i += 4)
    {
        vec4d u, v;
        u = vec4d_load_aligned(wk_bak + i);
        v = vec4d_reduce_to_0n(u, p4, pinv4);
        vec4d_store_aligned(wk_bak + i, v);
    }
    for (; i < N; i++)
    {
        x[i] = vec1d_reduce_to_0n(wk_bak[i], p, pinv); // [0..p)
    }
    for (i = 0; i < N; i++) // we redo a few here, does not matter
    {
        x[i] = wk_bak[i];
    }
#else
    for (i = 0; i < N; i++)
    {
        x[i] = (mp_limb_t) vec1d_reduce_to_0n(wk_bak[i], p, pinv);
    }
#endif
    flint_free(wk_bak);
}
