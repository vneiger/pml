#include "nmod_poly_extra.h"

#ifdef HAS_INT128


/*------------------------------------------------------------*/
/*  in-place size 1 butterfly                                 */
/*  o1, o2 = i1+i2, i1-i2 mod p                               */
/*------------------------------------------------------------*/
static inline void _fft_1(mp_ptr x, const nmod_t mod)
{
    mp_limb_t u0, u1, t0, t1;

    u0 = x[0];
    u1 = x[1];
    BUTTERFLY(t0, t1, u0, u1, mod);
    x[0] = t0;
    x[1] = t1;
}

/*------------------------------------------------------------*/
/*  in-place size 2 butterfly                                 */
/*  t1, t2, t3, t4 = i1+i3, i2+i4, i1-i3, w(i2-i4) mod p      */
/*  o1, o2, o3, o4 = t1+t2, t1-t2, t3+t4, t3-t4 mod p         */
/*------------------------------------------------------------*/
static inline void _fft_2(mp_ptr x, const mp_limb_t w, const nmod_t mod)
{
    mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;

    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];
    
    BUTTERFLY(v0, v2, u0, u2, mod);
    BUTTERFLY(v1, v3, u1, u3, mod);
    v3 = nmod_mul(v3, w, mod);
    BUTTERFLY(z0, z1, v0, v1, mod);
    BUTTERFLY(z2, z3, v2, v3, mod);
    
    x[0] = z0;
    x[1] = z1;
    x[2] = z2;
    x[3] = z3;
}


/*------------------------------------------------------------*/
/* in-place FFT in size 2^k, k >= 3                           */
/* powers_w_in are the roots of 1 needed                      */
/*------------------------------------------------------------*/
static void _fft_k(mp_ptr x, const mp_ptr powers_w_in, const mp_ptr i_powers_w_in, const nmod_t mod, const ulong k)
{
    // N=size of block, M=number of blocks
    ulong N, M;
    mp_limb_t w, iw, p, p2;
    ulong r, i;
    mp_ptr powers_w, i_powers_w, x0, x1;

    powers_w = powers_w_in;
    i_powers_w = i_powers_w_in;
    
    N = 1L<<k;
    M = 1;

    p = mod.n;
    p2 = p+p;
    
    // the first layer is different: no need to reduce the sums
    x0 = x;
    x1 = x + N/2;
    for (r = 0; r < M; r++, x0 += N, x1 += N)
    {
        for (i = 0; i < N/2; i+=2)
        {
            mp_limb_t u0, u1, t0, t1;

            u0 = x0[i];
            u1 = x1[i];
            t0 = u0 + u1;
            t1 = (u0 + p) - u1;
            
            x0[i] = t0;
            x1[i] = mul_mod_precon_unreduced(t1, powers_w[i], p, i_powers_w[i]);

            u0 = x0[i+1];
            u1 = x1[i+1];

            t0 = u0 + u1;
            t1 = (u0  + p) - u1;
            t1 = mul_mod_precon_unreduced(t1, powers_w[i+1], p, i_powers_w[i+1]);
            x0[i+1] = t0;
            x1[i+1] = t1;
        }
    }
    powers_w += N/2;
    i_powers_w += N/2;
    N /= 2;
    M *= 2;
    
    for (; N > 4; N /= 2, M *= 2)
    {
        x0 = x;
        x1 = x + N/2;

        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            ulong i;
            for (i = 0; i < N/2; i+=2)
            {
                mp_limb_t u0, u1, t1;
                mp_limb_signed_t t0;
                
                u0 = x0[i];
                u1 = x1[i];
                t0 = u0 + u1 - p2;
                x0[i] = (t0 < 0) ? (mp_limb_t) (t0+p2) : (mp_limb_t) t0;
                t1 = (u0 + p2) - u1;
                x1[i] = mul_mod_precon_unreduced(t1, powers_w[i], p, i_powers_w[i]);

                u0 = x0[i+1];
                u1 = x1[i+1];
                t0 = u0 + u1 - p2;
                x0[i+1] = (t0 < 0) ? (mp_limb_t) (t0+p2) : (mp_limb_t) t0;
                t1 = (u0 + p2) - u1;
                x1[i+1] = mul_mod_precon_unreduced(t1, powers_w[i+1], p, i_powers_w[i+1]);
            }
        }
        powers_w += N/2;
        i_powers_w += N/2;
    }
    
    // last two layers
    w = powers_w[1];
    iw = i_powers_w[1];

    for (r = 0; r < M; r++, x += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
        
        u0 = x[0];
        u1 = x[1];
        u2 = x[2];
        u3 = x[3];

        v0 = u0 + u2;
        v0 -= (v0 >= p2) ? (p2) : 0;
        
        v2 = u0 - u2;
        v2 += ((mp_limb_signed_t) v2 < 0) ? (p2) : 0;
        
        v1 = u1 + u3;
        v1 -= (v1 >= p2) ? (p2) : 0;

        v3 = u1 - u3 + p2;
        v3 = mul_mod_precon_unreduced(v3, w, p, iw);
        
        z0 = v0 + v1;
        z0 -= (z0 >= p2) ? (p2) : 0;
        z0 -= (z0 >= p) ? (p) : 0;

        z1 = v0 - v1;
        z1 += ((mp_limb_signed_t) z1 < 0) ? (p2) : 0;
        z1 -= (z1 >= p) ? (p) : 0;
        
        z2 = v2 + v3;
        z2 -= (z2 >= p2) ? (p2) : 0;
        z2 -= (z2 >= p) ? (p) : 0;

        z3 = v2 - v3;
        z3 += ((mp_limb_signed_t) z3 < 0) ? (p2) : 0;
        z3 -= (z3 >= p) ? (p) : 0;
        
        x[0] = z0;
        x[1] = z1;
        x[2] = z2;
        x[3] = z3;
    }

}

/*-----------------------------------------------------------*/
/* transpose CRT (= decomposition of sequences)              */
/* in place                                                  */
/* input/output has length N                                 */
/* tmp has length 2N, all zero                               */
/*-----------------------------------------------------------*/
static void CRT_t(mp_ptr x, mp_ptr tmp, ulong n, mp_limb_t p)
{
    ulong a, b, b2, n2, i, t, nn;
    mp_limb_t half;
    
    nn = n;
    half = 1 + (p >> 1);
    
    if (n <= 1)
    {
        return;
    }

    a = 1;
    n2 = (nn >> 1);
    while (a <= n2)
    {
        a <<= 1;
    }
    
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
        
        for (i = 0; i < b; i++)
        {
            mp_limb_t u;
            u = n_addmod(x[i], x[a+i], p);
            u = (u >> 1) + (u & 1)*half;   // u/2 mod p
            x[a+i] = u;
            tmp[i] = u;
        }

        for (i = b; i < nn; i++)
        {
            mp_limb_t u;
            u = n_addmod(x[i], x[a+i], p);
            x[a+i] = (u >> 1) + (u & 1)*half;   // u/2 mod p
        }
        
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
        
        for (i = 0; i < b; i++)
        {
            tmp[b + i] = n_submod(tmp[i], n_addmod(x[a + i], x[a + i], p), p);
        }
        t = 0;
        while (t < a)
        {
            for (i = 0; i < b2; i++, t++)
            {
                x[t] = n_submod(x[t], tmp[i], p);
            }
        }
        nn = nn+a;
    }
}


/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_64_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_64_fft_t F, const ulong k)
{
    ulong i, N;
    
    N = 1L << k;
    for (i = 0; i < N; i++)
    {
        x[i] = nmod_poly_get_coeff_ui(poly, i);
    }

    // special cases
    switch(k)
    {
        case 0:
            break;
        case 1:
            _fft_1(x, F->mod);
            break;
        case 2:
            _fft_2(x, F->powers_w[2][1], F->mod);
            break;
        default:
            _fft_k(x, F->powers_w[k], F->i_powers_w[k], F->mod, k);
            break;
    }
}


/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_64_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_64_fft_t F, const ulong N)
{
    ulong i, j, k, a, ell, b, t, b2, lambda;
    mp_limb_t p;
    mp_ptr x2, x2_bak, powers_rho, i_powers_rho;
    mp_limb_t u, v;

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
    p = F->mod.n;
    x2 = _nmod_vec_init(2*N);
    x2_bak = x2;
    
    for (i = 0; i < N; i++)
    {
        x2[i] = nmod_poly_get_coeff_ui(poly, i);
        x2[i + N] = 0;
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

        powers_rho = F->powers_w[k+1];
        i_powers_rho = F->i_powers_w[k+1];
        switch(b)
        {
            case 0:
                for (i = 0; i < a; i++)
                {
                    x2[i] = mul_mod_precon(n_submod(x2[i], x2[i+a], p), powers_rho[i], p, i_powers_rho[i]);
                }
                break;
                
                // b = 1: unroll the inner loop
            case 1:
                u = x2[0];
                v = x2[a+0];
                x2[0] = mul_mod_precon(n_submod(u, v, p), powers_rho[0], p, i_powers_rho[0]);
                x2[a+0] = n_addmod(u, v, p);
                u = x2[1];
                v = x2[a+1];
                x2[1] = mul_mod_precon(n_submod(u, v, p), powers_rho[1], p, i_powers_rho[1]);
                x2[a+1] = n_addmod(u, v, p);
                
                t = 2;
                lambda = 1L << (k-1);
                for (i = 1; i < lambda; i++)
                {
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a] = n_addmod(n_addmod(u, v, p), x2[a], p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a+1] = n_addmod(n_addmod(u, v, p), x2[a+1], p);
                    t++;
                }
                break;

                // b = 2: unroll the inner loop
            case 2:
                u = x2[0];
                v = x2[a+0];
                x2[0] = mul_mod_precon(n_submod(u, v, p), powers_rho[0], p, i_powers_rho[0]);
                x2[a+0] = n_addmod(u, v, p);
                u = x2[1];
                v = x2[a+1];
                x2[1] = mul_mod_precon(n_submod(u, v, p), powers_rho[1], p, i_powers_rho[1]);
                x2[a+1] = n_addmod(u, v, p);
                u = x2[2];
                v = x2[a+2];
                x2[2] = mul_mod_precon(n_submod(u, v, p), powers_rho[2], p, i_powers_rho[2]);
                x2[a+2] = n_addmod(u, v, p);
                u = x2[3];
                v = x2[a+3];
                x2[3] = mul_mod_precon(n_submod(u, v, p), powers_rho[3], p, i_powers_rho[3]);
                x2[a+3] = n_addmod(u, v, p);
                
                t = 4;
                lambda = 1L << (k-2);
                for (i = 1; i < lambda; i++)
                {
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a] = n_addmod(n_addmod(u, v, p), x2[a], p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a+1] = n_addmod(n_addmod(u, v, p), x2[a+1], p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a+2] = n_addmod(n_addmod(u, v, p), x2[a+2], p);
                    t++;
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a+3] = n_addmod(n_addmod(u, v, p), x2[a+3], p);
                    t++;
                }
                break;

            default:
                // now we know b >= 4 so b2 >= 8
                // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
                // i = 0 : special case
                for (; t < b2; t++)
                {
                    u = x2[t];
                    v = x2[a+t];
                    x2[t] = mul_mod_precon(n_submod(u, v, p), powers_rho[t], p, i_powers_rho[t]);
                    x2[a+t] = n_addmod(u, v, p);
                }
                
                lambda = 1L << (k-ell-1);
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j++)
                    {
                        u = x2[t];
                        v = x2[a+t];
                        x2[t] = mul_mod_precon(n_submod(u,  v, p), powers_rho[t], p, i_powers_rho[t]);
                        x2[a+j] = n_addmod(n_addmod(u, v, p), x2[a+j], p);
                        t++;
                    }
                }
        }

        switch(k)
        {
            case 0:
                break;
            case 1:
                _fft_1(x2, F->mod);
                break;
            case 2:
                _fft_2(x2, F->powers_w[2][1], F->mod);
                break;
            default:
                _fft_k(x2, F->powers_w[k], F->i_powers_w[k], F->mod, k);
                break;
        }
        
        x2 += a;
    }
    while (b != 0);
    
    x2 = x2_bak;
    for (i = 0; i < N; i++)
    {
        x[i] = x2[i];
    }
    _nmod_vec_clear(x2);
}


/*-----------------------------------------------------------*/
/* transpose inverse tft in length N                         */
/*-----------------------------------------------------------*/
void nmod_64_tft_interpolate_t(mp_ptr x, mp_srcptr A, const nmod_64_fft_t F, const ulong N)
{
    ulong i, n, a, k;
    mp_limb_t p;
    mp_ptr wk, wk_bak, wk2, powers, i_powers;
    
    if (N == 0)
    {
        return;
    }
    
    if (N == 1)
    {
        x[0] = A[0];
        return;
    }

    wk = _nmod_vec_init(3*N);
    wk2 = wk + N;
    wk_bak = wk;
    
    for (i = 0; i < N; i++)
    {
        wk[i] = A[i];
    }
    for (i = 0; i < 2*N; i++)
    {
        wk2[i] = 0;
    }

    p = F->mod.n;
    CRT_t(wk, wk2, N, p);
    
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
        powers = F->powers_inv_w_over_2[k];
        i_powers = F->i_powers_inv_w_over_2[k];
        for (i = 0; i < a; i++)
        {
            wk[i] = mul_mod_precon(wk[i], powers[i], p, i_powers[i]);
        }
        switch(k)
        {
            case 0:
                break;
            case 1:
                _fft_1(wk, F->mod);
                break;
            case 2:
                _fft_2(wk, F->powers_inv_w_t[2][1], F->mod);
                break;
            default:
                _fft_k(wk, F->powers_inv_w_t[k], F->i_powers_inv_w_t[k], F->mod, k);
                break;
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
    
    for (i = 0; i < N; i++)
    {
        x[i] = wk_bak[i];
    }
}


#endif
