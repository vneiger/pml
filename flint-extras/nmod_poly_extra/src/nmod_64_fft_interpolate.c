#include "nmod_poly_extra.h"

#ifdef HAS_INT128

/*------------------------------------------------------------*/
/*  in-place size 1 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void _inv_fft_1(mp_ptr x, const nmod_t mod)
{
    mp_limb_t u0, u1, t0, t1;
    
    u0 = x[0];
    u1 = x[1];
    BUTTERFLY(t0, t1, u0, u1, mod);
    x[0] = t0;
    x[1] = t1;
}

/*------------------------------------------------------------*/
/*  in-place size 2 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void _inv_fft_2(mp_ptr x, const mp_limb_t w, const nmod_t mod)
{
    mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
    
    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];

    BUTTERFLY(v0, v1, u0, u1, mod);
    BUTTERFLY(v2, v3, u2, u3, mod);
    BUTTERFLY(z0, z1, v0, v2, mod);
    v3 = nmod_mul(v3, w, mod);
    BUTTERFLY(z2, z3, v1, v3, mod);

    x[0] = z0;
    x[1] = z2;
    x[2] = z1;
    x[3] = z3;
}

/*------------------------------------------------------------*/
/*  in-place size 2^k inverse butterfly, k >= 3               */
/*------------------------------------------------------------*/
static void _inv_fft_k(mp_ptr x, const mp_ptr powers_inv_w_in, const mp_ptr i_powers_inv_w_in,
                       const nmod_t mod, const ulong k)
{
    mp_limb_t w, iw, p, p2;
    ulong N, M;
    ulong r, i;
    mp_ptr powers_inv_w, i_powers_inv_w, x0, x1, x2;

    N = 4;
    M = 1L << (k-2);    // number of blocks of size 4 = n/4
    powers_inv_w = powers_inv_w_in;
    i_powers_inv_w = i_powers_inv_w_in;
    w = powers_inv_w[1];
    iw = i_powers_inv_w[1];

    p = mod.n;
    p2 = p+p;
    
    // first two layers: no reduction
    x2 = x;
    for (r = 0; r < M; r++, x2 += 4)
    {
        mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3;
        u0 = x2[0];
        u1 = x2[1];
        u2 = x2[2];
        u3 = x2[3];
        
        v0 = u0 + u1; // [0,2p)
        v2 = u2 + u3; // [0,2p)
        v1 = u0 - u1 + p; // [0,2p)
        v3 = u2 - u3 + p; // [0,2p)
        v3 = mul_mod_precon_unreduced(v3, w, p, iw); // [0,2p)
        
        x2[0] = v0 + v2; // [0,4p)
        x2[1] = v1 + v3; // [0,4p)
        x2[2] = (v0 + p2) - v2; // [0,4p)
        x2[3] = (v1 + p2) - v3; // [0, 4p) 
    }
    
    powers_inv_w += N/2;
    i_powers_inv_w += N/2;
    N *= 2;
    M /= 2;
    
    // middle layers
    for (; M >= 2; N *= 2, M /= 2)
    {
        x0 = x;
        x1 = x + N/2;
        for (r = 0; r < M; r++, x0 += N, x1 += N)
        {
            for (i = 0; i < N/2; i+=2)
            {
                mp_limb_t u0, u1;
                mp_limb_signed_t z;
                
                z = x0[i] - p2; // (-2p,2p)
                u0 = (z < 0) ? (mp_limb_t)(z + p2) : (mp_limb_t) z; // [0,2p)
                u1 = mul_mod_precon_unreduced(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
                x0[i] = u0 + u1; // [0,4p)
                x1[i] = (u0 + p2) - u1; // [0,4p)

                z = x0[i+1] - p2; // (-2p,2p)
                u0 = (z < 0) ? (mp_limb_t) (z + p2) : (mp_limb_t) z; // [0,2p)
                u1 = mul_mod_precon_unreduced(x1[i+1], powers_inv_w[i+1], p, i_powers_inv_w[i+1]);  // [0,2p)
                x0[i+1] = u0 + u1; // [0,4p)
                x1[i+1] = (u0 + p2) - u1; // [0,4p)
            }
        }
        powers_inv_w += N/2;
        i_powers_inv_w += N/2;
    }

    // last layer, with full reduction
    x0 = x;
    x1 = x + N/2;

    for (i = 0; i < N/2; i+=2)
    {
        mp_limb_t u0, u1, q;
        mp_limb_signed_t z;
        
        z = x0[i] - p2;
        u0 = (z < 0) ? (mp_limb_t) (z + p2) : (mp_limb_t) z;
        u1 = mul_mod_precon_unreduced(x1[i], powers_inv_w[i], p, i_powers_inv_w[i]);  // [0,2p)
        
        q = u0 + u1; 
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x0[i] = q;
        
        q = (u0 + p2) - u1;
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x1[i] = q;
        
        z = x0[i+1] - p2;
        u0 = (z < 0) ? (mp_limb_t) (z + p2) : (mp_limb_t) z; // [0,2p)
        u1 = mul_mod_precon_unreduced(x1[i+1], powers_inv_w[i+1], p, i_powers_inv_w[i+1]);  // [0,2p)
        
        q = u0 + u1; // [0,4p)
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x0[i+1] = q;
        
        q = (u0 + p2) - u1; // [0,4p)
        q -= (q >= p2) ? p2 : 0;
        q -= (q >= p) ? p : 0;
        x1[i+1] = q;
    }
}


/*---------------------------------------------------------*/
/* reduces s mod (X^sz-1), assuming s has length N         */
/* assumes sz divides N                                    */
/* all calculations are mod p                              */
/*---------------------------------------------------------*/
static inline
void fold_minus(mp_ptr x, mp_ptr s, ulong sz, ulong N, mp_limb_t p)
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
                x[0] = n_addmod(x[0], s[i], p);
                i++;
                x[1] = n_addmod(x[1], s[i], p);
                i++;
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
                    x[j] = n_addmod(x[j], s[i], p);
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
void CRT(mp_ptr x, mp_ptr tmp, ulong n, mp_limb_t p)
{
    ulong a, b, b2, n2, i, nn;
    mp_limb_t half;

    nn = n;
    if (nn <= 1)
    {
        return;
    }

    half = 1 + (p >> 1);

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
            x[a+i] = n_addmod(n_addmod(tmp[b+i], tmp[b+i], p), x[a+i], p);
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
            mp_limb_t u;
            u = n_addmod(tmp[i], tmp[b+i], p);
            u = n_submod(x[a+i], u, p);
            x[a+i] = (u >> 1) + (u & 1)*half;   // u/2 mod p
            x[i] = n_addmod(x[i], x[a+i], p);
        }
        
        for (i = b; i < n; i++)
        {
            x[a+i] = (x[a+i] >> 1) + (x[a+i] & 1)*half;   // u/2 mod p
            x[i] = n_addmod(x[i], x[a+i], p);
        }
        
        n = n + a;
    }
}


/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_64_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_64_fft_t F, const ulong k)
{
    ulong i, N;
    
    N = 1L << k;
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = x[i];
    }

    // special cases
    switch(k)
    {
        case 0:
            break;
        case 1:
            _inv_fft_1(poly->coeffs, F->mod);
            break;
        case 2:
            _inv_fft_2(poly->coeffs, F->powers_inv_w[2][3], F->mod);
            break;
        default:
            _inv_fft_k(poly->coeffs, F->powers_inv_w[k]+2, F->i_powers_inv_w[k]+2, F->mod, k);
            break;
    }
 
    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = mul_mod_precon(poly->coeffs[i], F->powers_inv_2[k], F->mod.n, F->i_powers_inv_2[k]);
    }

    _nmod_poly_normalise(poly);
}


/*------------------------------------------------------------*/
/* inverse tft                                                */
/*------------------------------------------------------------*/
void nmod_64_tft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_64_fft_t F, const ulong N)
{
    ulong i, nn, k, aa;
    mp_ptr wk, wk2, powers, i_powers;
    mp_limb_t p;
    
    nmod_poly_fit_length(poly, N);
    poly->length = N;

    wk = _nmod_vec_init(3*N);
    for (i = 0; i < N; i++)
    {
        wk[i] = x[i];
        wk[i+N] = 0;
        wk[i+2*N] = 0;
    }
    wk2 = wk;
    p = F->mod.n;
    
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
        switch(k)
        {
            case 0:
                break;
            case 1:
                _inv_fft_1(wk, F->mod);
                break;
            case 2:
                _inv_fft_2(wk, F->powers_inv_w[2][3], F->mod);
                break;
            default:
                _inv_fft_k(wk, F->powers_inv_w[k]+2, F->i_powers_inv_w[k]+2, F->mod, k);
                break;
        }
        
        powers = F->powers_inv_w_over_2[k];
        i_powers = F->i_powers_inv_w_over_2[k];

        for (i = 0; i < aa; i++)
        {
            wk[i] = mul_mod_precon(wk[i], powers[i], p, i_powers[i]);
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
    CRT(wk, wk + N, N, p);

    for (i = 0; i < N; i++)
    {
        poly->coeffs[i] = wk[i];
    }
    
    _nmod_vec_clear(wk);
    _nmod_poly_normalise(poly);
}


/*-----------------------------------------------------------*/
/* transpose tft in length N                                 */
/*-----------------------------------------------------------*/
void nmod_64_tft_evaluate_t(mp_ptr x, mp_srcptr A, const nmod_64_fft_t F, const ulong N)
{
    ulong n, a, b, b2, t, i, j, k, ell, lambda;
    mp_limb_t p;
    mp_ptr wk, wk2, powers_rho, i_powers_rho;
    
    if (N == 0)
    {
        return;
    }
    
    if (N == 1)
    {
        x[0] = A[0];
        return;
    }

    wk = _nmod_vec_init(2*N);
    
    for (i = 0; i < N; i++)
    {
        wk[i] = A[i];
        wk[i + N] = 0;
    }

    p = F->mod.n;
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
        // special cases
        switch(k)
        {
            case 0:
                break;
            case 1:
                _inv_fft_1(wk, F->mod);
                break;
            case 2:
                _inv_fft_2(wk, F->powers_w_t[2][3], F->mod);
                break;
            default:
                _inv_fft_k(wk, F->powers_w_t[k]+2, F->i_powers_w_t[k]+2, F->mod, k);
                break;
        }
        
        powers_rho = F->powers_w[k+1];
        i_powers_rho = F->i_powers_w[k+1];
        for (i = 0; i < a; i++)
        {
            wk[i] = mul_mod_precon(wk[i], powers_rho[i], p, i_powers_rho[i]);
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
        if (wk[N-a+i] != 0) // neg_mod
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
                    mp_limb_t u, v;
                    u = wk[N-n+t];
                    v = wk[N-n+a];
                    wk[N-n+t] = n_addmod(v, u, p);
                    wk[N-n+a+t] = n_submod(v, u, p);
                    t++;
                    u = wk[N-n+t];
                    v = wk[N-n+a+1];
                    wk[N-n+t] = n_addmod(v, u, p);
                    wk[N-n+a+t] = n_submod(v, u, p);
                    t++;
                }
                
                mp_limb_t u, v;
                u = wk[N-n];
                v = wk[N-n+a];
                wk[N-n] = n_addmod(v, u, p);
                wk[N-n+a] = n_submod(v, u, p);
                u = wk[N-n+1];
                v = wk[N-n+a+1];
                wk[N-n+1] = n_addmod(v, u, p);
                wk[N-n+a+1] = n_submod(v, u, p);
                break;
            case 2:
            default:
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j++)
                    {
                        mp_limb_t u, v;
                        u = wk[N-n+t];
                        v = wk[N-n+a+j];
                        wk[N-n+t] = n_addmod(v, u, p);
                        wk[N-n+a+t] = n_submod(v, u, p);
                        t++;
                    }
                }
                for (j = 0; j < b2; j++)
                {
                    mp_limb_t u, v;
                    u = wk[N-n+j];
                    v = wk[N-n+a+j];
                    wk[N-n+j] = n_addmod(v, u, p);
                    wk[N-n+a+j] = n_submod(v, u, p);
                }
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = wk[i];
    }

    _nmod_vec_clear(wk);
}







#endif
