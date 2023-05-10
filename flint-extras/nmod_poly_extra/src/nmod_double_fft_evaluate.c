#include <flint/fft_small.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

#ifdef HAS_AVX2
#include <immintrin.h>

/*------------------------------------------------------------*/
/*  in-place size 1 butterfly                                 */
/*  o1, o2 = i1+i2, i1-i2 mod p                               */
/*  input in [0..p), output in [0..p)                        */
/*------------------------------------------------------------*/
static inline void _fft_1(vec1d *x, vec1d p)
{
    vec1d u0, u1, v;

    u0 = x[0];
    u1 = x[1];
    
    v = u0 + u1 - p;
    x[0] = v < 0 ? v + p : v;
    v = u0 - u1;
    x[1] = v < 0 ? v + p : v;
}

/*------------------------------------------------------------*/
/*  in-place size 2 butterfly                                 */
/*  t1, t2, t3, t4 = i1+i3, i2+i4, i1-i3, w(i2-i4) mod p      */
/*  o1, o2, o3, o4 = t1+t2, t1-t2, t3+t4, t3-t4 mod p         */
/*  input in [0..p), output in [0..p)                        */
/*------------------------------------------------------------*/
static inline void _fft_2(vec1d *x, const vec1d w, const vec1d p, const vec1d pinv)
{
    vec1d u0, u1, u2, u3, v0, v1, v2, v3;

    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];

    v0 = u0 + u2;  // [0..2p)
    v2 = u0 - u2;  // (-p..p)
    v1 = u1 + u3;  // [0..2p)
    v3 = vec1d_mulmod(u1 - u3, w, p, pinv);      // (-p..p)

    x[0] = vec1d_reduce_to_0n(v0 + v1, p, pinv);  // [0..4p) -> [0..p)
    x[1] = vec1d_reduce_to_0n(v0 - v1, p, pinv);  // (-2p..2p) -> [0..p)
    x[2] = vec1d_reduce_to_0n(v2 + v3, p, pinv);  // (-2p..2p) -> [0..p)
    x[3] = vec1d_reduce_to_0n(v2 - v3, p, pinv);  // (-2p..2p) -> [0..p)
}


/*------------------------------------------------------------*/
/* in-place FFT in size 2^k, k >= 2                           */
/* powers_w_in are the roots of 1 needed                      */
/* output is in bit-reverse order                             */
/*------------------------------------------------------------*/
 void _fft_k(vec1d *x, const vec1d *powers_w_in, const vec1d p, const vec1d pinv, const ulong k)
{
    ulong N, M;
    vec1d w;
    vec4d w4, p4, pinv4;
    ulong r;
    const vec1d *powers_w;
    vec1d *x0, *x1, *x2, *x3;


    powers_w = powers_w_in;
    N = 1L<<k;
    M = 1;

    p4 = vec4d_set_d(p);
    pinv4 = vec4d_set_d(pinv);

    
    for (; N > 8; )
    {
        ulong N2, N4, M2;
        
        N2 = N >> 1L;
        N4 = N2 >> 1L;
        M2 = M << 1L;
        
        x0 = x;
        x1 = x + N2;
        x2 = x0 + N4;
        x3 = x1 + N4;
        
        
        for (r = 0; r < M; r++, x0 += N, x1 += N, x2 += N, x3 += N)
        {
            ulong i;
            for (i = 0; i < N4; i += 4) // we do two layers
            {
                vec4d u0, u1, u, v, x, y, z;
                
                u0 = vec4d_load(x0 + i);
                u1 = vec4d_load(x1 + i);
                w4 = vec4d_load(powers_w + i);
                u = vec4d_add(u0, u1);   // (-2p..2p)
                z = vec4d_sub(u0, u1);
                v = vec4d_mulmod(z, w4, p4, pinv4);
                
                u0 = vec4d_load(x2 + i);
                u1 = vec4d_load(x3 + i);
                w4 = vec4d_load(powers_w + i + N4);
                x = vec4d_add(u0, u1);  // (-2p..2p)
                z = vec4d_sub(u0, u1);
                y = vec4d_mulmod(z, w4, p4, pinv4);
                
                w4 = vec4d_load(powers_w + i + N2);
                
                z = vec4d_add(u, x);
                u0 = vec4d_reduce_to_pm1n(z, p4, pinv4);
                vec4d_store(x0 + i, u0);
                z = vec4d_sub(u, x); // (-4p..4p)
                u1 = vec4d_mulmod(z, w4, p4, pinv4);  // w4 is (-p/2..p/2) so product (-2p^2..2p^2), reduced to (-p..p)
                vec4d_store(x2 + i, u1);
                
                z = vec4d_add(v, y);
                u0 = vec4d_reduce_to_pm1n(z, p4, pinv4);
                vec4d_store(x1 + i, u0);
                z = vec4d_sub(v, y);
                u1 = vec4d_mulmod(z, w4, p4, pinv4);
                vec4d_store(x3 + i, u1);
            }
        }
        powers_w += N2 + N4;
        N = N4;
        M = M2 << 1L;
    }
    
    /* N = 8: 3 layers at once */
    if (N == 8)
    {
        x0 = x;
        w = powers_w[5];
        
        for (r = 0; r < M; r++, x0 += 8)
        {
            vec4d u0, u1, u, v;
            vec1d z0, z1, z2, z3, v0, v1, v2, v3;
            
            u0 = vec4d_load(x0);
            u1 = vec4d_load(x0 + 4);
            w4 = vec4d_load(powers_w);
            
            u = vec4d_add(u0, u1); // (-2p..2p)
            z0 = u[0];
            z1 = u[1];
            z2 = u[2];
            z3 = u[3];
            v0 = z0 + z2;  // (-4p..4p)
            v2 = z0 - z2;  // (-4p..4p)
            v1 = z1 + z3;  // (-4p..4p)
            v3 = vec1d_mulmod(z1 - z3, w, p, pinv);      // (-p..p)
            x0[0] = vec1d_reduce_to_0n(v0 + v1, p, pinv);  // [0..p)
            x0[1] = vec1d_reduce_to_0n(v0 - v1, p, pinv);      // [0..p)
            x0[2] = vec1d_reduce_to_0n(v2 + v3, p, pinv);  // [0..p)
            x0[3] = vec1d_reduce_to_0n(v2 - v3, p, pinv);      // [0..p)
            
            u = vec4d_sub(u0, u1);
            v = vec4d_mulmod(u, w4, p4, pinv4);
            z0 = v[0];
            z1 = v[1];
            z2 = v[2];
            z3 = v[3];
            v0 = z0 + z2;  // (-2p..2p)
            v2 = z0 - z2;  // (-2p..2p)
            v1 = z1 + z3;  // (-2p..2p)
            v3 = vec1d_mulmod(z1 - z3, w, p, pinv);      // (-p..p)
            x0[4] = vec1d_reduce_to_0n(v0 + v1, p, pinv);  // [0..p)
            x0[5] = vec1d_reduce_to_0n(v0 - v1, p, pinv);      // [0..p)
            x0[6] = vec1d_reduce_to_0n(v2 + v3, p, pinv);  // [0..p)
            x0[7] = vec1d_reduce_to_0n(v2 - v3, p, pinv);      // [0..p)
        }
        
    }
    /* N = 4: 2 layers at once */
    else
    {
        w = powers_w[1];
        
        for (r = 0; r < M; r++, x += 4)
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
            
            x[0] = vec1d_reduce_to_0n(v0 + v1, p, pinv);  // [0..p)
            x[1] = vec1d_reduce_to_0n(v0 - v1, p, pinv);      // [0..p)
            x[2] = vec1d_reduce_to_0n(v2 + v3, p, pinv);  // [0..p)
            x[3] = vec1d_reduce_to_0n(v2 - v3, p, pinv);      // [0..p)
        }
    }
}


/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_double_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_double_fft_t F, const ulong N)
{
    ulong i, j, k, a, ell, b, t, b2, lambda, N2;
    mp_limb_t u, v;
    vec1d *x2, *x2_bak, *powers_rho;
    vec1d p, pinv;
    
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
    p = F->p;
    pinv = F->pinv;
    
    N2 = n_next_pow2m1(2*N-1) + 1;

    x2 = (vec1d *) aligned_alloc(32, FLINT_MAX(N2, 4) * sizeof(vec1d));
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
        switch(b)
        {
            case 0:
                for (i = 0; i < a; i++)
                {
                    x2[i] = vec1d_mulmod(vec1d_sub(x2[i], x2[i+a]), powers_rho[i], p, pinv); // (-p..p)
                }
                break;

            default:
                // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
                // i = 0 : special case
                for (; t < b2; t++)
                {
                    u = x2[t];
                    v = x2[a+t];

                    x2[t] = vec1d_mulmod(vec1d_sub(u, v), powers_rho[t], p, pinv); // (-p..p)
                    x2[a+t] = vec1d_reduce_2n_to_n(u+v, p);
                }
                
                lambda = (1L << k) >> (ell+1);
                for (i = 1; i < lambda; i++)
                {
                    for (j = 0; j < b2; j++)
                    {
                        double tmp;
                        u = x2[t];
                        v = x2[a+t];
                        x2[t] = vec1d_mulmod(vec1d_sub(u, v), powers_rho[t], p, pinv);
                        tmp = vec1d_reduce_2n_to_n(u+v, p);
                        x2[a+j] = vec1d_reduce_2n_to_n(tmp + x2[a+j], p);
                        t++;
                    }
                }
        }

        switch(k)
        {
            case 0:
                if (x2[0] < 0)
                {
                    x2[0] += p;
                }
                break;
            case 1:
                if (x2[0] < 0)
                {
                    x2[0] += p;
                }
                if (x2[1] < 0)
                {
                    x2[1] += p;
                }
                _fft_1(x2, p);
                break;
            case 2:
                if (x2[0] < 0)
                {
                    x2[0] += p;
                }
                if (x2[1] < 0)
                {
                    x2[1] += p;
                }
                if (x2[2] < 0)
                {
                    x2[2] += p;
                }
                if (x2[3] < 0)
                {
                    x2[3] += p;
                }
                _fft_2(x2, F->powers_w[2][1], p, pinv);
                break;
            default:
                _fft_k(x2, F->powers_w[k], p, pinv, k);
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

    free(x2);
}




/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_double_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_double_fft_t F, const ulong k)
{
    ulong i, N;
    vec1d p;
    vec1d *dx;
    N = 1L << k;
    p = F->p;

    dx = (vec1d *) aligned_alloc(32, FLINT_MAX(N, 4) * sizeof(vec1d));

    
    for (i = 0; i < N; i++)
    {
        dx[i] = (vec1d) nmod_poly_get_coeff_ui(poly, i);
    }
    /* special cases */
    switch(k)
    {
        case 0:
            break;
        case 1:
            _fft_1(dx, p);
            break;
        case 2:
            _fft_2(dx, F->powers_w[2][1], F->p, F->pinv);
            break;
        default:
            _fft_k(dx, F->powers_w[k], F->p, F->pinv, k);
            break;
    }

    for (i = 0; i < N; i++)
    {
        x[i] = dx[i];
    }
    
    free(dx);

}

#endif
