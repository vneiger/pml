#include "nmod_poly_extra.h"

void nmod_cfft_reduce(mp_ptr x, const mp_limb_t p, const slong N)
{
    slong nn, k, a;
    
    if (N <= 1)
    {
        return;
    }
  
    nn = N;
    k = 0;
    a = 1;
    while (a <= nn/2)
    {
        k++;
        a <<= 1;
    }
  
    // if N is a power of 2, nothing to do
    if (a == nn)
    {
        return;
    }
    
    // else, loop through the reduction process
    while (a != nn)
    {
        slong ell, lambda, b, b2, t;
        
        b = a >> 1;
        ell = k-1;
        while (! (nn & b))
        {
            b >>= 1;
            ell--;
        }
        
        t = 0;
        b2 = b << 1;

        // b = 1: unroll the loop
        if (b == 1)
        {
            mp_limb_t u, v;
            u = x[0];
            v = x[a+0];
            x[0] = n_submod(u, v, p);
            x[a+0] = n_addmod(u, v, p);
            u = x[1];
            v = x[a+1];
            x[1] = n_submod(u, v, p);
            x[a+1] = n_addmod(u, v, p);
            
            long t = 2;
            long lambda = 1L << (k-ell-1);
            for (long i = 1; i < lambda; i++)
            {
                long u, v;
                u = x[t];
                v = x[a+t];
                x[t] = u - v + p;
                x[a] = n_addmod(n_addmod(u, v, p), x[a], p);
                t++;
                u = x[t];
                v = x[a+t];
                x[t] = u - v + p;
                x[a+1] = n_addmod(n_addmod(u, v, p), x[a+1], p);
                t++;
            }
        }
        else
        {
            // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
            // i = 0 : special case
            for (; t < b2; t++)
            {
                mp_limb_t u, v;
                u = x[t];
                v = x[a+t];
                x[t] = n_submod(u, v, p);
                x[a+t] = n_addmod(u, v, p);
            }
      
            lambda = 1L << (k-ell-1);
            for (long i = 1; i < lambda; i++)
            {
                for (long j = 0; j < b2; j++)
                {
                    mp_limb_t u, v;
                    u = x[t];
                    v = x[a+t];
                    x[t] = n_submod(u,  v, p);
                    x[a+j] = n_addmod(n_addmod(u, v, p), x[a+j], p);
                    t++;
                }
            }
        }
        
        x += a;
        nn -= a;
        a = b;
        k = ell;
    }
    
    for (long t = 0; t < a; t++)
    {
        x[t] = n_submod(x[t], x[a+t], p);
    }
}
