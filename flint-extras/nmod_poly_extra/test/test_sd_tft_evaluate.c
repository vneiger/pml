#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* returns the number of bits = 1 in n                        */
/*------------------------------------------------------------*/
static long find_length(const long nn)
{
    long n, a, i;
    n = nn;
    a = 1;
    i = 1;
    while (a <= n/2)
    {
        a <<= 1;
    }
    
    while (a != n)
    {
        long b = a >> 1;
        while (! (n & b))
        {
            b >>= 1;
        }
        i++;
        n -= a;
        a = b;
    }
    
    return i;
}

/*------------------------------------------------------------*/
/* returns the powers of 2 ni such that                       */
/*    n = sum_i n_i                                           */
/* decreasing order, the last entry is 0                      */
/*------------------------------------------------------------*/
void split_degrees(long * vec, const long nn)
{
    long i, j, n, a;
    
    i = find_length(nn);
    for (j = 0; j <= i; j++)
    {
        vec[j] = 0;
    }
    
    n = nn;
    a = 1;
    i = 0;
    while (a <= n/2)
    {
        a <<= 1;
    }
    vec[i++] = a;
    
    while (a != n)
    {
        long b;
        b = a >> 1;
        while (! (n & b))
        {
            b >>= 1;
        }
        n -= a;
        a = b;
        vec[i++] = a;
    }
}

/*------------------------------------------------------------*/
/* returns the exponents ki such that                         */
/*    n = sum_i 2^k_i                                         */
/*------------------------------------------------------------*/
void split_exponents(long * vec, const long nn)
{
    long i, j, n, a;

    i = find_length(nn);
    for (j = 0; j <= i; j++)
    {
        vec[j] = 0;
    }
    vec[i] = -1;
    
    n = nn;
    a = 1;
    i = 0;
    j = 0;
    while (a <= n/2)
    {
        a <<= 1;
        j++;
    }
    vec[i++] = j;
    
    while (a != n)
    {
        long k, b;
        k = j-1;
        b = a >> 1;
        while (! (n & b))
        {
            k--;
            b >>= 1;
        }
        n -= a;
        a = b;
        vec[i++] = k;
        j = k;
    }
}

/*------------------------------------------------------------*/
/* returns reverse_nb(u)                                      */
/*------------------------------------------------------------*/
ulong bit_reverse(ulong u, ulong nb)
{
    ulong res, i;
    res = 0;
    for (i = 0; i < nb; i++)
    {
        res <<= 1;
        res |= (u & 1);
        u >>= 1;
    }
    return res;
}

/*------------------------------------------------------------*/
/* TFT evaluation, compared to naive evaluation               */
/*------------------------------------------------------------*/
void check()
{
    ulong nmin, nmax;
    flint_rand_t state;
    ulong w0, w, p;
    nmod_t mod;
    sd_fft_ctx_t Q;
    nmod_sd_fft_t F;
    nn_ptr val, val2;
    nmod_poly_t P;
    
    flint_rand_init(state);

    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    nmod_init(&mod, p);

    
    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);
    w = nmod_pow_ui(w0, 1L<<(16-16), mod);
    nmod_sd_fft_init_set(F, w, 16, mod);

    nmin = 1;
    nmax = 10000;

    for (long n = nmin; n < nmax+1; n++)
    {
        sd_fft_lctx_t QL;
        long * degrees, * exponents;
        ulong i, j, idx, len, shift, OK, found;
        
        len = find_length(n);

        degrees = flint_malloc((len + 1) * sizeof(long));
        exponents = flint_malloc((len + 1) * sizeof(long));
        split_degrees(degrees, n);
        split_exponents(exponents, n);
        
        sd_fft_lctx_init(QL, Q, 16);

        nmod_poly_init2(P, p, n);
        for (i = 0; i < n; i++)
        {
            ulong c;
            c = n_randtest(state);
            c = nmod_mul(c, c, mod);
            c = nmod_mul(c, c, mod);
            c = nmod_mul(c, c, mod);
            nmod_poly_set_coeff_ui(P, i, c);
        }
        
        val = _nmod_vec_init(n);
        nmod_sd_tft_evaluate(val, P, QL, F, n);

        val2 = _nmod_vec_init(n);
        i = 0;
        idx = 0;
        do
        {
            slong e;
            ulong rho, g;
            
            e = exponents[i];
            rho = nmod_pow_ui(w0, 1L << (16-e-1), Q->mod);
            g = nmod_pow_ui(w0, 1L << (16-e), Q->mod);
            
            for (j = 0; j < (1L << e); j++)
            {
                val2[idx] = nmod_poly_evaluate_nmod(P, nmod_mul(rho, nmod_pow_ui(g, bit_reverse(j, e), Q->mod), Q->mod));
                idx++;
            }
            shift += (1L << e);
            i++;
        } while (exponents[i] != -1);

        
        OK = 1;
        for (i = 0; i < n; i++)
        {
            found = 0;
            for (j = 0; j < n; j++)
                if (val2[j] == val[i])
                    found = 1;
            if (found == 0)
                OK = 0;
            found = 0;
            for (j = 0; j < n; j++)
                if (val2[i] == val[j])
                    found = 1;
            if (found == 0)
                OK = 0;
        }

        nmod_poly_clear(P);
        _nmod_vec_clear(val);
        _nmod_vec_clear(val2);
        flint_free(exponents);
        flint_free(degrees);

        if (OK == 0)
        {
            printf("%ld\n", n);
            continue;
        }

        if (n > 2000)
            n += 100;
    }
    sd_fft_ctx_clear(Q);
    nmod_sd_fft_clear(F);
    flint_rand_clear(state);
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
