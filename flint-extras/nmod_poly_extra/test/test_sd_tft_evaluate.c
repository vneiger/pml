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
    /* ulong order_w0, order_max; */
    ulong nmin, nmax;
    flint_rand_t state;
    mp_limb_t w0, p;
    nmod_t mod;
    sd_fft_ctx_t Q;
    mp_ptr val;
    nmod_poly_t P;
    
    flint_randinit(state);

    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    nmod_init(&mod, p);

    
    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);

    nmin = 1;
    nmax = 2000;

    for (long n = nmin; n < nmax+1; n++)
    {
        sd_fft_lctx_t QL;
        long * degrees, * exponents;
        ulong i;
        ulong len, shift;
        
        len = find_length(n);

        degrees = flint_malloc((len + 1) * sizeof(long));
        exponents = flint_malloc((len + 1) * sizeof(long));
        split_degrees(degrees, n);
        split_exponents(exponents, n);
        
        sd_fft_lctx_init(QL, Q, 16);

        nmod_poly_init2(P, p, n);
        for (i = 0; i < n; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        val = _nmod_vec_init(n);
        nmod_sd_tft_evaluate(val, P, QL, n);

        printf("S:=[");
        for (i = 0; i < n; i++)
        {
            printf("%ld", val[i]);
            if (i < n-1)
            {
                printf(", ");
            }
        }
        printf("];");
        printf("\n");

        printf("T:=[");
        i = 0;
        do
        {
            slong e, j;
            mp_limb_t rho, g;
            
            e = exponents[i];
            rho = nmod_pow_ui(w0, 1L << (16-e-1), Q->mod);
            g = nmod_pow_ui(w0, 1L << (16-e), Q->mod);
            
            for (j = 0; j < (1L << e); j++)
            {
                printf("%lu, ", nmod_poly_evaluate_nmod(P, nmod_mul(rho, nmod_pow_ui(g, bit_reverse(j, e), Q->mod), Q->mod)));
            }
            shift += (1L << e);
            i++;
        } while (exponents[i] != -1);
        printf("0];\n");
        printf("T:=T[1..#T-1];\n");
        printf("SequenceToSet(S) eq SequenceToSet(T);\n");
        
        nmod_poly_clear(P);
        _nmod_vec_clear(val);
        flint_free(exponents);
        flint_free(degrees);
    }
    sd_fft_ctx_clear(Q);
    flint_randclear(state);
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
