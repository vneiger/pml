#include <flint/nmod.h>
#include <flint/fft_small.h>

void sd_fft_lctx_inverse(sd_fft_lctx_t Qt, sd_fft_lctx_t Q)
{
    ulong m, l, NN, i;
    vec1d p, pinv;
    vec1d *t;
    nmod_t mmod;

    
    p = Q->p;
    pinv = Q->pinv;
    nmod_init(&mmod, p);

    Qt->p = p;
    Qt->pinv = pinv;

    NN = n_pow2(SD_FFT_CTX_INIT_DEPTH - 1);
    t = (double*) flint_aligned_alloc(4096, n_round_up(NN*sizeof(double), 4096));
    Qt->w2tab[0] = t;
    t[0] = 1;

    for (m = 1, l = 1; m < SD_FFT_CTX_INIT_DEPTH; m++, l *= 2)
    {
        double* fwd;
        const double *bck;
        fwd = t + l;
        Qt->w2tab[m] = fwd;
        bck = Q->w2tab[0] + l;
        for (i = 0; i < l; i++)
        {
            fwd[i] = -bck[l - 1 - i];
        }
    }

    /* mp_limb_t p0 = (mp_limb_t) p; */
    /* for (m = 1, l = 1; m < SD_FFT_CTX_INIT_DEPTH; m++, l *= 2) */
    /* { */
    /*     ulong w0 = n_primitive_root_prime(p0); */
    /*     ulong ww = nmod_pow_ui(w0, (p0 - 1)>>(m + 1), mmod); */
    /*     ww = nmod_inv(ww, mmod); */
    /*     double w = vec1d_set_d(vec1d_reduce_0n_to_pmhn(ww, p)); */
    /*     double* curr = t + l; */
    /*     Qt->w2tab[m] = curr; */
    /*     i = 0; */
    /*     do */
    /*     { */
    /*         vec1d x = vec1d_load(t + i); */
    /*         x = vec1d_mulmod(x, w, p, pinv); */
    /*         x = vec1d_reduce_pm1n_to_pmhn(x, p); */
    /*         vec1d_store(curr + i, x); */
    /*     } */
    /*     while (i += 1, i < l); */
    /* } */
}
