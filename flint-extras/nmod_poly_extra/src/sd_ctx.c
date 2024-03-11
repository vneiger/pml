#include <flint/nmod.h>
#include <flint/fft_small.h>

void sd_fft_ctx_init_inverse(sd_fft_ctx_t Qt, sd_fft_ctx_t Q)
{
    ulong m, l, NN, i, k;
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
    
    Qt->w2tab_depth = SD_FFT_CTX_INIT_DEPTH;
    
    for (k = SD_FFT_CTX_INIT_DEPTH; k < FLINT_BITS; k++)
        Qt->w2tab[k] = NULL;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
