#include <flint/ulong_extras.h>
#include "nmod_mat_poly.h"

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* INIT                                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


void nmod_mat_poly_init_preinv(nmod_mat_poly_t matp,
                               slong r,
                               slong c,
                               ulong n,
                               ulong ninv)
{
    matp->coeffs = NULL;

    matp->alloc = 0;
    matp->length = 0;

    matp->r = r;
    matp->c = c;

    matp->mod.n = n;
    matp->mod.ninv = ninv;
    matp->mod.norm = flint_clz(n);
}

void nmod_mat_poly_init(nmod_mat_poly_t matp,
                        slong r,
                        slong c,
                        ulong n)
{
    nmod_mat_poly_init_preinv(matp, r, c, n, n_preinvert_limb(n));
}


void nmod_mat_poly_init2_preinv(nmod_mat_poly_t matp,
                                slong r,
                                slong c,
                                ulong n,
                                ulong ninv,
                                slong alloc)
{
    if (alloc)
        matp->coeffs = (nmod_mat_struct *) flint_malloc(alloc * sizeof(nmod_mat_struct));
    else
        matp->coeffs = NULL;

    matp->alloc = alloc;
    matp->length = 0;

    matp->r = r;
    matp->c = c;

    matp->mod.n = n;
    matp->mod.ninv = ninv;

    matp->mod.norm = flint_clz(n);
}

void nmod_mat_poly_init2(nmod_mat_poly_t matp,
                     slong r,
                     slong c,
                     ulong n,
                     slong alloc)
{
    nmod_mat_poly_init2_preinv(matp, r, c, n, n_preinvert_limb(n), alloc);
}

/*------------------------------------------------------------*/
/* SET / INIT SET                                             */
/*------------------------------------------------------------*/

void nmod_mat_poly_set(nmod_mat_poly_t matp1, const nmod_mat_poly_t matp2)
{
    if (matp1 != matp2)         /* Aliasing is trivial */
    {
        const slong len = matp2->length;

        nmod_mat_poly_fit_length(matp1, len);
        _nmod_mat_poly_set_length(matp1, len);

        for (slong i = 0; i < len; i++)
            nmod_mat_set(matp1->coeffs + i, matp2->coeffs + i);
    }
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CLEAR / REALLOC / FIT LENGTH                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_mat_poly_clear(nmod_mat_poly_t matp)
{
    // clear any allocated matrix coefficient
    for (slong i=0; i < matp->length; i++)
        nmod_mat_clear(matp->coeffs + i);
    // free coeffs
    if (matp->coeffs)
        flint_free(matp->coeffs);
}

void nmod_mat_poly_realloc(nmod_mat_poly_t matp, slong alloc)
{
    // clear and set to 0
    if (alloc == 0)
    {
        nmod_mat_poly_clear(matp);
        matp->length = 0;
        matp->alloc = 0;
        matp->coeffs = NULL;
        return;
    }

    if (matp->alloc) // realloc
    {
        // truncate at order `alloc`
        nmod_mat_poly_truncate(matp, alloc);
        matp->coeffs = (nmod_mat_struct *) flint_realloc(matp->coeffs, alloc * sizeof(nmod_mat_struct));
    }
    else // not allocated yet, do it now
        matp->coeffs = (nmod_mat_struct *) flint_malloc(alloc * sizeof(nmod_mat_struct));

    matp->alloc = alloc;
}

void nmod_mat_poly_fit_length(nmod_mat_poly_t matp, slong alloc)
{
    // realloc to the maximum of `alloc` and 
    // the double of the current `matp->alloc`
    if (alloc > matp->alloc)
    {
        if (alloc < 2 * matp->alloc)
            alloc = 2 * matp->alloc;

        nmod_mat_poly_realloc(matp, alloc);
    }
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
