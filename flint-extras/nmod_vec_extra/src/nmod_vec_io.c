#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PRETTY PRINTING THE VECTOR                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void _nmod_vec_print_pretty(mp_ptr vec, slong len, nmod_t mod)
{
    flint_printf("<size-%wd vector over Z/nZ with n = %ld>\n", len, mod.n);
    flint_printf("[");
    for (slong i = 0; i < len-1; i++)
        flint_printf("%ld, ", vec[i]);
    flint_printf("%ld]\n", vec[len-1]);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
