/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

void
_nmod_mat_krylov_degrees(slong * d, slong m, const slong * degrees,
        slong n, int clamp)
{
    slong i;

    for (i = 0; i < m; i++)
    {
        if (degrees == NULL)
            d[i] = n;
        else if (degrees[i] < 0)
            flint_throw(FLINT_ERROR, "nmod_mat_krylov: negative degree bound\n");
        else
            d[i] = (clamp && degrees[i] > n) ? n : degrees[i];
    }
}

typedef struct
{
    slong s;
    slong i;
}
_krylov_shift_index_struct;

static int
_krylov_shift_index_cmp(const void * a, const void * b)
{
    const _krylov_shift_index_struct * x = a;
    const _krylov_shift_index_struct * y = b;

    if (x->s != y->s)
        return (x->s < y->s) ? -1 : 1;
    return (x->i < y->i) ? -1 : (x->i > y->i);
}

void
_nmod_mat_krylov_normalize_shift(slong * t, const slong * shift,
        slong m, slong D)
{
    slong k, cur;
    _krylov_shift_index_struct * v;

    if (shift == NULL)
    {
        for (k = 0; k < m; k++)
            t[k] = 0;
        return;
    }

    if (m == 0)
        return;

    /* shift already nondecreasing (for instance uniform): no sorting */
    for (k = 1; k < m && shift[k] >= shift[k - 1]; k++)
        ;

    if (k == m)
    {
        cur = 0;
        t[0] = 0;
        for (k = 1; k < m; k++)
        {
            ulong gap = (ulong) shift[k] - (ulong) shift[k - 1];

            if (gap > (ulong) D + 1)
                gap = D + 1;
            cur += gap;
            t[k] = cur;
        }
        return;
    }

    v = flint_malloc(m * sizeof(_krylov_shift_index_struct));
    for (k = 0; k < m; k++)
    {
        v[k].s = shift[k];
        v[k].i = k;
    }

    qsort(v, m, sizeof(_krylov_shift_index_struct), _krylov_shift_index_cmp);

    cur = 0;
    t[v[0].i] = 0;
    for (k = 1; k < m; k++)
    {
        /* exact difference of two slongs, as an unsigned word */
        ulong gap = (ulong) v[k].s - (ulong) v[k - 1].s;

        if (gap > (ulong) D + 1)
            gap = D + 1;
        cur += gap;
        t[v[k].i] = cur;
    }

    flint_free(v);
}

slong
_nmod_mat_krylov_nrows(slong m, const slong * degrees)
{
    slong i, N = 0;

    for (i = 0; i < m; i++)
        N += degrees[i] + 1;

    return N;
}

slong
_nmod_mat_krylov_order(slong * rows, slong * exps, slong m,
        const slong * shift, const slong * degrees)
{
    slong i, k, len, D;
    slong * key, * nextj, * heap;

    key = flint_malloc((3 * m + 1) * sizeof(slong));
    nextj = key + m;
    heap = nextj + m;

    D = 0;
    for (i = 0; i < m; i++)
        D = FLINT_MAX(D, degrees[i]);

    _nmod_mat_krylov_normalize_shift(key, shift, m, D);

    len = 0;
    for (i = 0; i < m; i++)
    {
        nextj[i] = 0;
        _nmod_mat_krylov_heap_push(heap, &len, i, key);
    }

    k = 0;
    while (len > 0)
    {
        i = _nmod_mat_krylov_heap_pop(heap, &len, key);
        rows[k] = i;
        exps[k] = nextj[i];
        k++;
        nextj[i]++;
        if (nextj[i] <= degrees[i])
        {
            key[i]++;
            _nmod_mat_krylov_heap_push(heap, &len, i, key);
        }
    }

    flint_free(key);

    return k;
}
