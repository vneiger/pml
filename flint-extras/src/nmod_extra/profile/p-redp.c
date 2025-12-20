/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/profiler.h>
#include <flint/ulong_extras.h>
#include "nmod_extra.h"
#include "nmod32_vec.h"

/* ulong pbits_arr[] = { 4, 9, 14, FLINT_BITS / 2 - 2, FLINT_BITS / 2 - 1, FLINT_BITS / 2, FLINT_BITS - 8, FLINT_BITS - 1, FLINT_BITS, 0 }; */
ulong pbits_arr[] = { FLINT_BITS / 2,  FLINT_BITS / 2 + 1, FLINT_BITS - 8, FLINT_BITS - 1, FLINT_BITS, 0 };

FLINT_STATIC_NOINLINE void _loop_n_mod2_preinv(nn_ptr xr, nn_srcptr x, ulong n, ulong ninv, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod2_preinv(x[i], n, ninv);
}

FLINT_STATIC_NOINLINE void _loop_n_mod_shoup(nn_ptr xr, nn_srcptr x, ulong n, ulong n_pre, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mulmod_shoup(UWORD(1), x[i], n_pre, n);
}

FLINT_STATIC_NOINLINE void _loop_n_mod_redp_fast(nn_ptr xr, nn_srcptr x, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod_redp_fast(x[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_n_mod_redp(nn_ptr xr, nn_srcptr x, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod_redp(x[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_n32_mod_redp_fast(nn_ptr xr, nn_srcptr x, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n32_mod_redp_fast(x[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_n32_mod_redp(nn_ptr xr, nn_srcptr x, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n32_mod_redp(x[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_n32s_mod_redp_fast(n32_ptr xr_32, n32_srcptr x_32, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr_32[i] = n32s_mod_redp_fast(x_32[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_n32s_mod_redp(n32_ptr xr_32, n32_srcptr x_32, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr_32[i] = n32s_mod_redp(x_32[i], n, redp, shift);
}


/****************
*  EXPERIMENT  *
****************/

FLINT_STATIC_NOINLINE void _loop_n_mod_redp_lazy(nn_ptr xr, nn_srcptr x, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod_redp_lazy(x[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_n_mod_redp_lazy_correct(nn_ptr xr, nn_srcptr x, ulong n, ulong redp, ulong shift, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod_redp_lazy_correct(x[i], n, redp, shift);
}

FLINT_STATIC_NOINLINE void _loop_nmod_mul(nn_ptr xr, nn_srcptr x, nn_srcptr y, nmod_t mod, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = nmod_mul(x[i], y[i], mod);
}

FLINT_STATIC_NOINLINE void _loop_n_mod_mulmod_redp_lazy(nn_ptr xr, nn_srcptr x, nn_srcptr y, ulong n, ulong redp, ulong nbits, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod_mulmod_redp_lazy(x[i], y[i], n, redp, nbits);
}

FLINT_STATIC_NOINLINE void _loop_n_mod_mulmod_redp(nn_ptr xr, nn_srcptr x, nn_srcptr y, ulong n, ulong redp, ulong nbits, slong N)
{
    slong i;
    for (i = 0; i < N; i++)
        xr[i] = n_mod_mulmod_redp(x[i], y[i], n, redp, nbits);
}

int main()
{
    flint_rand_t state;
    slong N;
    ulong n;
    ulong n_pre;
    nmod_t mod;
    nmod_redp_t modp;
    nmod32_redp_t modp_32;
    nmod32s_redp_t modp_32s;
    nn_ptr x, y, xr;
    n32_ptr x_32, xr_32;
    ulong t_tmp;
    ulong pbits;
    slong i;
    slong ipbits;
    double t1, FLINT_SET_BUT_UNUSED(tcpu);

    flint_rand_init(state);

    flint_printf("All times in nanoseconds\n");

    N = 1000;

    x = _nmod_vec_init(N);
    y = _nmod_vec_init(N);
    xr = _nmod_vec_init(N);
    x_32 = _nmod32_vec_init(N);
    xr_32 = _nmod32_vec_init(N);

#define PRINTF(f) flint_printf("%30s    ", f)
#define PRINTF2(f, t) flint_printf("%30s    %g\n", f, t / N / 1e-9)

    for (ipbits = 0; (pbits = pbits_arr[ipbits]) != 0; ipbits++)
    {
        n = n_nextprime(UWORD(1) << (pbits - 1), 1);
        nmod_init(&mod, n);
        n_pre = n_mulmod_precomp_shoup(UWORD(1), n);

        nmod_redp_init(&modp, n);
        if (pbits <= 32)
        {
            nmod32_redp_init(&modp_32, n);
            nmod32s_redp_init(&modp_32s, (uint32_t)n);
        }

        flint_printf("\nn = %wu\n\n", n);

        for (i = 0; i < N; i++)
        {
            x[i] = n_randlimb(state);
        }

        TIMEIT_START;
        for (i = 0; i < N; i++)
            nmod_redp_init(&modp, n);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_redp_init", t1);

        TIMEIT_START;
        _loop_n_mod2_preinv(xr, x, n, mod.ninv, N);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("ulong_mod2_preinv", t1);
        t_tmp = xr[0];

        if (pbits < FLINT_BITS)
        {
            TIMEIT_START;
            _loop_n_mod_shoup(xr, x, n, n_pre, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("ulong_mod_shoup", t1);
            if (t_tmp != xr[0])
                flint_printf("!! Warning1 !!\n");
        }

        TIMEIT_START;
        _loop_n_mod_redp(xr, x, n, modp.redp, modp.shift, N);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("n_mod_redp", t1);
        if (t_tmp != xr[0])
            flint_printf("!! Warning2 !!\n");

        TIMEIT_START;
        _loop_n_mod_redp_fast(xr, x, n, modp.redp, modp.shift, N);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("n_mod_redp_fast", t1);
        if (t_tmp != xr[0])
            flint_printf("!! Warning3 !!\n");

        if (pbits <= 32)
        {
            for (i = 0; i < N; i++)
            {
                x[i] = n_randlimb(state);
                x[i] >>= 32;
                x_32[i] = (uint32_t)x[i];
            }

            t_tmp = x[0] % n;

            TIMEIT_START;
            _loop_n32_mod_redp(xr, x, n, modp_32.redp, modp_32.shift, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("n32_mod_redp", t1);
            if (t_tmp != xr[0])
                flint_printf("!! Warning4 !!\n");

            TIMEIT_START;
            _loop_n32_mod_redp_fast(xr, x, n, modp_32.redp, modp_32.shift, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("n32_mod_redp_fast", t1);
            if (t_tmp != xr[0])
                flint_printf("!! Warning5 !!\n");

            TIMEIT_START;
            _loop_n32s_mod_redp(xr_32, x_32, n, modp_32s.redp, modp_32s.shift, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("n32s_mod_redp", t1);
            if (t_tmp != xr_32[0])
                flint_printf("!! Warning6 !!\n");

            TIMEIT_START;
            _loop_n32s_mod_redp_fast(xr_32, x_32, n, modp_32s.redp, modp_32s.shift, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("n32s_mod_redp_fast", t1);
            if (t_tmp != xr_32[0])
                flint_printf("!! Warning7 !!\n");
        }
        else  /* pbits >= 33 */
        {
            TIMEIT_START;
            _loop_n_mod_redp_lazy(xr, x, n, modp.redp, modp.shift, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("n_mod_redp_lazy", t1);
            if (t_tmp != xr[0])
                flint_printf("!! Warning8 !!\n");

            TIMEIT_START;
            _loop_n_mod_redp_lazy_correct(xr, x, n, modp.redp, modp.shift, N);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("n_mod_redp_lazy_correct", t1);
            if (t_tmp != xr[0])
                flint_printf("!! Warning9 !!\n");
        }

        for (i = 0; i < N; i++)
        {
            x[i] = n_randlimb(state) % n;
            y[i] = n_randlimb(state) % n;
        }

        t_tmp = n_mulmod2(x[0], y[0], n);

        TIMEIT_START;
        _loop_nmod_mul(xr, x, y, mod, N);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_mul", t1);
        if ((t_tmp != xr[0]) && (t_tmp != xr[0] - n) && (t_tmp != xr[0] - 2*n))
            flint_printf("!! Warning10 !!\n");

        TIMEIT_START;
        _loop_n_mod_mulmod_redp_lazy(xr, x, y, n, modp.redp, pbits, N);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("n_mod_mulmod_redp_lazy", t1);
        if ((t_tmp != xr[0]) && (t_tmp != xr[0] - n) && (t_tmp != xr[0] - 2*n))
            flint_printf("!! Warning10 !!\n");

        TIMEIT_START;
        _loop_n_mod_mulmod_redp(xr, x, y, n, modp.redp, pbits, N);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("n_mod_mulmod_redp", t1);
        if (t_tmp != xr[0])
            flint_printf("!! Warning11 !!\n");

        flint_printf("\n");
    }

    _nmod_vec_clear(x);
    _nmod_vec_clear(y);
    _nmod_vec_clear(xr);
    _nmod32_vec_clear(x_32);
    _nmod32_vec_clear(xr_32);

    flint_rand_clear(state);
}
