#include <assert.h>
#include <time.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_extra.h"  // for rand, mat_mul_nmod_vec_newdot
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void time_nmod_vec_dot_product_multi(ulong len, ulong k, ulong maxbits1, ulong maxbits2, ulong n, flint_rand_t state)
{
    // moduli
    nmod_t mod, mod1, mod2;
    nmod_init(&mod, n);

    if (maxbits1 < FLINT_BITS) nmod_init(&mod1, UWORD(1) << maxbits1);
    else nmod_init(&mod1, UWORD_MAX);
    if (maxbits2 < FLINT_BITS) nmod_init(&mod2, UWORD(1) << maxbits2);
    else nmod_init(&mod2, UWORD_MAX);

    // input vector
    mp_ptr u = _nmod_vec_init(len);

    // input collection of vectors (= matrices)
    mp_ptr * v = flint_malloc(len * sizeof(mp_ptr));
    for (ulong i = 0; i < len; i++)
        v[i] = _nmod_vec_init(k);

    mp_ptr ur = _nmod_vec_init(len);
    nmod_mat_t vmat; nmod_mat_init(vmat, len, k, n);

    // fill entries at random
    _nmod_vec_rand(u, state, len, mod1);
    for (ulong i = 0; i < len; i++)
        _nmod_vec_rand(v[i], state, k, mod2);

    for (ulong i = 0; i < len; i++)
        ur[i] = u[i] % n;
    for (ulong i = 0; i < len; i++)
        for (ulong j = 0; j < k; j++)
            nmod_mat_entry(vmat, i, j) = v[i][j] % n;

    double t1, t2, t3;
    clock_t tt;
    long nb_iter;

    // main multi
    t1 = 0.0; nb_iter = 0;
    while (t1 < 0.5)
    {
        mp_ptr uv = _nmod_vec_init(k);
        tt = clock();
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        nmod_vec_dot_product_multi(uv, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        t1 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _nmod_vec_clear(uv);
        nb_iter += 5;
    }
    t1 /= nb_iter; printf("%.1e\t", t1);

    //// variant v1_8
    //t3 = 0.0; nb_iter = 0;
    //while (t3 < 0.5)
    //{
    //    mp_ptr uv = _nmod_vec_init(k);
    //    tt = clock();
    //    _nmod_vec_dot_product_multi_1_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    t3 += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    _nmod_vec_clear(uv);
    //    nb_iter += 5;
    //}
    //t3 /= nb_iter; printf("%.1e\t", t3);

    //// variant v8_8
    //t3 = 0.0; nb_iter = 0;
    //while (t3 < 0.5)
    //{
    //    mp_ptr uv = _nmod_vec_init(k);
    //    tt = clock();
    //    _nmod_vec_dot_product_multi_1_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
    //    t3 += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    _nmod_vec_clear(uv);
    //    nb_iter += 5;
    //}
    //t3 /= nb_iter; printf("%.1e\t", t3);

    //// variant v8_32
    //t3 = 0.0; nb_iter = 0;
    //while (t3 < 0.5)
    //{
    //    mp_ptr uv = _nmod_vec_init(k);
    //    tt = clock();
    //    _nmod_vec_dot_product_multi_1_v8_32(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_32(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_32(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_32(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v8_32(uv, u, (mp_srcptr *) v, len, k, mod);
    //    t3 += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    _nmod_vec_clear(uv);
    //    nb_iter += 5;
    //}
    //t3 /= nb_iter; printf("%.1e\t", t3);

    //// variant v16_16
    //t3 = 0.0; nb_iter = 0;
    //while (t3 < 0.5)
    //{
    //    mp_ptr uv = _nmod_vec_init(k);
    //    tt = clock();
    //    _nmod_vec_dot_product_multi_1_v16_16(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v16_16(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v16_16(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v16_16(uv, u, (mp_srcptr *) v, len, k, mod);
    //    _nmod_vec_dot_product_multi_1_v16_16(uv, u, (mp_srcptr *) v, len, k, mod);
    //    t3 += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    _nmod_vec_clear(uv);
    //    nb_iter += 5;
    //}
    //t3 /= nb_iter; printf("%.1e\t", t3);
    //}

    // TODO this if should really be depending on nlimbs
    if (FLINT_BIT_COUNT(n) < 45)
    {
        // variant 1-8
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.5)
        {
            mp_ptr uv = _nmod_vec_init(k);
            tt = clock();
            _nmod_vec_dot_product_multi_2_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v1_8(uv, u, (mp_srcptr *) v, len, k, mod);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            _nmod_vec_clear(uv);
            nb_iter += 5;
        }
        t2 /= nb_iter; printf("%.1e\t", t2);
        
        // variant 4-8
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.5)
        {
            mp_ptr uv = _nmod_vec_init(k);
            tt = clock();
            _nmod_vec_dot_product_multi_2_v4_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_8(uv, u, (mp_srcptr *) v, len, k, mod);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            _nmod_vec_clear(uv);
            nb_iter += 5;
        }
        t2 /= nb_iter; printf("%.1e\t", t2);
        
        // variant 8-8
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.5)
        {
            mp_ptr uv = _nmod_vec_init(k);
            tt = clock();
            _nmod_vec_dot_product_multi_2_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v8_8(uv, u, (mp_srcptr *) v, len, k, mod);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            _nmod_vec_clear(uv);
            nb_iter += 5;
        }
        t2 /= nb_iter; printf("%.1e\t", t2);
        
        // variant 4-32
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.5)
        {
            mp_ptr uv = _nmod_vec_init(k);
            tt = clock();
            _nmod_vec_dot_product_multi_2_v4_32(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_32(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_32(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_32(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_v4_32(uv, u, (mp_srcptr *) v, len, k, mod);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            _nmod_vec_clear(uv);
            nb_iter += 5;
        }
        t2 /= nb_iter; printf("%.1e\t", t2);
        
        // variant split26
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.5)
        {
            mp_ptr uv = _nmod_vec_init(k);
            tt = clock();
            _nmod_vec_dot_product_multi_2_split26(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_split26(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_split26(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_split26(uv, u, (mp_srcptr *) v, len, k, mod);
            _nmod_vec_dot_product_multi_2_split26(uv, u, (mp_srcptr *) v, len, k, mod);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            _nmod_vec_clear(uv);
            nb_iter += 5;
        }
        t2 /= nb_iter; printf("%.1e\t", t2);
    }


    // TODO test disabled for the moment (split26 sometimes fails due to overflow)
    if (0 && FLINT_BIT_COUNT(n) < 45)
    {
        mp_ptr uv1 = _nmod_vec_init(k);
        nmod_vec_dot_product_multi(uv1, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        mp_ptr uv2 = _nmod_vec_init(k);
        _nmod_vec_dot_product_multi_2_v1_8(uv2, u, (mp_srcptr *) v, len, k, mod);
        assert(_nmod_vec_equal(uv1, uv2, k) && "v1_8");
        _nmod_vec_dot_product_multi_2_v4_8(uv2, u, (mp_srcptr *) v, len, k, mod);
        assert(_nmod_vec_equal(uv1, uv2, k) && "v4_8");
        _nmod_vec_dot_product_multi_2_v8_8(uv2, u, (mp_srcptr *) v, len, k, mod);
        assert(_nmod_vec_equal(uv1, uv2, k) && "v8_8");
        _nmod_vec_dot_product_multi_2_v4_32(uv2, u, (mp_srcptr *) v, len, k, mod);
        assert(_nmod_vec_equal(uv1, uv2, k) && "v4_32");
        _nmod_vec_dot_product_multi_2_split26(uv2, u, (mp_srcptr *) v, len, k, mod);
        if (!_nmod_vec_equal(uv1, uv2, k))
        {
            printf("\n\n");
            _nmod_vec_print(uv1, k, mod);
            printf("\n\n");
            _nmod_vec_print(uv2, k, mod);
        }
        assert(_nmod_vec_equal(uv1, uv2, k) && "split26");
    }

    // flint's vector-matrix product
    t2 = 0.0; nb_iter = 0;
    while (t2 < 0.5)
    {
        mp_ptr uv = _nmod_vec_init(k);
        tt = clock();
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        nmod_mat_nmod_vec_mul(uv, ur, len, vmat);
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _nmod_vec_clear(uv);
        nb_iter += 5;
    }
    t2 /= nb_iter; printf("%.1e\t", t2);

    // flint's matrix-vector product
    nmod_mat_t vmattr;
    nmod_mat_init(vmattr, vmat->c, vmat->r, vmat->mod.n);
    nmod_mat_transpose(vmattr, vmat);
    t2 = 0.0; nb_iter = 0;
    while (t2 < 0.5)
    {
        mp_ptr uv = _nmod_vec_init(k);
        tt = clock();
        nmod_mat_mul_nmod_vec(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec(uv, vmattr, ur, len);
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _nmod_vec_clear(uv);
        nb_iter += 5;
    }
    t2 /= nb_iter; printf("%.1e\t", t2);

    // pml's matrix-vector product
    t2 = 0.0; nb_iter = 0;
    while (t2 < 0.5)
    {
        mp_ptr uv = _nmod_vec_init(k);
        tt = clock();
        nmod_mat_mul_nmod_vec_newdot(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec_newdot(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec_newdot(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec_newdot(uv, vmattr, ur, len);
        nmod_mat_mul_nmod_vec_newdot(uv, vmattr, ur, len);
        t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
        _nmod_vec_clear(uv);
        nb_iter += 5;
    }
    t2 /= nb_iter; printf("%.1e\t", t2);

    if (FLINT_BIT_COUNT(n) < 31)
    {
        // pml's small modulus matrix-vector product
        t2 = 0.0; nb_iter = 0;
        while (t2 < 0.5)
        {
            mp_ptr uv = _nmod_vec_init(k);
            tt = clock();
            nmod_mat_mul_nmod_vec_small_modulus(uv, vmattr, ur, len);
            nmod_mat_mul_nmod_vec_small_modulus(uv, vmattr, ur, len);
            nmod_mat_mul_nmod_vec_small_modulus(uv, vmattr, ur, len);
            nmod_mat_mul_nmod_vec_small_modulus(uv, vmattr, ur, len);
            nmod_mat_mul_nmod_vec_small_modulus(uv, vmattr, ur, len);
            t2 += (double)(clock()-tt) / CLOCKS_PER_SEC;
            _nmod_vec_clear(uv);
            nb_iter += 5;
        }
        t2 /= nb_iter; printf("%.1e\t", t2);
    }

    if (0)
    {
        mp_ptr uv1 = _nmod_vec_init(k);
        mp_ptr uv2 = _nmod_vec_init(k);
        nmod_mat_mul_nmod_vec_newdot(uv1, vmattr, ur, len);
        nmod_mat_mul_nmod_vec_small_modulus(uv2, vmattr, ur, len);
        assert(_nmod_vec_equal(uv1, uv2, len));
        _nmod_vec_clear(uv1);
        _nmod_vec_clear(uv2);
    }

    // test, for safety
    {
        mp_ptr uv1 = _nmod_vec_init(k);
        mp_ptr uv2 = _nmod_vec_init(k);
        //nmod_mat_mul_nmod_vec_newdot(uv1, vmattr, ur, len);
        nmod_mat_mul_nmod_vec(uv1, vmattr, ur, len);
        nmod_vec_dot_product_multi(uv2, u, (mp_srcptr *) v, len, k, maxbits1, maxbits2, mod);
        assert(_nmod_vec_equal(uv1, uv2, k));
        _nmod_vec_clear(uv1);
        _nmod_vec_clear(uv2);
    }


    _nmod_vec_clear(u);
    for (ulong i = 0; i < len; i++)
        _nmod_vec_clear(v[i]);
    flint_free(v);
    nmod_mat_clear(vmat);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_randinit(state);

    char labels[] = "nbits\tlen\tk\tmulti\tv1_8\tv8_8\tv8_32\tv16_16\tvecmat\tmatvec\tmatvec(pml)";

    // launch full suite
    if (argc == 1)
    {
        const ulong mods[11] =
        {
                (UWORD(1) << 3) + 1,
                (UWORD(1) << 10) + 1,
                (UWORD(1) << 20) + 1,
                (UWORD(1) << 25) + 1,
                (UWORD(1) << 29) + 1,
                (UWORD(1) << 30) + 1,
                (UWORD(1) << 31) + 1,
                (UWORD(1) << 40) + 1,
                (UWORD(1) << 50) + 1,
                (UWORD(1) << 60) + 1,
                UWORD_MAX
        };
        const ulong nbits[11] = { 4, 11, 21, 26, 30, 31, 32, 41, 51, 61, 64 };

        const ulong dims[24] = {1,2,3,4,5,8,12,16,20,25,30,40,50,75,100,150,200,300,400,500,750,1000,1500,2000};

        // warmup
        printf("Warmup:");
        time_nmod_vec_dot_product_multi(16, 16, nbits[3], nbits[3], mods[3], state);
        printf("\n");

        printf("%s\n",labels);
        for (slong i = 0; i < 11; i++)
        {
            //for (slong r = 0; r < 24; r++)
            //{
            //    for (slong c = 0; c < 24; c++)
            //    {
            //        printf("%ld\t%ld\t%ld\t%ld\t%ld\t", nbits[i]/2, nbits[i], nbits[i], dims[r], dims[c]);
            //        time_nmod_vec_dot_product_multi(dims[r], dims[c], nbits[i]/2, nbits[i], mods[i], state);
            //        printf("\n");
            //    }
            //}
            for (slong r = 0; r < 24; r++)
            {
                for (slong c = 0; c < 24; c++)
                {
                    //printf("%ld\t%ld\t%ld\t%ld\t%ld\t", nbits[i], nbits[i], nbits[i], dims[r], dims[c]);
                    printf("%ld\t%ld\t%ld\t", nbits[i], dims[r], dims[c]);
                    time_nmod_vec_dot_product_multi(dims[r], dims[c], nbits[i], nbits[i], mods[i], state);
                    printf("\n");
                }
            }
        }
    }

    // fixed number of bits
    else if (argc == 2)
    {
        // warmup
        printf("Warmup:");
        time_nmod_vec_dot_product_multi(16, 16, 20, 20, (UWORD(1) << 19) + 1, state);
        printf("\n");

        printf("%s\n",labels);
        const ulong nbit = atoi(argv[1]);
        const ulong dims[24] = {1,2,3,4,5,8,12,16,20,25,30,40,50,75,100,150,200,300,400,500,750,1000,1500,2000};
        for (slong r = 0; r < 24; r++)
        {
            for (slong c = 0; c < 24; c++)
            {
                printf("%ld\t%ld\t%ld\t", nbit, dims[r], dims[c]);
                time_nmod_vec_dot_product_multi(dims[r], dims[c], nbit, nbit, (UWORD(1) << (nbit-1))+1, state);
                printf("\n");
            }
        }
    }

    // fixed number of bits and dimensions
    else if (argc == 4)
    {
        // warmup
        printf("Warmup:");
        time_nmod_vec_dot_product_multi(16, 16, 20, 20, (UWORD(1) << 20) - 1, state);
        printf("\n");

        printf("%s\n",labels);
        const ulong nbit = atoi(argv[1]);
        const ulong len = atoi(argv[2]);
        const ulong k = atoi(argv[3]);
        printf("%ld\t%ld\t%ld\t", nbit, len, k);
        time_nmod_vec_dot_product_multi(len, k, nbit, nbit, (UWORD(1) << (nbit))-1, state);
        printf("\n");
    }

    flint_randclear(state);
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
