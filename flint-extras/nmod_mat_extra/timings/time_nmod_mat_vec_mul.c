#include <flint/nmod_vec.h>
#include <time.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_extra.h"
#include "nmod_vec_extra.h"

/*--------------------------------------------------------------*/
/* computes a matrix-vector product in size r x c modulo n      */
/*--------------------------------------------------------------*/
void time_nmod_mat_vec_mul(slong r, slong c, ulong n)
{
    flint_rand_t state;
    flint_rand_init(state);

    nmod_t mod;
    nmod_init(&mod, n);

    double t;
    clock_t tt;
    long nb_iter;

    nmod_mat_t A;
    nmod_mat_init(A, r, c, mod.n);
    nmod_mat_rand(A, state);

    ulong * u = flint_malloc(c * sizeof(ulong));
    ulong * v = flint_malloc(r * sizeof(ulong));
    nmod_mat_t umat; nmod_mat_init(umat, c, 1, n);
    nmod_mat_t vmat1; nmod_mat_init(vmat1, r, 1, n);
    nmod_mat_t vmat2; nmod_mat_init(vmat2, 1, r, n);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        _nmod_vec_rand(u, state, c, mod);
        tt = clock();
        nmod_mat_mul_nmod_vec(v, A, u, c);
        nmod_mat_mul_nmod_vec(v, A, u, c);
        nmod_mat_mul_nmod_vec(v, A, u, c);
        nmod_mat_mul_nmod_vec(v, A, u, c);
        nmod_mat_mul_nmod_vec(v, A, u, c);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_rand(umat, state);
        tt = clock();
        nmod_mat_mul(vmat1, A, umat);
        nmod_mat_mul(vmat1, A, umat);
        nmod_mat_mul(vmat1, A, umat);
        nmod_mat_mul(vmat1, A, umat);
        nmod_mat_mul(vmat1, A, umat);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g\t", t);

    //t = 0.0;
    //nb_iter = 0;
    //while (t < 0.5)
    //{
    //    nmod_mat_rand(umat, state);
    //    tt = clock();
    //    nmod_mat_mul_blas(vmat1, A, umat);
    //    nmod_mat_mul_blas(vmat1, A, umat);
    //    nmod_mat_mul_blas(vmat1, A, umat);
    //    nmod_mat_mul_blas(vmat1, A, umat);
    //    nmod_mat_mul_blas(vmat1, A, umat);
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    nb_iter += 5;
    //}
    //t /= nb_iter;
    //printf("%4g\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        _nmod_vec_rand(v, state, c, mod);
        tt = clock();
        nmod_mat_nmod_vec_mul(u, v, r, A);
        nmod_mat_nmod_vec_mul(u, v, r, A);
        nmod_mat_nmod_vec_mul(u, v, r, A);
        nmod_mat_nmod_vec_mul(u, v, r, A);
        nmod_mat_nmod_vec_mul(u, v, r, A);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_rand(vmat2, state);
        tt = clock();
        nmod_mat_mul_classical(umat, vmat2, A);
        nmod_mat_mul_classical(umat, vmat2, A);
        nmod_mat_mul_classical(umat, vmat2, A);
        nmod_mat_mul_classical(umat, vmat2, A);
        nmod_mat_mul_classical(umat, vmat2, A);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g\t", t);

    //t = 0.0;
    //nb_iter = 0;
    //while (t < 0.5)
    //{
    //    nmod_mat_rand(vmat2, state);
    //    tt = clock();
    //    nmod_mat_mul_blas(umat, vmat2, A);
    //    nmod_mat_mul_blas(umat, vmat2, A);
    //    nmod_mat_mul_blas(umat, vmat2, A);
    //    nmod_mat_mul_blas(umat, vmat2, A);
    //    nmod_mat_mul_blas(umat, vmat2, A);
    //    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
    //    nb_iter += 5;
    //}
    //t /= nb_iter;
    //printf("%4g\t", t);

    printf("\n");

    nmod_mat_clear(A);
    nmod_mat_clear(umat);
    nmod_mat_clear(vmat1);
    nmod_mat_clear(vmat2);
    flint_free(u);
    flint_free(v);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main()
{
    //printf("nbits\tr\tc\tmul_nmod_vec\tmatrix mul\tmatrix mul_blas\tnmod_vec_mul\tmatrix mul\tmatrix mul_blas\n");
    printf("nbits\tr\tc\tmul_nmod_vec\tmatrix mul\tnmod_vec_mul\tmatrix mul\n");
    for (slong r = 1; r < 200; r += 5)
    {
        printf("%d\t%ld\t%ld\t", 30, r, r);
        time_nmod_mat_vec_mul(r, r, (1L << 29) + 1);
    }
    for (slong r = 1; r < 200; r += 5)
    {
        printf("%d\t%ld\t%ld\t", 61, r, r);
        time_nmod_mat_vec_mul(r, r, (1L << 60) + 1);
    }

    return 0;
}
