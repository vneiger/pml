#include <stdlib.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_extra.h"

int main(int argc, char ** argv)
{
    if (0)
    {
        nmod_mat_t mat;
        nmod_mat_init(mat, 3, 4, 7);

        mat->rows[0][0] = 2;
        mat->rows[0][1] = 4;
        mat->rows[0][2] = 5;
        mat->rows[0][3] = 1;
        mat->rows[1][0] = 2;
        mat->rows[1][1] = 3;
        //mat->rows[1][2] = 2; --> rank==2
        mat->rows[1][2] = 5;
        mat->rows[1][3] = 6;
        mat->rows[2][0] = 0;
        mat->rows[2][1] = 1;
        mat->rows[2][2] = 3;
        mat->rows[2][3] = 2;

        nmod_mat_print_pretty(mat);

        slong * P = _perm_init(mat->r);
        slong * Q = _perm_init(mat->c);
        nmod_mat_pluq(mat, P, Q);

        nmod_mat_print_pretty(mat);
    }

    if (1)
    {
        srand(time(NULL));
        flint_rand_t state;
        flint_randinit(state);
        flint_randseed(state, rand(), rand());

        nmod_mat_t mat;
        //nmod_mat_init(mat, atoi(argv[1]), atoi(argv[2]), 8388617);  // 24 bits
        nmod_mat_init(mat, atoi(argv[1]), atoi(argv[2]), 1099511627791);  // 41 bits
        //nmod_mat_init(mat, atoi(argv[1]), atoi(argv[2]), 36028797018963971);  // 56 bits
        //nmod_mat_init(mat, atoi(argv[1]), atoi(argv[2]), 2305843009213693967);  // 62 bits

        double t;
        clock_t tt;
        long nb_iter;

        // warmup
        t = 0.0;
        nb_iter = 0;
        while (t < 2)
        {
            nmod_mat_rand(mat, state);
            slong * P = _perm_init(mat->r);
            slong * Q = _perm_init(mat->c);
            tt = clock();
            nmod_mat_pluq(mat, P, Q);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 1;
        }
        t /= nb_iter;
        printf("warmup done\n");

        printf("new pluq\tlu_classical\tlu_delayed\tlu_recursive\n");

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            nmod_mat_rand(mat, state);
            slong * P = _perm_init(mat->r);
            slong * Q = _perm_init(mat->c);
            tt = clock();
            nmod_mat_pluq(mat, P, Q);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 1;
        }
        t /= nb_iter;
        printf("%4e\t", t);

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            nmod_mat_rand(mat, state);
            slong * P = _perm_init(mat->r);
            tt = clock();
            nmod_mat_lu_classical(P, mat, 0);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 1;
        }
        t /= nb_iter;
        printf("%4e\t", t);

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            nmod_mat_rand(mat, state);
            slong * P = _perm_init(mat->r);
            tt = clock();
            nmod_mat_lu_classical_delayed(P, mat, 0);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 1;
        }
        t /= nb_iter;
        printf("%4e\t", t);

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            nmod_mat_rand(mat, state);
            slong * P = _perm_init(mat->r);
            tt = clock();
            nmod_mat_lu_recursive(P, mat, 0);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 1;
        }
        t /= nb_iter;
        printf("%4e\n", t);
    }

    return 0;
}
