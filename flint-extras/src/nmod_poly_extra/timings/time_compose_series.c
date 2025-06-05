#include <stdlib.h>

#include <flint/nmod_poly.h>
#include <flint/flint.h>
#include <time.h>

#include "nmod_poly_extra.h"

int main(int argc, char *argv[])
{
    const long p = 1048583;
    const long prec = atoi(argv[1]);

    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, time(NULL), time(NULL)+4);

    nmod_poly_t f;
    nmod_poly_init(f, p);
    nmod_poly_rand(f, state, prec);
    nmod_poly_t g;
    nmod_poly_init(g, p);
    nmod_poly_rand(g, state, prec);
    nmod_poly_set_coeff_ui(g, 0, 0);

    double t = 0.0;
    double tt;
    long nb_iter = 0;
    while (t < 1)
    {
        nmod_poly_t res;
        nmod_poly_init(res, p);
        tt = clock();
        nmod_poly_compose_series(res, f, g, prec);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter++;
    }
    printf("%ld -- %.1e\n", prec, t/nb_iter);
    return 0;
}
