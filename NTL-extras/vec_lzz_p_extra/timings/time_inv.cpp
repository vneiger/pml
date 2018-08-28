#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "util.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* runs timings up to size 1000                               */
/*------------------------------------------------------------*/
void check(int opt)
{
    zz_p::FFTInit(0);
    long p = zz_p::modulus();

    for (long i = 0; i < 1000; i += 1)
    {
        vec_zz_p A, invA1, invA2;
        random_vec_zz_p(A, i);

        cout << p << " ";
        cout << i << " ";
        double t;
        long nb;
        const double thresh = 0.05;

        nb = 0;
        t = get_time();
        do
        {
            inv(invA1, A);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        nb = 0;
        t = get_time();
        do
        {
            inv_naive(invA1, A);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << " ";

        cout << endl;
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    int opt = 0;
    if (argc > 1)
        opt = atoi(argv[1]);
    check(opt);

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
