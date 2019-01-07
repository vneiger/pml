#include <vector>
#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "util.h"
#include "vec_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* runs timings up to size 1000                               */
/*------------------------------------------------------------*/
void bench()
{
    zz_p::FFTInit(0);
    long p = zz_p::modulus();

    std::vector<long> lengths = {0, 1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    std::cout << "prime\t\tlength\t\tfast\t\tnaive" << std::endl;

    for (long i : lengths)
    {
        vec_zz_p A, invA1, invA2;
        random(A, i);

        cout << p << "\t";
        cout << i << "\t";
        double t;
        long nb;
        const double thresh = 0.1;

        nb = 0;
        t = get_time();
        do
        {
            inv(invA1, A);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << "\t";

        nb = 0;
        t = get_time();
        do
        {
            inv_naive(invA1, A);
            nb++;
        }
        while ( (get_time() - t) < thresh);
        t = (get_time() - t) / nb;
        cout << t << "\t";

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
    bench();

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
