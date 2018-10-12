#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"


NTL_CLIENT

/*------------------------------------------------------------*/
/* checks one prime                                           */
/*------------------------------------------------------------*/
void check(long p)
{
    zz_p::init(p);

    vector<long> sizes = {10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 200, 220, 240, 260, 280, 300};
    vector<long> degrees = {10, 50, 100, 150, 200, 250, 300, 350, 400};
    const double thres = 0.01;

    if (NumBits(p) <= SMALL_PRIME_SIZE)
        cout << "#define MATRIX_LMULTIPLIER_DEGREE_THRESHOLD_SMALL ";
    else
        cout << "#define MATRIX_LMULTIPLIER_DEGREE_THRESHOLD_LARGE ";

    long done = 0;
    size_t i;
    long sz;
    for (i = 0; done == 0 && i < sizes.size(); i++)
    {
        double t_3, t_g;
        size_t j;
        long nb_3 = 0;
        sz = sizes[i];
        for (j = 0; j < degrees.size(); j++)
        {
            Mat<zz_pX> a, b, c;
            long deg, nb;

            deg = degrees[j];
            random(a, sz, sz, deg);
            
            mat_lzz_pX_lmultiplier_geometric mulg;
            mulg = mat_lzz_pX_lmultiplier_geometric(a, deg-1);
            mat_lzz_pX_lmultiplier_3_primes mul3;
            mul3 = mat_lzz_pX_lmultiplier_3_primes(a, deg-1);
            
            random(b, sz, 1, deg);

            t_3 = get_time();
            nb = 0;
            do
            {
                mul3.multiply(c, b);
                nb++;
            }
            while ((get_time()-t_3) <= thres);
            t_3 = (get_time()-t_3) / nb;


            t_g = get_time();
            nb = 0;
            do
            {
                mulg.multiply(c, b);
                nb++;
            }
            while ((get_time()-t_g) <= thres);
            t_g = (get_time()-t_g) / nb;

            if (t_3 < 1.05 * t_g)
                nb_3++;
        }             
        if (nb_3 <= 4)
            done = 1;
    }
    cout << sz << endl; 
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup();
    check(23068673);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
