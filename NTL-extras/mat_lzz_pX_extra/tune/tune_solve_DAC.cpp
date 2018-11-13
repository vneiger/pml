#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* index of the smallest element in v                         */
/*------------------------------------------------------------*/
long index_min(const vector<double>& v)
{
    long idx = 0;
    for (size_t i = 1; i < v.size(); i++)
    {
        if (v[i] < v[idx])
            idx = i;
    }
    return idx;
}

/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    VecLong szs =
        {
            1, 10, 20, 30, 40, 50, 100
        };

    VecLong degs =
        {
            100, 150, 200, 250, 300, 350, 400
        };

    long sum = 0;
    size_t start = 3;
    for (size_t si = 0; si < szs.size(); si++)
    {
        long sz = szs[si];
        // cout << sz << " ";
        double best = 0;
        double len = 0;
        for (size_t di = 0; di < degs.size(); di++)
        {
            long deg = degs[di];
            Mat<zz_pX> A, b, u, res;
            const double thresh = 0.01;
            random(A, sz, sz, deg);
            random(b, sz, 1, deg);
	    
            vector<double> vec;
            for (long i = (1L << start); i <= 64; i = 2*i)
            {
                double t;
                long nb;
                t = get_time();
                nb = 0;
                do
                {
                    solve_series_low_precision(u, A, b, deg, i);
                    nb++;
                }
                while ((get_time()-t) <= thresh);
                t = (get_time()-t) / nb;
                vec.push_back(t);
            }
            long min_t = start + index_min(vec);
            best += min_t;
            len += 1;
        }
        sum += (int) round(((double)best) / len);
    }
    long t = type_of_prime();

    cout << "#define ";

    switch(t)
    {
    case TYPE_FFT_PRIME:
        cout << "THRESHOLDS_SOLVE_LOW_PRECISION_FFT ";
        break;
    case TYPE_SMALL_PRIME:
        cout << "THRESHOLDS_SOLVE_LOW_PRECISION_SMALL ";
        break;
    case TYPE_LARGE_PRIME:
        cout << "THRESHOLDS_SOLVE_LOW_PRECISION_LARGE ";
        break;
    default:
        LogicError("Unknown prime type in linear solving.");
    }
    cout << (1L << (int) round (((double)sum) / szs.size())) << endl;
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup();
    check(0);
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
