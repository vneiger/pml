#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

vector<long> sizes = {5, 10, 50, 100};
vector<long> degrees = {50, 150};

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
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    // cout << "p=" << p << endl;
    const double thres = 0.01;

    long middle = 0;
    long geometric = 0;
    long FFT = 0;
    long num_runs = 0;

    for (size_t i = 0; i < sizes.size(); i++)
    {
        long sz = sizes[i];
        size_t j;
        for (j = 0; j < degrees.size(); j++)
        {
            num_runs++;
            long deg = degrees[j];
            
            Mat<zz_pX> a, x;
            Mat<zz_p> a0;
            do
            {
                random_mat_zz_pX(a, sz, sz, deg);
                GetCoeff(a0, a, 0);
            }
            while (determinant(a0) == 0);
        
            // cout << "size = " << sz << " deg = " << deg << endl;
            if (p != 0)
            {
                vector<double> vec_middle;
                for (long r = 1; r < 8 && (1L << r) < deg; r++)
                {
                    double t_middle = get_time();
                    long nb = 0;
                    do
                    {
                        newton_inv_trunc_middle_product(x, a, deg, r);
                        nb++;
                    }
                    while ((get_time()-t_middle) <= thres);
                    t_middle = (get_time()-t_middle) / nb;
                    vec_middle.push_back(t_middle);
                }                
                long min_middle = 1 + index_min(vec_middle);
                middle += min_middle;
                
                vector<double> vec_geometric;
                for (long r = 1; r < 8 && (1L << r) < deg; r++)
                {
                    double t_geometric = get_time();
                    long nb = 0;
                    do
                    {
                        newton_inv_trunc_geometric(x, a, deg, r);
                        nb++;
                    }
                    while ((get_time()-t_geometric) <= thres);
                    t_geometric = (get_time()-t_geometric) / nb;
                    vec_geometric.push_back(t_geometric);
                }
                long min_geometric = 1 + index_min(vec_geometric);
                geometric += min_geometric;
                // cout << min_middle << " ";
                // cout << min_geometric << " ";
                // cout << endl;
            }
            else
            {
                vector<double> vec_FFT;
                for (long r = 1; r < 8 && (1L << r) < deg; r++)
                {
                    double t_FFT = get_time();
                    long nb = 0;
                    do
                    {
                        newton_inv_trunc_FFT(x, a, deg, r);
                        nb++;
                    }
                    while ((get_time()-t_FFT) <= thres);
                    t_FFT = (get_time()-t_FFT) / nb;
                    vec_FFT.push_back(t_FFT);
                }
                long min_FFT = 1 + index_min(vec_FFT);
                FFT += min_FFT;
                // cout << min_FFT << " ";
                // cout << endl;
            }
        }
    }


    if (p == 0)
    {
        cout << "#define MATRIX_INV_TRUNC_PLAIN_THRESHOLD_FFT " << (long) (FFT/num_runs) << endl;
    }
    else
    {
        if (type_of_prime() == TYPE_SMALL_PRIME)
        {
            cout << "#define MATRIX_INV_TRUNC_PLAIN_THRESHOLD_MIDDLE_SMALL " << (long) (middle/num_runs) << endl;
            cout << "#define MATRIX_INV_TRUNC_PLAIN_THRESHOLD_GEOMETRIC_SMALL " << (long) (geometric/num_runs) << endl;
        }
        else
        {
            cout << "#define MATRIX_INV_TRUNC_PLAIN_THRESHOLD_MIDDLE_LARGE " << (long) (middle/num_runs) << endl;
            cout << "#define MATRIX_INV_TRUNC_PLAIN_THRESHOLD_GEOMETRIC_LARGE " << (long) (geometric/num_runs) << endl;
        }
    }
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
