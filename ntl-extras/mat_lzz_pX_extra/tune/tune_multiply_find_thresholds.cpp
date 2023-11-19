#include <iomanip>
#include <limits.h>

#include "util.h"
#include "mat_lzz_pX_multiply.h"

vector<long> sizes = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250};

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    const double thres = 0.1;

    // first, find the crossover between evaluate and 3 primes 
    // need it only for non-FFT primes
    vector<long> degrees = {1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250};
    vector<long> thresholds;

    if (p != 0)
    {
        if (NumBits(p) <= SMALL_PRIME_SIZE)
            cout << "static long MATRIX_DEGREE_THRESHOLDS_SMALL[" << sizes.size() << "] = {";
        else
            cout << "static long MATRIX_DEGREE_THRESHOLDS_LARGE[" << sizes.size() << "] = {";

        for (size_t i = 0; i < sizes.size(); i++)
        {
            long done = 0;
            long nb_in_a_row = 0;
            long sz = sizes[i];
            size_t j;

            for (j = 0; done == 0 && j < degrees.size(); j++)
            {
                long deg, nb;
                Mat<zz_pX> a, b, c;
                double t_eval, t_3primes;
                deg = degrees[j];
                random(a, sz, sz, deg);
                random(b, sz, sz, deg);

                t_eval = get_time();
                nb = 0;
                do
                {
                    multiply_evaluate_dense(c, a, b);
                    nb++;
                }
                while ((get_time()-t_eval) <= thres);
                t_eval = (get_time()-t_eval) / nb;

                t_3primes = get_time();
                nb = 0;
                do
                {
                    multiply_3_primes(c, a, b);
                    nb++;
                }
                while ((get_time()-t_3primes) <= thres);
                t_3primes = (get_time()-t_3primes) / nb;

                if (t_3primes < t_eval)
                    nb_in_a_row++;
                else 
                    nb_in_a_row = 0;

                if (nb_in_a_row == 3)
                    done = 1;
            }

            if (nb_in_a_row == 0)
                cout << LONG_MAX;
            else
                cout << degrees[j - nb_in_a_row];

            if (i < sizes.size() - 1)
                cout << ", ";
            thresholds.push_back(degrees[j]);
        }
        cout << "};\n";
    }

    if (p == 0)
    {
        cout << "static long MATRIX_WAKSMAN_THRESHOLDS_FFT[" << sizes.size() << "] = {";

        // find degree threshold for waksman 
        for (size_t i = 0; i < sizes.size(); i++)
        {
            double t_fft, t_waksman;
            size_t j;
            bool done = false;
            long nb_in_a_row = 0;
            long sz = sizes[i];
            for (j = 0; (not done) && j < degrees.size(); j++)
            {
                Mat<zz_pX> a, b, c;
                long deg, nb;

                deg = degrees[j];
                random(a, sz, sz, deg);
                random(b, sz, sz, deg);

                t_fft = get_time();
                nb = 0;
                do
                {
                    multiply_evaluate_FFT(c, a, b);
                    nb++;
                }
                while ((get_time()-t_fft) <= thres);
                t_fft = (get_time()-t_fft) / nb;

                t_waksman = get_time();
                nb = 0;
                do
                {
                    multiply_waksman(c, a, b);
                    nb++;
                }
                while ((get_time()-t_waksman) <= thres);
                t_waksman = (get_time()-t_waksman) / nb;

                if (t_waksman > t_fft)
                    nb_in_a_row++;
                else 
                    nb_in_a_row = 0;
                if (nb_in_a_row == 3)
                    done = true;
            }             
            if (nb_in_a_row == 0)
                cout << LONG_MAX;
            else
                cout << degrees[j - nb_in_a_row];

            if (i < sizes.size() - 1)
                cout << ", ";
        }
    }
    else 
    {
        if (NumBits(p) <= SMALL_PRIME_SIZE)
            cout << "static long MATRIX_WAKSMAN_THRESHOLDS_SMALL[" << sizes.size() << "] = {";
        else
            cout << "static long MATRIX_WAKSMAN_THRESHOLDS_LARGE[" << sizes.size() << "] = {";

        // find degree threshold for waksman 
        for (size_t i = 0; i < sizes.size(); i++)
        {
            double t_old, t_waksman;
            size_t j;
            long done = 0;
            long nb_in_a_row = 0;
            long sz = sizes[i];
            for (j = 0; done == 0 && j < degrees.size(); j++)
            {
                Mat<zz_pX> a, b, c;
                long deg, nb;

                deg = degrees[j];
                random(a, sz, sz, deg);
                random(b, sz, sz, deg);

                if (deg <= thresholds[i])
                {
                    t_old = get_time();
                    nb = 0;
                    do
                    {
                        multiply_evaluate_dense(c, a, b);
                        nb++;
                    }
                    while ((get_time()-t_old) <= thres);
                    t_old = (get_time()-t_old) / nb;
                }
                else
                {
                    t_old = get_time();
                    nb = 0;
                    do
                    {
                        multiply_3_primes(c, a, b);
                        nb++;
                    }
                    while ((get_time()-t_old) <= thres);
                    t_old = (get_time()-t_old) / nb;
                }

                t_waksman = get_time();
                nb = 0;
                do
                {
                    multiply_waksman(c, a, b);
                    nb++;
                }
                while ((get_time()-t_waksman) <= thres);
                t_waksman = (get_time()-t_waksman) / nb;

                if (t_waksman > t_old)
                    nb_in_a_row++;
                else 
                    nb_in_a_row = 0;
                if (nb_in_a_row == 3)
                    done = 1;
            }             
            if (nb_in_a_row == 0)
                cout << LONG_MAX;
            else
                cout << degrees[j - nb_in_a_row];

            if (i < sizes.size() - 1)
                cout << ", ";
        }
    }
    cout << "};\n";
}




/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    cout << "#define MATRIX_THRESHOLDS_LEN " << sizes.size() << endl;
    cout << "static long MATRIX_THRESHOLDS_SIZES[" << sizes.size() << "] = {";
    for (size_t i = 0; i < sizes.size(); i++)
    {
        cout << sizes[i];
        if (i < sizes.size() - 1)
            cout << ", ";
    }
    cout << "};\n";

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
