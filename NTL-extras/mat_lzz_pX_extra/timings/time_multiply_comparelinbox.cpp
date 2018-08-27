#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_bench_fft(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    double t;

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    if (sz<=4)
    {
        t = GetWallTime();
        multiply_waksman(c1, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (naive)" << endl;

        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_FFT(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (evaluate)" << endl;
    }
    else
    {
        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_FFT(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << endl;
    }
}

/*------------------------------------------------------------*/
void one_bench_3primes(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;
    double t;
		
    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    if (sz<=4)
    {
        t = GetWallTime();
        multiply_waksman(c1, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (naive)" << endl;

        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_3_primes(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (3 primes)" << endl;
    }
    else
    {
        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_3_primes(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << endl;
    }
}

/*------------------------------------------------------------*/
void one_bench_geometric(long sz, long deg)
{
    Mat<zz_pX> a, b, c0, c1, c2, c4, c5;
    double t;

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    if (sz<=4)
    {
        t = GetWallTime();
        multiply_waksman(c1, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (naive)" << endl;

        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_geometric(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << "  (geometric)" << endl;
    }
    else
    {
        // evaluation -- should be done only if feasible
        t = GetWallTime();
        multiply_evaluate_geometric(c2, a, b);
        t = GetWallTime()-t;
        cout << sz << "," << deg << "," << t << endl;
    }
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long test, long nbits)
{
		std::vector<long> szs =
				{
            2,2,2,2,2,2,2,2,2,2,2,2,
            4,4,4,4,4,4,4,4,4,4,4,4,
            8,8,8,8,8,8,8,8,8,8,8,8,
            16,16,16,16,16,16,16,16,16,16,16,16,
            32,32,32,32,32,32,32,32,32,32,32,
            64,64,64,64,64,64,64,64,64,
            128,128,128,128,128,128,128,
            256,256,256,256,256,
            512,512,512,
            1024,1024,
				};
		std::vector<long> degs =
				{
            32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
            32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
            32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
            32,64,128,256,512,1024,2048,4096,8192,16384,32768,131072,
            32,64,128,256,512,1024,2048,4096,8192,16384,32768,
            32,64,128,256,512,1024,2048,4096,8192,
            32,64,128,256,512,1024,2048,
            32,64,128,256,512,
            32,64,128,
            32,64,
				};

		if (test==0 || test==1)
		{
			std::cout << "Bench polynomial matrix multiplication (FFT prime)" << std::endl;
			if (nbits < 25)
			{
				zz_p::UserFFTInit(786433); // 20 bits
				cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 20 << ")" << endl;
			}
			else if (nbits < 35)
			{
				zz_p::UserFFTInit(2013265921); // 31 bits
				cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 31 << ")" << endl;
			}
			else if (nbits < 45)
			{
				zz_p::UserFFTInit(2748779069441); // 42 bits
				cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 42 << ")" << endl;
			}
			else if (nbits < 65)
			{
				zz_p::UserFFTInit(1139410705724735489); // 60 bits
				cout << "p = " << zz_p::modulus() << "  (FFT prime, bit length = " << 60 << ")" << endl;
			}
			for (size_t i=0;i<szs.size();i++)
				one_bench_fft(szs[i],degs[i]);
			cout << endl;
		}

		if (test==0 || test==2)
		{
			long prime = NTL::GenPrime_long(nbits);
			zz_p::init(prime);
			std::cout << "Bench polynomial matrix multiplication (3 primes)" << std::endl;
			cout << "p = " << zz_p::modulus() << "  (normal prime, bit length = " << nbits << ")" << endl;
			for (size_t i=0;i<szs.size();i++)
				one_bench_3primes(szs[i],degs[i]);
			cout << endl;
		}
		
		if (test==0 || test==3)
		{
			long prime = NTL::GenPrime_long(nbits);
			zz_p::init(prime);
			std::cout << "Bench polynomial matrix multiplication (geometric)" << std::endl;
			cout << "p = " << zz_p::modulus() << "  (normal prime, bit length = " << nbits << ")" << endl;
			for (size_t i=0;i<szs.size();i++)
				one_bench_geometric(szs[i],degs[i]);
			cout << endl;
		}

}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    long test = 0; // default: run all benchs
    long nbits = 60;
    
    if (argc>=2)
        nbits = atoi(argv[1]);
    if (argc==3)
        test = atoi(argv[2]);
    if (argc>3)
        throw std::invalid_argument("Usage: ./time_multiply_comparelinbox OR ./time_multiply_comparelinbox nbits OR ./time_multiply_comparelinbox nbits test");
    
    warmup();
    run_bench(test,nbits);

    return 0;
}
