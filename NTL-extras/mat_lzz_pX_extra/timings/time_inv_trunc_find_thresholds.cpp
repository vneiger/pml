#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

vector<long> sizes = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150};
vector<long> degrees = {50, 100, 150};

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    cout << "p=" << p << endl;
    const double thres = 0.01;

    for (size_t i = 0; i < sizes.size(); i++)
    {
	long sz = sizes[i];
	size_t j;
	for (j = 0; j < degrees.size(); j++)
	{
	    long deg = degrees[j];
	    
	    double t_middle, t_geometric, t_FFT;
	    long nb;
	    
	    Mat<zz_pX> a, x;
	    Mat<zz_p> a0;
	    do
	    {
		random_mat_zz_pX(a, sz, sz, deg);
		GetCoeff(a0, a, 0);
	    }
	    while (determinant(a0) == 0);
	    
	    cout << "size = " << sz << " deg = " << deg << endl;
	    for (long i = 1; i < 8 && (1L << i) < deg; i++)
	    {
		t_middle = get_time();
		nb = 0;
		do
		{
		    newton_inv_trunc_middle_product(x, a, deg, i);
		    nb++;
		}
		while ((get_time()-t_middle) <= thres);
		t_middle = (get_time()-t_middle) / nb;
		cout << t_middle << " ";
	    }
	    cout << endl;

 	    for (long i = 1; i < 8 && (1L << i) < deg; i++)
	    {
		t_geometric = get_time();
		nb = 0;
		do
		{
		    newton_inv_trunc_geometric(x, a, deg, i);
		    nb++;
		}
		while ((get_time()-t_geometric) <= thres);
		t_geometric = (get_time()-t_geometric) / nb;
		cout << t_geometric << " ";
	    }
	    cout << endl;

	    if (is_FFT_prime())
	    {
		for (long i = 1; i < 8 && (1L << i) < deg; i++)
		{
		    t_FFT = get_time();
		    nb = 0;
		    do
		    {
			newton_inv_trunc_FFT(x, a, deg, i);
			nb++;
		    }
		    while ((get_time()-t_FFT) <= thres);
		    t_FFT = (get_time()-t_FFT) / nb;
		    cout << t_FFT << " ";
		}
		cout << endl;
	    }
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
