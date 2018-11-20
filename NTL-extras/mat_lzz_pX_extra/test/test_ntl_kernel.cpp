#include <NTL/mat_lzz_p.h>
#include <vector>

using namespace NTL;

int main(int argc, char * argv[])
{
    if (argc != 4)
        throw std::invalid_argument("Usage: ./test_ntl_kernel rdim cdim bitsize");
    long m = atoi(argv[1]);
    long n = atoi(argv[2]);
	zz_p::init(GenPrime_long(atoi(argv[3])));
    //zz_p::init(3);

	double t1=0.0, t2=0.0;
	double tt1, tt2;
	long nb_iter=0;

	while (nb_iter<10000000)
	//while (t1<1)
	{
		Mat<zz_p> mat;
		random(mat,m,n);

		tt1 = GetWallTime();
		Mat<zz_p> ker;
		kernel(ker,mat);
		tt2 = GetWallTime();
		t1 += tt2 - tt1;

        //Mat<zz_p> im;
        //image(im, mat);
        //std::cout << "Matrix: " << std::endl;
        //std::cout << mat << std::endl;
        //std::cout << "Kernel: " << std::endl;
        //std::cout << ker << std::endl;
        //std::cout << "Image: " << std::endl;
        //std::cout << im << std::endl;
        //image(im, ker);
        //std::cout << "Image of kernel: " << std::endl;
        //std::cout << im << std::endl;

		tt1 = GetWallTime();
        long k = ker.NumRows();
        std::vector<long> pcol(k,m-1);
        Vec<zz_p> pval;
        pval.SetLength(k);
        bool flag = false;

		for (long i = 0; i < k; ++i)
        {
            while (IsZero(ker[i][pcol[i]]))
                --pcol[i];
            pval[i] = ker[i][pcol[i]];

            //--> this does not happen ==> if pivot is at expected location, it is 1
            if (pval[i] != zz_p(1) and pcol[i] == m-k+i)
            {
                std::cout << i << "," << pcol[i] << "," << pval[i] << std::endl;
                flag = true;
            }

            //--> this does not happen ==> if pivot is at expected location, the rest of column is zero
            if (pcol[i] == m-k+i)
            {
                for (long ii = 0; ii < k; ++ii)
                {
                    if (ii != i)
                    {
                        if (not IsZero(ker[ii][pcol[i]]))
                        {
                            flag = true;
                        }
                    }
                }
                if (flag)
                    std::cout << i << "," << pcol[i] << "," << pval[i] << std::endl;
            }
        }
        if (flag)
        {
            std::cout << ker << std::endl;
            std::cout << mat << std::endl;
        }
		tt2 = GetWallTime();
		t2 += tt2-tt1;

		++nb_iter;
	}

	std::cout << zz_p::modulus() << std::endl;
	std::cout << nb_iter << std::endl;
	std::cout << (t1/nb_iter) << "," << (t2/nb_iter) << std::endl;
	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
