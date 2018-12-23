#include <iomanip>

#include <NTL/BasicThreadPool.h>
#include "util.h"
#include "mat_lzz_pX_interpolant.h"

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_bench_pmbasis(long sz, long order)
{
    long rdim = sz*2;
    long cdim = sz;

    std::cout << rdim << "\t" << cdim << "\t" << order << "\t";

    VecLong shift(rdim,0);

    double t,tt;
    long nb_iter;

    // build random matrix
    Vec<Mat<zz_p>> evals;
    Vec<zz_p> pts;

    nb_iter=0; t=0.0;
    while (t<0.5)
    {
        evals.SetLength(order);
        for (long pt = 0; pt < order; ++pt)
            random(evals[pt], rdim, cdim);
        random(pts, order);
        tt = GetWallTime();
        Mat<zz_pX> intbas;
        pmbasis(intbas,evals,pts,shift,0,pts.length());
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << t/nb_iter << "\t";

    std::cout << std::endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void run_bench(long nbits)
{
    VecLong szs =
    {
        1,1,1,1,1,1,1,1,1,1,//1,1,
        2,2,2,2,2,2,2,2,2,2,//2,1,
        4,4,4,4,4,4,4,4,4,//4,4,
        8,8,8,8,8,8,8,8,8,//8,8,
        16,16,16,16,16,16,16,16,16,//16,16,
        32,32,32,32,32,32,32,32,32,//32,32,
        64,64,64,64,64,64,64,//64,64,
        128,128,128,128,128,//128,128,
        256,256,256,//256,256,
        512,//512,512,
        //1024,1024,
    };
    VecLong degs =
    {
        32,64,128,256,512,1024,2048,4096,8192,16384,//32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,16384,//32768,131072,
        32,64,128,256,512,1024,2048,4096,8192,//16384,32768,
        32,64,128,256,512,1024,2048,4096,8192,//16384,32768,
        32,64,128,256,512,1024,2048,4096,8192,//16384,32768,
        32,64,128,256,512,1024,2048,4096,//8192,16384,
        32,64,128,256,512,1024,2048,//4096,8192,
        32,64,128,256,512,//1024,2048,
        32,64,128,//256,512,
        32,//64,128,
        //32,64,
    };

    std::cout << "Bench pm-basis (FFT prime)" << std::endl;
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
    std::cout << "rdim\tcdim\torder\tpmbasis\t" << std::endl;
    for (size_t i=0;i<szs.size();++i)
        one_bench_pmbasis(szs[i],degs[i]);
    cout << endl;
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(1);

    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    if (argc==1)
    {
        warmup();
        run_bench(60);
    }
    if (argc==2)
    {
        warmup();
        run_bench(atoi(argv[1]));
    }
    if (argc>3)
        throw std::invalid_argument("Usage: ./time_pmbasis OR ./time_pmbasis nbits");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
