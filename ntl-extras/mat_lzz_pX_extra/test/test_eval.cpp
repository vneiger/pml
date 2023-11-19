#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "mat_lzz_pX_extra.h"

/*
* how does it compare, in terms of speed, evaluation when points are random, and geometric
* when base field is FFT prime, how does it compare: random / geometric / FFTpoints
* when evaluating a bunch of polynomials at the same points, in the polynomial evaluation code (at least for random & geometric points), is there nothing that could be gathered as pre-computation, done once for all polynomials and then apply the rest to each?
* Mat.elts?
*/

NTL_CLIENT

int main(int argc, char *argv[]){
    long nbits = 60;
    long deg = 5000;
    long npts = 5000;
    
    if (argc==4)
    {
        deg = atoi(argv[1]);
        npts = atoi(argv[2]);
        nbits = atoi(argv[3]);
    }else{
        cout << "usage: deg, npts, nbits" << endl;
        return 1;
    }
    
    zz_p::init(NTL::GenPrime_long(nbits));
    cout << "Evaluating over " <<
            nbits << " bit prime, degree " << deg <<
            ", and pts " << npts << endl;
    
    zz_pX rand_pol;
    random(rand_pol, deg+1);
    Vec<zz_p> val;
    
    double t1w, t2w;
    
    // random points
    {
        auto eval_general = get_general_points(npts);
        t1w = GetWallTime();
        eval_general.evaluate(val, rand_pol);
        t2w = GetWallTime();
        cout << "Time (rand): " << t2w - t1w << endl;
    
        // geometric points
        zz_p r;
        random(r);
        zz_pX_Multipoint_Geometric eval_geo(r,npts);
        t1w = GetWallTime();
        eval_geo.evaluate(val, rand_pol);
        t2w = GetWallTime();
        cout << "Time (geo): " << t2w-t1w << endl;
    }
    // FFT prime
    zz_p::FFTInit(0);
    cout << "\nFFT prime" << endl;
    {
        auto eval_general = get_general_points(npts);
        t1w = GetWallTime();
        eval_general.evaluate(val, rand_pol);
        t2w = GetWallTime();
        cout << "Time (rand): " << t2w - t1w << endl;
    
        // geometric points
        zz_p r;
        random(r);
        zz_pX_Multipoint_Geometric eval_geo(r,npts);
        t1w = GetWallTime();
        eval_geo.evaluate(val, rand_pol);
        t2w = GetWallTime();
        cout << "Time (geo): " << t2w-t1w << endl;
        
        eval_geo.set_FFT_evaluate();
        t1w = GetWallTime();
        eval_geo.evaluate(val, rand_pol);
        t2w = GetWallTime();
        cout << "Time (geo FFT): " << t2w-t1w << endl;
    }
    
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
