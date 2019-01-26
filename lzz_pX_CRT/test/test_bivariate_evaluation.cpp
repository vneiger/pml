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

NTL_CLIENT

int main(int argc, char *argv[]){
    long nbits = atoi(argv[1]);
    long dx = atoi(argv[2]);
    long dy = atoi(argv[3]);
    bool check = (atoi(argv[4])==1) ? true : false;

    long nbpoints = dx*dy;

    std::cout << "~~~BIVARIATE EVALUATION (semi-naive algo)~~~" << std::endl;

    if (nbits==0)
    {
        zz_p::FFTInit(0);
        std::cout << "Working over prime field, p = " << zz_p::modulus() << "  (FFT prime)" << std::endl;
    }
    else
    {
        zz_p::init(NTL::GenPrime_long(nbits));
        std::cout << "Working over prime field, p = " << zz_p::modulus() << std::endl;
    }

    // declare dx polynomials in y, of degree < dy
    // (they are the coefficients of f(x,y) in x)
    Vec<zz_pX> ycoeffs;
    ycoeffs.SetLength(dx);
    for (long i = 0; i < dx; ++i)
        random(ycoeffs[i], dy);

    // choose random points
    Vec<zz_p> xpoints = random_vec_zz_p(nbpoints);
    Vec<zz_p> ypoints = random_vec_zz_p(nbpoints);

    // build multipoint evaluators/interpolators
    zz_pX_Multipoint_General eval_x(xpoints);
    zz_pX_Multipoint_General eval_y(ypoints);

    // product of x-xpoints[i] for i in 1...nboints
    zz_pX modulus;

    double t1w, t2w, t1, t2;
    t1 = GetWallTime();

    // construct modulus and interpolant
    t1w = GetWallTime();
    BuildFromRoots(modulus, xpoints);
    t2w = GetWallTime();
    std::cout << "modulus: " << (t2w-t1w) << std::endl;

    // polynomial that will contain f(x,ell) mod g
    zz_pX pol;

    // compute all ycoeffs(lagrange(x)) by evaluation-interpolation
    Vec<zz_p> evals;
    zz_pX interpolant;
    t1w = GetWallTime();
    for (long i = 0; i < dx; ++i)
    {
        // evaluate ycoeffs[i] at all ypoints
        eval_y.evaluate(evals, ycoeffs[i]);
        // interpolate
        eval_x.interpolate(interpolant, evals);
        //random(interpolant, nbpoints); // FIXME
        // add to result
        interpolant <<= i;
        pol += interpolant;
    }
    t2w = GetWallTime();
    std::cout << "main loop: " << (t2w-t1w) << std::endl;

    t1w = GetWallTime();
    rem(pol, pol, modulus);
    t2w = GetWallTime();
    std::cout << "rem: " << (t2w-t1w) << std::endl;

    t1w = GetWallTime();
    eval_x.evaluate(evals, pol);
    t2w = GetWallTime();
    std::cout << "final eval: " << (t2w-t1w) << std::endl;

    t2 = GetWallTime();
    cout << "Time (via-eval): " << t2-t1 << endl;

    if (check)
    {
        t1w = GetWallTime();
        Vec<zz_p> evs;
        evs.SetLength(nbpoints);
        for (long k = 0; k < nbpoints; ++k)
            for (long i = 0; i < dx; ++i)
            {
                zz_p evi = eval(ycoeffs[i], ypoints[k]);
                evs[k] += evi * power(xpoints[k], i);
            }
        t2w = GetWallTime();

        std::cout << "time naive algo: " << (t2w-t1w) << std::endl;

        for (long k = 0; k < nbpoints; ++k)
            if (evs[k] != evals[k])
            {
                cout << "wrong!!" << endl;
                return 0;
            }
    }

    return 0;

}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
