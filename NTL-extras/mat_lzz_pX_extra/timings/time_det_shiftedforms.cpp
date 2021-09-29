// timing charpoly of matrix in "shifted forms" [Pernet-Storjohann]
#include <functional>
#include <iomanip>
#include <fstream>

#include <NTL/BasicThreadPool.h>
#include <numeric>

#include "mat_lzz_pX_forms.h"
#include "util.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_determinant.h"

NTL_CLIENT


// stores the degree matrix and column degree
void retrieve_degree_matrix(Mat<long> & degree_matrix, Vec<long> & cdeg, const char *filename)
{
    // open file and check it went fine
    ifstream file; file.open(filename);
    if(not file.good()) {
        cerr << "Error in opening file stream" << endl;
        return;
    }

    // read matrix dimension and set up degree matrix
    long dim;
    file >> dim;
    degree_matrix.SetDims(dim,dim);

    long i=0;
    long j=0;
    while ( file >> degree_matrix[i][j] )
    {
        ++degree_matrix[i][j];
        if (j < dim-1)
            ++j;
        else
        {
            ++i;
            j=0;
        }
    }
    file.close();

    // deduce column degree and current leading positions
    cdeg.SetLength(dim,-1);
    Vec<long> lpos(INIT_SIZE,dim);
    for (long j = 0; j < dim; ++j)
        for (long i = 0; i < dim; ++i)
            if (degree_matrix[i][j]-1 > cdeg[j])
            {
                cdeg[j] = degree_matrix[i][j]-1;
                lpos[j] = i;
            }

    // permute rows of degree matrix to make it diagonal-dominant
    // (leading positions will be on the diagonal)
    Mat<long> tmp_mat(degree_matrix);
    for (long i = 0; i < dim; ++i)
        degree_matrix[i] = tmp_mat[lpos[i]];

    //std::cout << tmp_mat << std::endl;
    //std::cout << degree_matrix << std::endl;
}

/*---------------------------------*/
/* runs one bench on a given file  */
/*---------------------------------*/
void run_one_bench(long nthreads, bool fftprime, long nbits, const char* filename)
{
    SetNumThreads(nthreads);

    // set base field
    if (fftprime)
    {
        //cout << "Bench determinant of shifted form, FFT prime p = ";
        if (nbits <= 20)
        {
            zz_p::UserFFTInit(786433); // 20 bits
            nbits=20;
            //cout << zz_p::modulus() << ", bit length = " << 20 << endl;
        }
        else if (nbits <= 31)
        {
            zz_p::UserFFTInit(2013265921); // 31 bits
            nbits=31;
            //cout << zz_p::modulus() << ", bit length = " << 31 << endl;
        }
        else if (nbits <= 42)
        {
            zz_p::UserFFTInit(2748779069441); // 42 bits
            nbits=42;
            //cout << zz_p::modulus() << ", bit length = " << 42 << endl;
        }
        else if (nbits <= 60)
        {
            zz_p::UserFFTInit(1139410705724735489); // 60 bits
            nbits=60;
            //cout << zz_p::modulus() << ", bit length = " << 60 << endl;
        }
        else
        {
            std::cout << "FFT prime with more than 60 bits is not supported" << std::endl;
            return;
        }
    }
    else
    {
        //cout << "Bench determinant of shifted form, random prime p = ";
        zz_p::init(NTL::GenPrime_long(nbits));
        //cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    // retrieve degree matrix from file (this one has decreasing diag degrees)
    Mat<long> dmat;
    Vec<long> cdeg;
    retrieve_degree_matrix(dmat,cdeg,filename);
    const long dim = dmat.NumRows();
    const long degdet = std::accumulate(cdeg.begin(),cdeg.end(),0);

    // compute mirrored degree matrix (increasing diag degrees)
    Mat<long> dmat2(INIT_SIZE, dim, dim);
    Vec<long> cdeg2(INIT_SIZE, dim);
    for (long i = 0; i < dim; ++i)
    {
        cdeg2[i] = cdeg[dim-1-i];
        for (long j = 0; j < dim; ++j)
            dmat2[i][j] = dmat[dim-1-i][dim-1-j];
    }
    //std::cout << dmat2 << std::endl;
    //std::cout << cdeg2 << std::endl;

    //std::cout << "Matrix dimension: " << dim << ", degree of determinant: " << degdet << std::endl;

    double t,tt,t_naivetri,t_naivetri2,t_linsolve,t_linsolve2;
    long nb_iter;

    { // generic case
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_generic_knowing_degree(det, pmat, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        t_naivetri = t/nb_iter;
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in naive triangular approach" << std::endl;
    }

    std::cout << std::endl;

    { // generic case, mirrored
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat2);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_generic_knowing_degree(det, pmat, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        t_naivetri2 = t/nb_iter;
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in naive triangular(mirror) approach" << std::endl;
    }

    std::cout << std::endl;

    { // via random linear system
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_linsolve(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        t_linsolve = t/nb_iter;
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in linsolve approach" << std::endl;
    }

    { // via random linear system, mirrored
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat2);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_linsolve(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        t_linsolve2 = t/nb_iter;
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in linsolve(mirror) approach" << std::endl;
    }


    // the three following methods are very slow (unsurprisingly) for low degree matrices
    if (false)
    { // via evaluation, general points
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_evaluation_general(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(ev-general):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (false)
    { // via evaluation, geometric points
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_evaluation_geometric(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(ev-geometric):\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }

    if (false)
    { // via evaluation, FFT points
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat);
            tt = GetWallTime();
            zz_pX det;
            determinant_via_evaluation_FFT(det, pmat);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = verify_determinant(det, pmat, true, true);
        }
        std::cout << "Time(ev-FFT):\t\t" << t/nb_iter << (ok ? "\t(ok)":"  (notok)") << std::endl;
    }
    std::cout << nthreads << "\t" << fftprime << "\t" << nbits << "\t" << dim << "\t" << degdet << "\t" << t_naivetri << "\t" << t_naivetri2 << "\t" << t_linsolve << "\t" << t_linsolve2 << std::endl;
}

/*---------------------------*/
/* runs benchs on all files  */
/*---------------------------*/
void run_bench()
{
    std::vector<long> nthreads = {1,2,3,4,8,12,16,24,32};
    std::vector<long> fftprimes = {0,1};
    std::vector<long> nbits = {20,31,42,60};
    for (long nthread : nthreads)
        for (long fftprime : fftprimes)
            for (long nbit : nbits)
            {
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-6.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-7.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-8.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-9.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-10.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-11.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-12.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-13.txt");
                run_one_bench(nthread,fftprime,nbit,"degree-pattern-random-2-14.txt");
            }
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    if (argc == 1)
    {
        cout << "threads\tfftp\tnbits\tdim\tdegdet\tnaivetri\tnaivetri-mirror\tlinsolve\tlinsolve-mirror\t" << endl;
        warmup();
        run_bench();
    }

    else if (argc!= 4 and argc!=5)
        throw std::invalid_argument("Usage: ./time_det_shiftedforms [nbits fftprime degreefile [nthreads]]");

    else
    {
        cout << "threads\tfftp\tnbits\tdim\tdegdet\tnaivetri\tnaivetri-mirror\tlinsolve\tlinsolve-mirror\t" << endl;
        const long nbits = atoi(argv[1]);
        const bool fftprime = (atoi(argv[2])==1);
        long nthreads=1;
        if (argc == 5)
            nthreads=atoi(argv[4]);

        warmup();
        run_one_bench(nthreads,fftprime,nbits,argv[3]);
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
