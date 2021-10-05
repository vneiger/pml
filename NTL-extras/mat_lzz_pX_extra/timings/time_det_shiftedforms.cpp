// timing charpoly of matrix in "shifted forms" [Pernet-Storjohann]
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iterator>
#include <fstream>

#include <NTL/BasicThreadPool.h>
#include <numeric>
#include <string>

#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_forms.h"
#include "util.h"
//#include "mat_lzz_pX_forms.h"
//#include "mat_lzz_pX_utils.h"
//#include "mat_lzz_pX_determinant.h"
//#define GENERIC_DETSALL_PROFILE
//#define GENERIC_DETSONE_PROFILE
//#define SLOW

NTL_CLIENT

void conv_cdeg_uniquemult(std::vector<long> & unique_cdeg, std::vector<long> & mult_cdeg, const Vec<long> & cdeg)
{
    long current = -1;
    for (long i = 0; i < cdeg.length(); ++i)
    {
        if (cdeg[i] != current)
        {
            current = cdeg[i];
            unique_cdeg.push_back(current);
            mult_cdeg.push_back(1);
        }
        else
            ++mult_cdeg.back();
    }
}

// TODO threshold to investigate
// TODO there are many copies in this algo... avoidable?
// TODO for initial low degrees, better exploit Vec<Mat<zz_p>>
// TODO what about starting with another set of columns, which is more numerous
//      than the first one? this would better reduce the matrix size... still
//      unclear whether it would be faster...
// TODO at some point (e.g. with n=11 example) we manipulate rectangular wide pmbasis,
//      this will make the kernel degree grow
//      --> it might be a good idea to limit cdim1 to cdim/2
//      --> or actually to always take cdim/2, but the computed thing
//      will be a kernel of first columns and approx of the rest?
bool determinant_shifted_form_degaware_updateall(zz_pX & det, const Mat<zz_pX> & pmat, const VecLong & mult_cdeg, VecLong & diag_deg, VecLong & shifted_rdeg, long index, long threshold, long target_degdet, char prefix=' ')
{
#ifdef GENERIC_DETSALL_PROFILE
    double t;
#endif // GENERIC_DETSALL_PROFILE

    // det: output determinant
    // pmat: input square matrix, assumed in "shiftedform" (X^{s} Id + R, cdeg(R) < s)
    // --> degree of determinant is sum(cdeg) = sum of diagonal degrees
    // mult_cdeg: list of lengths, sum should be pmat row (and column) dimension
    // diag_deg: list of diagonal degrees of pmat (= its s-pivot-degree)
    // shifted_rdeg: row degree shifted by some "s"
    //     (initially, row degree 0, and all along, s = -initial column degree)
    // --> the matrix pmat remains s-owP all along, and shifted_rdeg gives its s-row deg

    const long dim = pmat.NumRows();

    //std::cout << "shifted_rdeg: ";
    //std::copy(shifted_rdeg.begin(), shifted_rdeg.end(), std::ostream_iterator<long>(std::cout, "\t"));
    //std::cout << std::endl;

    // retrieve the row degree and the determinant degree
    const long degdet = std::accumulate(diag_deg.begin(), diag_deg.end(), 0);

    // if degdet is not target_degree, then something went wrong in earlier steps
    // if they are equal, then we can proceed with the recursive calls
    if (degdet != target_degdet)
        return false;

    // if small dim, just run expansion by minors
    if (dim<=4)
    {
#ifdef GENERIC_DETSALL_PROFILE
        t=GetWallTime();
#endif // GENERIC_DETSALL_PROFILE
        determinant_expansion_by_minors(det, pmat);
#ifdef GENERIC_DETSALL_PROFILE
        std::cout << prefix << "\tbase case --> " << GetWallTime()-t << std::endl;
#endif // GENERIC_DETSALL_PROFILE
        return true;
    }

    // column dimension of the matrix that we will compute the kernel of
    long cdim1;

    // above some threshold (or if just one mult_cdeg left),
    // just run the usual algo splitting column dimension in two equal parts
    if (index >= std::min<long>(threshold,mult_cdeg.size()-1))
    {
#ifdef GENERIC_DETSALL_PROFILE
        std::cout << "\t-->Entering halving stage" << std::endl;
        prefix = '\t';
#endif
        cdim1 = (dim>>1);  // cdim1 ~ dim/2
    }
    // otherwise, handle leftmost mult_cdeg[index] columns and recurse with the rest
    else cdim1 = mult_cdeg[index];

    const long cdim2 = dim-cdim1;

    Mat<zz_pX> pmat_l;
    Mat<zz_pX> pmat_r;
    pmat_l.SetDims(dim,cdim1);
    pmat_r.SetDims(dim,cdim2);

    for (long i = 0; i < dim; ++i)
    {
        for (long j = 0; j < cdim1; ++j)
            pmat_l[i][j] = pmat[i][j];
        for (long j = 0; j < cdim2; ++j)
            pmat_r[i][j] = pmat[i][j+cdim1];
    }

    // compute the kernel via approximant basis at high order
    Mat<zz_pX> appbas;
    // degree of kernel basis will be (generically)  D = cdim1 * deg(pmat_l) / (dim - cdim1)
    // --> compute approximants at order deg(pmat_l) + D + 1
    // (cf for example Neiger-Rosenkilde-Solomatov ISSAC 2018, Lemma 4.3)
    const long degdet_ker = std::accumulate(diag_deg.begin(), diag_deg.begin()+cdim1, 0);
    const long deg_ker = ceil( degdet_ker / (double)(dim-cdim1) );
    const long order = *std::max_element(diag_deg.begin(), diag_deg.begin()+cdim1) + deg_ker + 1;

    VecLong shift(shifted_rdeg);
#ifdef GENERIC_DETSALL_PROFILE
    t = GetWallTime();
#endif // GENERIC_DETSALL_PROFILE
    pmbasis(appbas, pmat_l, order, shifted_rdeg);
#ifdef GENERIC_DETSALL_PROFILE
    t = GetWallTime()-t;
    std::cout << prefix << "\tpmbasis dims x order = " << dim << " x " << cdim1 << " x " << order << " || time " << t << std::endl;
#endif // GENERIC_DETSALL_PROFILE

    // minimal left kernel basis of pmat_r : last rows of app
    Mat<zz_pX> kerbas;
    kerbas.SetDims(cdim2,dim);
    for (long i = 0; i < cdim2; ++i)
        for (long j = 0; j < dim; ++j)
            kerbas[i][j] = appbas[i+cdim1][j];

    // update shifted_rdeg
    // shifted_rdeg-row degree of kerbas
    // = last entries of shift
    // = s-row degree of product kerbas*pmat_r (see description of function, for "s")
    std::vector<long>(shifted_rdeg.begin()+cdim1, shifted_rdeg.end()).swap(shifted_rdeg);

    // update diag_deg (keeping only the end of it)
    std::vector<long>(diag_deg.begin()+cdim1, diag_deg.end()).swap(diag_deg);
    for (long k = 0; k < cdim2; ++k)
        diag_deg[k] += shifted_rdeg[k] - shift[cdim1+k];

    //if (dim < 60)
    //{
    //    std::cout << degree_matrix(pmat_l) << std::endl;
    //    std::cout << degree_matrix(kerbas) << std::endl;
    //}

    // then compute the product
    Mat<zz_pX> pmatt;
#ifdef GENERIC_DETSALL_PROFILE
    t = GetWallTime();
#endif // GENERIC_DETSALL_PROFILE
    multiply(pmatt, kerbas, pmat_r);
#ifdef GENERIC_DETSALL_PROFILE
    t = GetWallTime()-t;
    std::cout << prefix << "\tmultiply degrees " << deg(kerbas) << "," << deg(pmat_r) << " || time " << t << std::endl;
    const long actual_deg_ker = deg(kerbas);
    std::cout << prefix << "\tker deg, actual deg: " << deg_ker << "\t" << actual_deg_ker << std::endl;
#endif // GENERIC_DETSALL_PROFILE

    return determinant_shifted_form_degaware_updateall(det,pmatt,mult_cdeg,diag_deg,shifted_rdeg,index+1,threshold,target_degdet,prefix);
}

// TODO threshold to investigate
// TODO there are many copies in this algo... avoidable?
// TODO for initial low degrees, better exploit Vec<Mat<zz_p>>
// TODO what about starting with another set of columns, which is more numerous
//      than the first one? this would better reduce the matrix size... still
//      unclear whether it would be faster...
// TODO at some point (e.g. with n=11 example) we manipulate rectangular wide pmbasis,
//      this will make the kernel degree grow
//      --> it might be a good idea to limit cdim1 to cdim/2
//      --> or actually to always take cdim/2, but the computed thing
//      will be a kernel of first columns and approx of the rest?
bool determinant_shifted_form_degaware_updateone(zz_pX & det, Mat<zz_pX> & kerbas, const Mat<zz_pX> & pmat, const VecLong & mult_cdeg, VecLong & diag_deg, VecLong & shifted_rdeg, long index, long threshold, long target_degdet, char prefix=' ')
{
#ifdef GENERIC_DETSONE_PROFILE
    double t;
#endif // GENERIC_DETSONE_PROFILE

    // det: output determinant
    // kerbas: intermediate kernel basis, stored to multiply remaining columns
    //      (initially, must be an "empty" matrix with 0x0 dimension)
    // pmat: input matrix,
    //     the matrix we want the determinant of is in fact kerbas*pmat,
    //     which is assumed in "shiftedform" (X^{s} Id + R, cdeg(R) < s)
    //     --> degree of determinant is sum(cdeg) = sum of diagonal degrees
    // mult_cdeg: list of lengths, sum should be pmat row (and column) dimension
    // shifted_rdeg: row degree shifted by some "s"
    //     (initially, row degree 0, and all along, s = -initial column degree)
    // --> the matrix pmat remains s-owP all along, and shifted_rdeg gives its s-row deg
    // --> 0-row-degree of kerbas is shifted_rdeg
    // diag_deg: diagonal degree of current matrix (the actual one we want the determinant of)

    const long rdim = pmat.NumRows();
    const long cdim = pmat.NumCols();
    const bool kerbas_empty = (kerbas.NumRows() == 0);

    //std::cout << rdim << " x " << cdim << std::endl;
    //std::cout << kerbas.NumRows() << " x " << kerbas.NumCols() << std::endl;
    //std::cout << pmat.NumRows() << " x " << pmat.NumCols() << std::endl;

    //std::cout << "shifted_rdeg: ";
    //std::copy(shifted_rdeg.begin(), shifted_rdeg.end(), std::ostream_iterator<long>(std::cout, "\t"));
    //std::cout << std::endl;

    // retrieve the row degree and the determinant degree
    const long degdet = std::accumulate(diag_deg.begin(), diag_deg.end(), 0);

    // if degdet is not target_degree, then something went wrong in earlier steps
    // if they are equal, then we can proceed with the recursive calls
    if (degdet != target_degdet)
        return false;

    // if small dim, just run expansion by minors
    if (cdim<=4)
    {
#ifdef GENERIC_DETSONE_PROFILE
        t=GetWallTime();
#endif // GENERIC_DETSONE_PROFILE
        if (not kerbas_empty)
        {
            Mat<zz_pX> pmatt;
            multiply(pmatt, kerbas, pmat);
            determinant_expansion_by_minors(det, pmatt);
        }
        else
            determinant_expansion_by_minors(det, pmat);
#ifdef GENERIC_DETSONE_PROFILE
        std::cout << prefix << "\tbase case --> " << GetWallTime()-t << std::endl;
#endif // GENERIC_DETSONE_PROFILE
        return true;
    }

    // column dimension of the matrix that we will compute the kernel of
    long cdim1;

    // above some threshold (or if just one mult_cdeg left),
    // just run the usual algo splitting column dimension in two equal parts
    const bool halving = (index >= std::min<long>(threshold,mult_cdeg.size()-1));
    if (halving && not kerbas_empty)
    {
        // kerbas nonempty => first time we entered the halving stage
        // --> recurse on pmatt = kerbas * pmat
        // --> make kerbas empty
        Mat<zz_pX> pmatt;
#ifdef GENERIC_DETSONE_PROFILE
        std::cout << "\t-->Entering halving stage" << std::endl;
        prefix = '\t';
        t = GetWallTime();
#endif // GENERIC_DETSONE_PROFILE
        multiply(pmatt, kerbas, pmat);
#ifdef GENERIC_DETSONE_PROFILE
        t = GetWallTime()-t;
        std::cout << prefix << "\tmultiply degrees " << deg(kerbas) << "," << deg(pmat) << " || time " << t << std::endl;
#endif // GENERIC_DETSONE_PROFILE
        kerbas.kill();
        return determinant_shifted_form_degaware_updateone(det,kerbas,pmatt,mult_cdeg,diag_deg,shifted_rdeg,index,threshold,target_degdet,prefix);
    }
    // if halving and kerbas is empty
    // --> either we directly entered halving from the initial input
    // --> or we have already applied one step of the halving stage 
    // in both case we do the "split in 2" divide and conquer
    // and in both cases at this point cdim=rdim and all other properties in
    // input of determinant update-all

    // if halving, take basically half dimension blocks
    // if not halving, handle leftmost mult_cdeg[index] columns and recurse with the rest
    if (halving)
        cdim1 = (cdim>>1);  // cdim1 ~ dim/2
    else
        cdim1 = mult_cdeg[index];

    const long cdim2 = cdim-cdim1;

    Mat<zz_pX> pmat_l;
    Mat<zz_pX> pmat_r;
    pmat_l.SetDims(rdim,cdim1);
    pmat_r.SetDims(rdim,cdim2);

    for (long i = 0; i < rdim; ++i)
    {
        for (long j = 0; j < cdim1; ++j)
            pmat_l[i][j] = pmat[i][j];
        for (long j = 0; j < cdim2; ++j)
            pmat_r[i][j] = pmat[i][j+cdim1];
    }

    // if not halving and not initial call, pmat_l must be multiplied by kerbas
    if (not halving && not kerbas_empty)
    {
#ifdef GENERIC_DETSONE_PROFILE
        std::cout << prefix << "\tkerbas*pmat_l multiply dimensions x degrees " << kerbas.NumRows() << "x" << rdim << "x" << cdim1 << "," << deg(kerbas) << "x" << deg(pmat_l) << " || time ";
        t = GetWallTime();
#endif // GENERIC_DETSONE_PROFILE
        multiply(pmat_l, kerbas, pmat_l);
#ifdef GENERIC_DETSONE_PROFILE
        t = GetWallTime()-t;
        std::cout << t << std::endl;
#endif // GENERIC_DETSONE_PROFILE
    }

    // compute the kernel via approximant basis at high order
    Mat<zz_pX> appbas;
    const long app_rdim = pmat_l.NumRows();
    // degree of kernel basis will be (generically)  D = cdim1 * deg(pmat_l) / (dim - cdim1)
    // --> compute approximants at order deg(pmat_l) + D + 1
    // (cf for example Neiger-Rosenkilde-Solomatov ISSAC 2018, Lemma 4.3)
    const long degdet_ker = std::accumulate(diag_deg.begin(), diag_deg.begin()+cdim1, 0);
    const long deg_ker = ceil( degdet_ker / (double)(app_rdim-cdim1) );
    const long order = *std::max_element(diag_deg.begin(), diag_deg.begin()+cdim1) + deg_ker + 1;

    // save shift to be able to update diag_deg
    VecLong shift(shifted_rdeg);
#ifdef GENERIC_DETSONE_PROFILE
    t = GetWallTime();
#endif // GENERIC_DETSONE_PROFILE
    pmbasis(appbas, pmat_l, order, shifted_rdeg);
#ifdef GENERIC_DETSONE_PROFILE
    t = GetWallTime()-t;
    std::cout << prefix << "\tpmbasis dims x order = " << app_rdim << " x " << cdim1 << " x " << order << " || time " << t << std::endl;
#endif // GENERIC_DETSONE_PROFILE

    // minimal left kernel basis of pmat_r : last rows of app
    Mat<zz_pX> kerbas2;
    kerbas2.SetDims(cdim2,app_rdim);
    for (long i = 0; i < cdim2; ++i)
        for (long j = 0; j < app_rdim; ++j)
            kerbas2[i][j] = appbas[i+cdim1][j];

    // update shifted_rdeg
    // shifted_rdeg-row degree of kerbas
    // = last entries of shift
    // = s-row degree of product kerbas*pmat_r (see description of function, for "s")
    std::vector<long>(shifted_rdeg.begin()+cdim1, shifted_rdeg.end()).swap(shifted_rdeg);

    // update diag_deg (keeping only the end of it)
    std::vector<long>(diag_deg.begin()+cdim1, diag_deg.end()).swap(diag_deg);
    for (long k = 0; k < cdim2; ++k)
        diag_deg[k] += shifted_rdeg[k] - shift[cdim1+k];

    //if (dim < 60)
    //{
    //    std::cout << degree_matrix(pmat_l) << std::endl;
    //    std::cout << degree_matrix(kerbas) << std::endl;
    //}

    // if halving: compute kerbas2*pmat_r
    if (halving)
    {
#ifdef GENERIC_DETSONE_PROFILE
        std::cout << prefix << "\tmultiply degrees " << deg(kerbas2) << "," << deg(pmat_r) << " || time ";
        t = GetWallTime();
#endif // GENERIC_DETSONE_PROFILE
        multiply(pmat_r, kerbas2, pmat_r);
#ifdef GENERIC_DETSONE_PROFILE
        t = GetWallTime()-t;
        std::cout << t << std::endl;
        const long actual_deg_ker = deg(kerbas2);
        std::cout << prefix << "\tker deg, actual deg: " << deg_ker << "\t" << actual_deg_ker << std::endl;
#endif // GENERIC_DETSONE_PROFILE
    }
    // else if non-initial non-halving step: update kerbas by multiplication
    else if (kerbas.NumRows() > 0)
    {
#ifdef GENERIC_DETSONE_PROFILE
        std::cout << prefix << "\tkerbas multiply degrees " << deg(kerbas2) << "x" << deg(kerbas) << " || time ";
        t = GetWallTime();
#endif // GENERIC_DETSONE_PROFILE
        multiply(kerbas, kerbas2, kerbas);
#ifdef GENERIC_DETSONE_PROFILE
        t = GetWallTime()-t;
        std::cout << t << std::endl;
        const long actual_deg_ker = deg(kerbas2);
        std::cout << prefix << "\tnew ker deg, actual deg: " << deg_ker << "\t" << actual_deg_ker << std::endl;
#endif // GENERIC_DETSONE_PROFILE
    }
    // else, initial non-halving step: kerbas is empty, meaning identity
    // --> make it equal to kerbas2
    else
        kerbas.swap(kerbas2);

    // call recursively
    return determinant_shifted_form_degaware_updateone(det,kerbas,pmat_r,mult_cdeg,diag_deg,shifted_rdeg,index+1,threshold,target_degdet,prefix);
}


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

    std::vector<long> unique_cdeg;
    std::vector<long> unique_cdeg2;
    std::vector<long> mult_cdeg;
    std::vector<long> mult_cdeg2;
    conv_cdeg_uniquemult(unique_cdeg,mult_cdeg,cdeg);
    conv_cdeg_uniquemult(unique_cdeg2,mult_cdeg2,cdeg2);
    //std::cout << "sequence of degrees:\t";
    //std::copy(unique_cdeg2.begin(), unique_cdeg2.end(), std::ostream_iterator<long>(std::cout, "\t"));
    //std::cout << std::endl;
    //std::cout << "ncols for each degree:\t";
    //std::copy(mult_cdeg2.begin(), mult_cdeg2.end(), std::ostream_iterator<long>(std::cout, "\t"));
    //std::cout << std::endl;
    //std::copy(unique_cdeg2.begin(), unique_cdeg2.end(), std::ostream_iterator<long>(std::cout, "\t"));
    //std::cout << std::endl;
    //std::copy(mult_cdeg2.begin(), mult_cdeg2.end(), std::ostream_iterator<long>(std::cout, "\t"));
    //std::cout << std::endl;

    //std::cout << dmat << std::endl;
    //std::cout << cdeg << std::endl;
    //std::copy(unique_cdeg.begin(), unique_cdeg.end(), std::ostream_iterator<long>(std::cout, ""));
    //std::cout << std::endl;
    //std::copy(mult_cdeg.begin(), mult_cdeg.end(), std::ostream_iterator<long>(std::cout, ""));
    //std::cout << std::endl;

    //std::cout << dmat2 << std::endl;
    //std::cout << cdeg2 << std::endl;
    //std::copy(unique_cdeg2.begin(), unique_cdeg2.end(), std::ostream_iterator<long>(std::cout, ""));
    //std::cout << std::endl;
    //std::copy(mult_cdeg2.begin(), mult_cdeg2.end(), std::ostream_iterator<long>(std::cout, ""));
    //std::cout << std::endl;

    //std::cout << "Matrix dimension: " << dim << ", degree of determinant: " << degdet << std::endl;

    std::vector<double> timings;
    double t,tt;
    long nb_iter;

#ifdef SLOW
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
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in naive triangular approach" << std::endl;
    }
#endif

#ifdef GENERIC_DETSALL_PROFILE
    std::cout << std::endl;
#endif

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
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in naive triangular(mirror) approach" << std::endl;
    }

#ifdef GENERIC_DETSALL_PROFILE
    std::cout << std::endl;
#endif

#ifdef SLOW
    { // shifted form specific, degree-aware, update all
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat);
            VecLong diag_deg = col_degree(pmat);
            VecLong shifted_rdeg(dim);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_shifted_form_degaware_updateall(det, pmat, mult_cdeg, diag_deg, shifted_rdeg, 0, 1, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in degree-aware triangular approach" << std::endl;
    }
#endif

#ifdef GENERIC_DETSALL_PROFILE
    std::cout << std::endl;
#endif

    for (long thres=1; thres < 8; ++thres)
    { // shifted form specific, degree-aware, update all
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat;
            random(pmat, dmat2);
            VecLong diag_deg = col_degree(pmat);
            VecLong shifted_rdeg(dim);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_shifted_form_degaware_updateall(det, pmat, mult_cdeg2, diag_deg, shifted_rdeg, 0, thres, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in degree-aware-all triangular(mirror) approach" << std::endl;
    }

#ifdef GENERIC_DETSONE_PROFILE
    std::cout << std::endl;
#endif

    for (long thres=1; thres < 8; ++thres)
    { // shifted form specific, degree-aware, update one
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            Mat<zz_pX> pmat,kerbas;
            random(pmat, dmat2);
            VecLong diag_deg = col_degree(pmat);
            VecLong shifted_rdeg(dim);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_shifted_form_degaware_updateone(det, kerbas, pmat, mult_cdeg2, diag_deg, shifted_rdeg, 0, thres, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in degree-aware-one triangular(mirror) approach" << std::endl;
    }

#ifdef GENERIC_DETSONE_PROFILE
    std::cout << std::endl;
#endif

#ifdef SLOW
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
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in linsolve approach" << std::endl;
    }
#endif

    //std::cout << nthreads << "\t" << fftprime << "\t" << nbits << "\t" << dim << "\t" << degdet << "\t";
    std::cout << nbits << "\t" << dim << "\t" << degdet << "\t" << std::setprecision(3) << (double)dim/degdet*100 << std::setprecision(8) << "\t";
    for (auto timing : timings)
        std::cout << timing << "\t";

    std::cout << std::endl;
}

/*---------------------------*/
/* runs benchs on all files  */
/*---------------------------*/
void run_bench()
{
    //std::vector<long> nthreads = {1,2,4,8,16};
    std::vector<long> nthreads = {1};
    //std::vector<long> fftprimes = {0,1};
    std::vector<long> fftprimes = {0};
    //std::vector<long> nbits = {23,28,31,60};
    std::vector<long> nbits = {28,60};
    std::vector<const char*> filenames = {
        "degree-pattern-random-2-6.txt",
        "degree-pattern-random-2-7.txt",
        "degree-pattern-random-2-8.txt",
        "degree-pattern-random-2-9.txt",
        "degree-pattern-random-2-10.txt",
        "degree-pattern-random-2-11.txt",
        "degree-pattern-random-2-12.txt",
        "degree-pattern-random-2-13.txt",
        "degree-pattern-random-2-14.txt",
        "degree-pattern-random-3-4.txt",
        "degree-pattern-random-3-5.txt",
        "degree-pattern-random-3-6.txt",
        "degree-pattern-random-3-7.txt",
        "degree-pattern-random-3-8.txt",
        "degree-pattern-random-3-9.txt",
        "degree-pattern-random-4-4.txt",
        "degree-pattern-random-4-5.txt",
        "degree-pattern-random-4-6.txt",
        "degree-pattern-random-4-7.txt",
        "degree-pattern-random-5-3.txt",
        "degree-pattern-random-5-4.txt",
        "degree-pattern-random-5-5.txt",
        "degree-pattern-random-5-6.txt",
    };
    std::vector<const char*> params = {
        "2-6",
        "2-7",
        "2-8",
        "2-9",
        "2-10",
        "2-11",
        "2-12",
        "2-13",
        "2-14",
        "3-4",
        "3-5",
        "3-6",
        "3-7",
        "3-8",
        "3-9",
        "4-4",
        "4-5",
        "4-6",
        "4-7",
        "5-3",
        "5-4",
        "5-5",
        "5-6",
    };
    for (long nthread : nthreads)
        for (long nbit : nbits)
            for (long fftprime : fftprimes)
            {
                long n = 6;
                for (size_t i=0; i < filenames.size(); ++i)
                {
                    std::cout << params[i] << "\t";
                    run_one_bench(nthread,fftprime,nbit,filenames[i]);
                    ++n;
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

    //std::vector<string> labels = {"naivetri","naivetri-mirror","linsolve"};
    std::vector<string> labels = {"naivetri-mirror","degaware-all-1","degaware-all-2","degaware-all-3","degaware-all-4","degaware-all-5","degaware-all-6","degaware-all-7","degaware-one-1","degaware-all-2","degaware-all-3","degaware-all-4","degaware-all-5","degaware-all-6","degaware-all-7"};

    if (argc == 1)
    {
        //cout << "n-deg\tthreads\tfftp\tnbits\tdim\tdegdet\t";
        cout << "n-deg\tnbits\tdim\tdegdet\tratio\t";
        for (auto lab : labels)
            std::cout << lab << "\t";
        std::cout << std::endl;
        warmup();
        run_bench();
    }

    else if (argc!= 4 and argc!=5)
        throw std::invalid_argument("Usage: ./time_det_shiftedforms [nbits fftprime degreefile [nthreads]]");

    else
    {
        cout << "threads\tfftp\tnbits\tdim\tdegdet\t";
        for (auto lab : labels)
            std::cout << lab << "\t";
        std::cout << std::endl;
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
