// timing charpoly of matrix in "shifted forms" [terminology by Pernet-Storjohann]
#include <NTL/ZZ.h>
#include <NTL/vec_lzz_p.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iterator>
#include <fstream>

#include <NTL/BasicThreadPool.h>
#include <NTL/mat_lzz_p.h>
#include <numeric>
#include <string>

#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_forms.h"
#include "mat_lzz_pX_kernel.h"
#include "mat_lzz_pX_utils.h"
#include "util.h"
//#include "mat_lzz_pX_forms.h"
//#include "mat_lzz_pX_utils.h"
//#include "mat_lzz_pX_determinant.h"
#define CACHE_FRIENDLY_SIZE 32
#define MATRIX_BLOCK_SIZE 8

//#define PROFILE_DEG_UPMAT
//#define PROFILE_DEG_UPKER
//#define PROFILE_SMART_UPMAT
//#define PROFILE_SMART_UPMAT_V2
//#define PROFILE_SMART_UPKER
//#define PROFILE_ZLS_AVG
//#define PROFILE_KER_AVG
//#define PROFILE_KERNEL_STEP1
//#define PROFILE_KERNEL_DEGREE1_CDEG
//#define PROFILE_KERNEL_STEP2
//#define PROFILE_KERNEL_DEGREE1
//#define DEBUGGING_NOW
//#define UNIFORM_KSHIFTED
//#define TIME_FMA

//#define TIME_LINSOLVE // via linear system solving with random rhs
//#define TIME_TRI_DECR // generic determinant on matrix with decreasing diagonal degrees
//#define TIME_TRI_INCR // generic determinant on matrix with increasing diagonal degrees
//#define TIME_DEG_UPMAT // going degree by degree, updating remaining matrix columns at each iteration
//#define TIME_DEG_UPKER // going degree by degree, updating kernel basis at each iteration
//#define TIME_SMART_UPMAT // going degree by degree, UPMAT, smart kernel computation
#define TIME_SMART_UPMAT_V2 // going degree by degree, UPMAT, smart kernel computation
//#define TIME_SMART_UPKER // going degree by degree, UPKER, smart kernel computation
//#define TIME_ZLS_AVG // ZLS style, but taking average degrees into account
//#define ESTIMATE_WIEDEMANN
//#define ESTIMATE_KELLERGEHRIG

#ifdef DEBUGGING_NOW
    static std::ostream &operator<<(std::ostream &out, const VecLong &s)
    {
        out << "[ ";
        for (auto &i: s)
            out << i << " ";
        return out << "]";
    }
#endif // DEBUGGING_NOW


NTL_CLIENT

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
#ifdef PROFILE_DEG_UPMAT
    double t;
#endif // PROFILE_DEG_UPMAT

    // det: output determinant
    // pmat: input square matrix, assumed in "shiftedform" (X^{s} Id + R, cdeg(R) < s)
    // --> degree of determinant is sum(cdeg) = sum of diagonal degrees
    // mult_cdeg: list of lengths, sum should be pmat row (and column) dimension
    // diag_deg: list of diagonal degrees of pmat (= its s-pivot-degree)
    // shifted_rdeg: row degree shifted by some "s"
    //     (initially, row degree 0, and all along, s = -initial column degree)
    // --> the matrix pmat remains s-owP all along, and shifted_rdeg gives its s-row deg

    const long dim = pmat.NumRows();

    // retrieve the row degree and the determinant degree
    const long degdet = std::accumulate(diag_deg.begin(), diag_deg.end(), 0);

    // if degdet is not target_degree, then something went wrong in earlier steps
    // if they are equal, then we can proceed with the recursive calls
    if (degdet != target_degdet)
        return false;

    // if small dim, just run expansion by minors
    if (dim<=4)
    {
#ifdef PROFILE_DEG_UPMAT
        t=GetWallTime();
#endif // PROFILE_DEG_UPMAT
        determinant_expansion_by_minors(det, pmat);
#ifdef PROFILE_DEG_UPMAT
        std::cout << prefix << "\tbase case --> " << GetWallTime()-t << std::endl;
#endif // PROFILE_DEG_UPMAT
        return true;
    }

    // column dimension of the matrix that we will compute the kernel of
    long cdim1;

    // above some threshold (or if just one mult_cdeg left),
    // just run the usual algo splitting column dimension in two equal parts
    if (index >= std::min<long>(threshold,mult_cdeg.size()-1))
    {
#ifdef PROFILE_DEG_UPMAT
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
#ifdef PROFILE_DEG_UPMAT
    t = GetWallTime();
#endif // PROFILE_DEG_UPMAT
    pmbasis(appbas, pmat_l, order, shifted_rdeg);
#ifdef PROFILE_DEG_UPMAT
    t = GetWallTime()-t;
    std::cout << prefix << "\tpmbasis dims x order = " << dim << " x " << cdim1 << " x " << order << " || time " << t << std::endl;
#endif // PROFILE_DEG_UPMAT

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
#ifdef PROFILE_DEG_UPMAT
    t = GetWallTime();
#endif // PROFILE_DEG_UPMAT
    multiply(pmatt, kerbas, pmat_r);
#ifdef PROFILE_DEG_UPMAT
    t = GetWallTime()-t;
    std::cout << prefix << "\tmultiply degrees " << deg(kerbas) << "," << deg(pmat_r) << " || time " << t << std::endl;
    const long actual_deg_ker = deg(kerbas);
    std::cout << prefix << "\tker deg, actual deg: " << deg_ker << "\t" << actual_deg_ker << std::endl;
#endif // PROFILE_DEG_UPMAT

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
#ifdef PROFILE_DEG_UPKER
    double t;
#endif // PROFILE_DEG_UPKER

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

    // retrieve the row degree and the determinant degree
    const long degdet = std::accumulate(diag_deg.begin(), diag_deg.end(), 0);

    // if degdet is not target_degree, then something went wrong in earlier steps
    // if they are equal, then we can proceed with the recursive calls
    if (degdet != target_degdet)
        return false;

    // if small dim, just run expansion by minors
    if (cdim<=4)
    {
#ifdef PROFILE_DEG_UPKER
        t=GetWallTime();
#endif // PROFILE_DEG_UPKER
        if (not kerbas_empty)
        {
            Mat<zz_pX> pmatt;
            multiply(pmatt, kerbas, pmat);
            determinant_expansion_by_minors(det, pmatt);
        }
        else
            determinant_expansion_by_minors(det, pmat);
#ifdef PROFILE_DEG_UPKER
        std::cout << prefix << "\tbase case --> " << GetWallTime()-t << std::endl;
#endif // PROFILE_DEG_UPKER
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
#ifdef PROFILE_DEG_UPKER
        std::cout << "\t-->Entering halving stage" << std::endl;
        prefix = '\t';
        t = GetWallTime();
#endif // PROFILE_DEG_UPKER
        multiply(pmatt, kerbas, pmat);
#ifdef PROFILE_DEG_UPKER
        t = GetWallTime()-t;
        std::cout << prefix << "\tmultiply degrees " << deg(kerbas) << "," << deg(pmat) << " || time " << t << std::endl;
#endif // PROFILE_DEG_UPKER
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
#ifdef PROFILE_DEG_UPKER
        std::cout << prefix << "\tkerbas*pmat_l multiply dimensions x degrees " << kerbas.NumRows() << "x" << rdim << "x" << cdim1 << "," << deg(kerbas) << "x" << deg(pmat_l) << " || time ";
        t = GetWallTime();
#endif // PROFILE_DEG_UPKER
        multiply(pmat_l, kerbas, pmat_l);
#ifdef PROFILE_DEG_UPKER
        t = GetWallTime()-t;
        std::cout << t << std::endl;
#endif // PROFILE_DEG_UPKER
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
#ifdef PROFILE_DEG_UPKER
    t = GetWallTime();
#endif // PROFILE_DEG_UPKER
    pmbasis(appbas, pmat_l, order, shifted_rdeg);
#ifdef PROFILE_DEG_UPKER
    t = GetWallTime()-t;
    std::cout << prefix << "\tpmbasis dims x order = " << app_rdim << " x " << cdim1 << " x " << order << " || time " << t << std::endl;
#endif // PROFILE_DEG_UPKER

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

    // if halving: compute kerbas2*pmat_r
    if (halving)
    {
#ifdef PROFILE_DEG_UPKER
        std::cout << prefix << "\tmultiply degrees " << deg(kerbas2) << "," << deg(pmat_r) << " || time ";
        t = GetWallTime();
#endif // PROFILE_DEG_UPKER
        multiply(pmat_r, kerbas2, pmat_r);
#ifdef PROFILE_DEG_UPKER
        t = GetWallTime()-t;
        std::cout << t << std::endl;
        const long actual_deg_ker = deg(kerbas2);
        std::cout << prefix << "\tker deg, actual deg: " << deg_ker << "\t" << actual_deg_ker << std::endl;
#endif // PROFILE_DEG_UPKER
    }
    // else if non-initial non-halving step: update kerbas by multiplication
    else if (kerbas.NumRows() > 0)
    {
#ifdef PROFILE_DEG_UPKER
        std::cout << prefix << "\tkerbas multiply degrees " << deg(kerbas2) << "x" << deg(kerbas) << " || time ";
        t = GetWallTime();
#endif // PROFILE_DEG_UPKER
        multiply(kerbas, kerbas2, kerbas);
#ifdef PROFILE_DEG_UPKER
        t = GetWallTime()-t;
        std::cout << t << std::endl;
        const long actual_deg_ker = deg(kerbas2);
        std::cout << prefix << "\tnew ker deg, actual deg: " << deg_ker << "\t" << actual_deg_ker << std::endl;
#endif // PROFILE_DEG_UPKER
    }
    // else, initial non-halving step: kerbas is empty, meaning identity
    // --> make it equal to kerbas2
    else
        kerbas.swap(kerbas2);

    // call recursively
    return determinant_shifted_form_degaware_updateone(det,kerbas,pmat_r,mult_cdeg,diag_deg,shifted_rdeg,index+1,threshold,target_degdet,prefix);
}

// kernel basis, ZLS style, of a full rank matrix of the form
//     [  x^d Ident ]
//     [  --------- ]  +   R
//     [      0     ]
// where cdeg(R) < d
// i.e. stacking a square row-wise -d-Popov form of -d-row degree 0, and a matrix reduced by it
// i.e. a column-wise 0-Popov form with column degree d and leading matrix the identity
// TODO verify mbasis is good when degrees is far from (or half of) order...
//      (since this is the main computational time spent here, at least in low degrees)
// TODO write formally that if matrix m x n with m > n has this shape, and d = d_1...d_n,
// then with the shift s=(d_1..d_n d_n..d_n) of length m,
// an s-weak Popov kernel basis K has all pivots in the right part
// and therefore |rdeg_s(K)| <= |s| = |d| + (m-n) d_n becomes |rdeg(K)| <= |d|
// --> generically deg(K) <= |d| / (m-n)   NO !! only the pivot part !!
// =====> nice observation but not so nice practice! the non-pivot part
// has higher degrees if d1<=d2<=...<=dn
// TODO more generally (outside of the targetted application context) improve
// ZLS kernel basis algorithm using this remark? note that we only need to
// *know* there exists a form of the input F with such a shape (or up to row
// permutation?) to apply this, indeed the left kernel of F is the same as that
// of F*U for nonsingular U.  Yet the shift may not be efficient if the actual
// degrees of F are far from those in its form having the above shape. For
// example, assume F already has this shape: can we efficiently apply a ZLS
// strategy to F directly? or do we need to apply some transformation to the
// columns of F to make sure its top rows are in fact row-wise 0-weak Popov
// (i.e. ensuring rdeg(top rows) = d)?
// TODO warning, does not work if numcols < 2numrows; and if numcols >= 2numrows
// the correctness may be only generic (to check)
void kernel_basis_zls_degaware_via_approximation(
                                                 Mat<zz_pX> & kerbas,
                                                 Mat<zz_pX> & pmat,
                                                 VecLong & degrees,
                                                 VecLong & shifted_rdeg,
                                                 VecLong & split_sizes
                                                )
{
    // kerbas[out]: the output kernel basis
    //     --> should be empty (dimensions 0) on initial call
    // pmat[in,out]: the input matrix, having the shape described above
    //     --> pmat is modified, used as a temporary
    // degrees[in,out]: the degrees d_1...d_m
    // shifted_rdeg: row degree shifted by some "s"
    //     --> initially, row degree 0, and all along, s = -initial column degree
    //     --> the matrix pmat remains s-owP all along, and shifted_rdeg gives its s-row deg
    //     --> 0-row-degree of kerbas is shifted_rdeg
    // split_sizes[in]: the indices at which we will split the matrix to
    //      perform the algorithm's recursion
    // index[in]: the index giving the split we are currently working on
    // threshold: after some threshold, just apply usual ZLS

#ifdef DEBUGGING_NOW
    std::cout << "degrees in kernel input matrix" << std::endl;
    std::cout << degree_matrix(pmat) << std::endl;
    std::cout << "degrees" << std::endl;
    std::cout << degrees << std::endl;
    std::cout << "shifted rdeg" << std::endl;
    std::cout << shifted_rdeg << std::endl;
#endif // DEBUGGING_NOW

    // row and column dimensions of pmat
    long rdim; long cdim;
    // size of block of columns we will handle in one given iteration
    long splitdim;
    // size of sought kernel
    long kdim;
    long degdet, order;
    // copy of shifted rdeg, for update purposes
    VecLong shift;
    // approximant basis for computing kernel
    Mat<zz_pX> appbas;
    // matrix to store local kernel
    Mat<zz_pX> kerbas2;
    // matrix to store copy of pmat for middle product with submatrix
    Mat<zz_pX> copy_pmat;
    // order accumulating the different approximation orders used
    long acc_order = 0;

    for (size_t index=0; index < split_sizes.size(); ++index)
    {
        rdim = pmat.NumRows();
        cdim = pmat.NumCols();
        splitdim = split_sizes[index];
        kdim = rdim-splitdim; // expected size of kerbas

        // find order for approximation
        degdet = std::accumulate(degrees.begin(), degrees.begin()+splitdim, 0) - splitdim*acc_order;
        long deg_ker = ceil(degdet / (double)kdim);
        order = *std::max_element(degrees.begin(), degrees.begin()+splitdim) - acc_order + deg_ker + 1;
        // TODO it seems sometimes we may take smaller order (e.g. -1 at first step)... reanalyze carefully...
        //order -= 1;
        acc_order += order;

        // save shift to be able to update diag_deg
        shift = shifted_rdeg;
#ifdef PROFILE_KER_AVG
        std::cout << index << " -> partial kernel: input matrix rdim, cdim, splitsize\n\t= (" << rdim << " x " << cdim << ") , " << splitdim  << std::endl;
        //std::cout << index << " -> column degree of input matrix" << std::endl << col_degree(pmat) << std::endl;
        double t;
        t = GetWallTime();
#endif // PROFILE_KER_AVG
        pmbasis(appbas, pmat, order, shifted_rdeg);
#ifdef PROFILE_KER_AVG
        t = GetWallTime()-t;
        std::cout << index << " -> expected kernel degree:" << deg_ker << std::endl;
        std::cout << index << " -> chosen appbas order:" << order << std::endl;
        VecLong rdeg_app = row_degree(appbas);
        std::cout << index << " -> found 'kernel' degree " << *max_element(rdeg_app.begin()+splitdim,rdeg_app.end()) << std::endl;
        std::cout << index << " -> pmbasis dims x order\n\t= (" << pmat.NumRows() << " x " << pmat.NumCols() << ") , " << order << " || time " << t << std::endl;
#endif // PROFILE_KER_AVG
        //std::cout << "\tchecking kernel does cancel..." << std::endl;
        //Mat<zz_pX> prod = appbas * pmat;
        //for (long i = splitdim; i < rdim; ++i)
        //    for (long j = 0; j < splitdim; ++j)
        //        if (not IsZero(prod[i][j]))
        //            std::cout << "KERNEL IS NOT KERNEL !!" << std::endl;
#ifdef DEBUGGING_NOW
        std::cout << "\trow degree of approx = " << std::endl << row_degree(appbas) << std::endl;
        std::cout << "\tdegree matrix of approx" << std::endl;
        std::cout << degree_matrix(appbas) << std::endl;
        std::cout << "\tdegree matrix of product" << std::endl;
        std::cout << degree_matrix(appbas*pmat) << std::endl;
#endif // DEBUGGING_NOW

        // update shifted_rdeg
        // shifted_rdeg-row degree of kerbas
        // = last entries of shift
        // = s-row degree of product kerbas*pmat_r (see description of function, for "s")
        std::vector<long>(shifted_rdeg.begin()+splitdim, shifted_rdeg.end()).swap(shifted_rdeg);

        // update diag_deg (keeping only the end of it)
        std::vector<long>(degrees.begin()+splitdim, degrees.end()).swap(degrees);
        for (long i = 0; i < kdim; ++i)
            degrees[i] += shifted_rdeg[i] - shift[splitdim+i];

        // extract the partial-kernel part of approx basis
        kerbas2.SetDims(kdim,rdim);
        for (long i = 0; i < kdim; ++i)
            for (long j = 0; j < rdim; ++j)
                kerbas2[i][j] = appbas[i+splitdim][j];

#ifdef DEBUGGING_NOW
        Mat<zz_pX> prod = kerbas2 * pmat;
        //std::cout << "degrees in kernel" << std::endl;
        //std::cout << degree_matrix(kerbas2) << std::endl;
        //std::cout << "degrees in prod" << std::endl;
        //std::cout << degree_matrix(prod) << std::endl;
#endif // DEBUGGING_NOW

        // update pmat by middle_product with appbas, only if necessary (i.e. not in last iteration)
        if (index<split_sizes.size()-1)
        {
            // -> ignoring top rows (i.e. top rows of appbas, i.e. only consider kerbas2*pmat)
            // -> ignoring left columns (kernel already computed)
#ifdef PROFILE_KER_AVG
            t = GetWallTime();
#endif
            copy_pmat.SetDims(rdim, cdim-splitdim);
            for (long i = 0; i < rdim; ++i)
                for (long j = 0; j < cdim-splitdim; ++j)
                    copy_pmat[i][j].swap(pmat[i][j+splitdim]);
            middle_product(pmat,kerbas2,copy_pmat,order,deg(copy_pmat)+deg(appbas)-order);
#ifdef PROFILE_KER_AVG
            t = GetWallTime()-t;
            std::cout << index << " -> updating middle product dims x degree x order || time\n\t= (" << kerbas2.NumRows() << " x " << kerbas2.NumCols() << ") * (" << copy_pmat.NumRows() << " x " << copy_pmat.NumCols() << ") , " << order << " , " << deg(copy_pmat)+deg(appbas)-order << " || " << t << std::endl;
#endif
        }

#ifdef PROFILE_KER_AVG
        std::cout << index << " -> updating multiply dims x degrees || time\n\t= (" << kerbas2.NumRows() << " x " << kerbas2.NumCols() << ") x (" << kerbas.NumRows() << " x " << kerbas.NumCols() << ") x " << deg(kerbas2) << " x " << deg(kerbas);
        t = GetWallTime();
#endif
        // update kerbas by multiplication
        if (index == 0) // initial call, kerbas was empty, corresponding to identity
            kerbas.swap(kerbas2);
        else
            multiply(kerbas,kerbas2,kerbas);
#ifdef PROFILE_KER_AVG
        t = GetWallTime()-t;
        std::cout << " || " << t << std::endl;
#endif
    }
}

bool determinant_shifted_form_zls_1(
                                    zz_pX & det,
                                    const Mat<zz_pX> & pmat,
                                    VecLong & diag_deg,
                                    VecLong & shifted_rdeg,
                                    VecLong & split_sizes,
                                    long target_degdet
                                   )
{
#ifdef PROFILE_ZLS_AVG
    double t;
#endif // PROFILE_ZLS_AVG

    // det: output determinant
    // pmat: input square matrix, assumed in "shiftedform" (X^{s} Id + R, cdeg(R) < s)
    // --> degree of determinant is sum(cdeg) = sum of diagonal degrees
    // mult_cdeg: list of lengths, sum should be pmat row (and column) dimension
    // diag_deg: list of diagonal degrees of pmat (= its s-pivot-degree)
    // shifted_rdeg: row degree shifted by some "s"
    //     (initially, row degree 0, and all along, s = -initial column degree)
    // --> the matrix pmat remains s-owP all along, and shifted_rdeg gives its s-row deg

    const long dim = pmat.NumRows();
    const long cdim = pmat.NumCols();

#ifdef DEBUGGING_NOW
    //std::cout << "shifted_rdeg: ";
    //std::copy(shifted_rdeg.begin(), shifted_rdeg.end(), std::ostream_iterator<long>(std::cout, " "));
    //std::cout << std::endl;
    //std::cout << "diag_deg: ";
    //std::copy(diag_deg.begin(), diag_deg.end(), std::ostream_iterator<long>(std::cout, " "));
    //std::cout << std::endl;
    //std::cout << "matrix degrees:" << std::endl;
    //std::cout << degree_matrix(pmat) << std::endl;
#endif // DEBUGGING_NOW

    // retrieve the row degree and the determinant degree
    const long degdet = std::accumulate(diag_deg.begin(), diag_deg.end(), 0);

#ifdef DEBUGGING_NOW
    std::cout << "verification:" << (degdet==target_degdet) << std::endl;
#endif // DEBUGGING_NOW

    // if degdet is not target_degree, then something went wrong in earlier steps
    // if they are equal, then we can proceed with the recursive calls
    if (degdet != target_degdet)
        return false;

    // if small dim, just run expansion by minors
    if (dim<=4)
    {
#ifdef PROFILE_ZLS_AVG
        t=GetWallTime();
#endif // PROFILE_ZLS_AVG
        determinant_expansion_by_minors(det, pmat);
#ifdef PROFILE_ZLS_AVG
        std::cout << "\tbase case --> " << GetWallTime()-t << std::endl;
#endif // PROFILE_ZLS_AVG
        return true;
    }

    // column dimension of the matrix that we will compute the kernel of
    // and number of splits this corresponds to
    long cdim1=0;
    long nb_splits=0;
    while (2*cdim1 < cdim)
    {
        cdim1 += split_sizes[nb_splits];
        ++nb_splits;
    }
    // if we exceeded cdim/2 (likely), set it to cdim/2
    // and modify split_sizes accordingly
    if (2*cdim1 > cdim)
    {
        long excess = cdim1 - (cdim>>1);
        split_sizes.insert(split_sizes.begin()+nb_splits, excess);
        split_sizes[nb_splits-1] -= excess;
        cdim1 = cdim>>1;
    }

    const long cdim2 = cdim-cdim1;

#ifdef DEBUGGING_NOW
    std::cout << "cdim1 , cdim2 , cdim" << std::endl;
    std::cout << cdim1 << "," << cdim2 << "," << cdim << std::endl;
#endif // DEBUGGING_NOW

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

    //VecLong copy_splits(split_sizes.begin(), split_sizes.begin()+nb_splits);
#ifdef PROFILE_ZLS_AVG
    t = GetWallTime();
#endif // PROFILE_ZLS_AVG
    Mat<zz_pX> kerbas;
    VecLong split_sz(split_sizes.begin(),split_sizes.begin()+nb_splits);
    kernel_basis_zls_degaware_via_approximation(kerbas, pmat_l, diag_deg, shifted_rdeg, split_sz);
#ifdef PROFILE_ZLS_AVG
    t = GetWallTime()-t;
    //std::cout << "\tkernel_basis dims x order = " << dim << " x " << cdim1 << " x " << order << " || time " << t << std::endl;
#endif // PROFILE_ZLS_AVG

    // update split_sizes (keeping only the end of it)
    std::vector<long>(split_sizes.begin()+nb_splits, split_sizes.end()).swap(split_sizes);

    //if (dim < 60)
    //{
    //    std::cout << degree_matrix(pmat_l) << std::endl;
    //    std::cout << degree_matrix(kerbas) << std::endl;
    //}

    // then compute the product
    Mat<zz_pX> pmatt;
#ifdef PROFILE_ZLS_AVG
    t = GetWallTime();
#endif // PROFILE_ZLS_AVG
    multiply(pmatt, kerbas, pmat_r);
#ifdef PROFILE_ZLS_AVG
    t = GetWallTime()-t;
    //std::cout << "\tmultiply degrees " << deg(kerbas) << "," << deg(pmat_r) << " || time " << t << std::endl;
#endif // PROFILE_ZLS_AVG
    //std::cout << degree_matrix(pmatt) << std::endl;
    // new degrees
    //diag_deg.resize(cdim2);
    //for (long j = 0; j < cdim2; ++j)
    //    diag_deg[j] = deg(pmat[j][j]);

#ifdef DEBUGGING_NOW
    std::cout << "degrees in updated pmat before recursive call" << std::endl;
    std::cout << degree_matrix(pmatt) << std::endl;
    std::cout << "diag_deg" << std::endl;
    std::cout << diag_deg << std::endl;
    std::cout << "shifted rdeg" << std::endl;
    std::cout << shifted_rdeg << std::endl;
    std::cout << "split_sizes" << std::endl;
    std::cout << split_sizes << std::endl;
#endif // DEBUGGING_NOW

    //return determinant_shifted_form_zls_1(det,pmatt,diag_deg,shifted_rdeg,split_sizes,index+1,target_degdet);
    // TODO call recursively
    return determinant_shifted_form_zls_1(det,pmatt,diag_deg,shifted_rdeg,split_sizes,target_degdet);
    //return determinant_generic_knowing_degree(det,pmatt,target_degdet);
    //kerbas.kill();
    //return determinant_shifted_form_degaware_updateall(det,pmatt,split_sizes,diag_deg,shifted_rdeg,0,3,target_degdet);
}

// ONE STEP degree 1
// kernel basis of a degree-1  (ell+2n) x n matrix of the form
//        [      F_0    ]
//        [  ---------- ]
//  F =   [      F_1    ]
//        [  ---------- ]
//        [  xI_n + F_2 ]
// where F_0 is   n x n, constant
//       F_1 is ell x n, constant
//       F_2 is   n x n, constant
// Assumption: ell >= 0, i.e. nrows >= 2 ncols
// Assumption: F_0 is invertible (not checked)
// Consequence 1: left nullspace of F mod x has the form
//   [  K_0  |  I_ell  |  0  ]
//   [  ----   -------   --- ]
//   [  K_1  |    0    | I_n ]
// where
//   [ K_0 ]       [ F_1 ]
//   [ --- ]  =    [ --- ] * iF_0
//   [ K_1 ]       [ F_2 ]
// where iF_0 = inverse(-F_0)
// Consequence 2: one reduced left kernel basis of F is
//   [  K_1 - x iF_0  |    0    |  I_n ]
//   [  -------------   -------   ---- ]
//   [       K_0      |  I_ell  |   0  ]
// (note that it is also in HNF up to row permutation)
// and the Popov left kernel basis of F is
//   [  -F_0 K_1 + x I_n  |    0    | -F_0 ]
//   [   ---------------    -------   ---- ]
//   [         K_0        |  I_ell  |   0  ]
// Input : (ell+2n) x n constant matrix
//        [  F_0  ]
//        [  ---- ]
//  F =   [  F_1  ]
//        [  ---- ]
//        [  F_2  ]
// Output : 
//   - kertop: n x n constant matrix kertop (above: -F_0 K_1)
//   - kerbot: ell x n constant matrix kerbot (above: K_0)
//   - cmat: replaced by the negation of its n top rows (above: -F_0)
void kernel_step1(Mat<zz_p> & kertop, Mat<zz_p> & kerbot, Mat<zz_p> & cmat)
{
#ifdef PROFILE_KERNEL_STEP1
    double t_total,t_gauss,t_mul;
    t_total = GetWallTime();
#endif
    // get dimensions, check numrows >= 2 numcols
    const long n = cmat.NumCols();
    const long ell = cmat.NumRows() - 2*n;
    if (ell < 0)
    {
        std::cout << "~~ERROR~~ numrows >= 2numcols required in kernel_step1; see docstring" << std::endl;
        return;
    }

#ifdef PROFILE_KERNEL_STEP1
    t_gauss = GetWallTime();
#endif
    // perform Gaussian elimination to retrieve kernel basis
    Mat<zz_p> kermat;
    kernel(kermat, cmat);
#ifdef PROFILE_KERNEL_STEP1
    t_gauss = GetWallTime() - t_gauss;
#endif

    // copy the left part
    kertop.SetDims(n,n);
    for (long i = 0; i < n; ++i)
        VectorCopy(kertop[i], kermat[ell+i], n);
    kerbot.SetDims(ell,n);
    for (long i = 0; i < ell; ++i)
        VectorCopy(kerbot[i], kermat[i], n);
    cmat.SetDims(n,n);
#ifdef PROFILE_KERNEL_STEP1
    t_mul = GetWallTime();
#endif
    NTL::negate(cmat,cmat);
    mul(kertop, cmat, kertop);
#ifdef PROFILE_KERNEL_STEP1
    t_mul = GetWallTime() - t_mul;
    t_total = GetWallTime() - t_total;
    std::cout << "\tKernel step1 profile:" << std::endl;
    std::cout << "\ttotal time:\t" << t_total << std::endl;
    std::cout << "\tGauss time:\t" << t_gauss << std::endl;
    std::cout << "\tMatMul time:\t" << t_mul << std::endl;
#endif
}

// same as above, but using direct computation of K via inv(.) instead of kernel(.)
void kernel_step1_direct(Mat<zz_p> & kertop, Mat<zz_p> & kerbot, Mat<zz_p> & cmat)
{
#ifdef PROFILE_KERNEL_STEP1
    double t_total,t_kernel,t_mul;
    t_total = GetWallTime();
#endif
    // get dimensions, check numrows >= 2 numcols
    const long n = cmat.NumCols();
    const long ell = cmat.NumRows() - 2*n;
    if (ell < 0)
    {
        std::cout << "~~ERROR~~ numrows >= 2numcols required in kernel_step1; see docstring" << std::endl;
        return;
    }

#ifdef PROFILE_KERNEL_STEP1
    t_kernel = GetWallTime();
#endif
    // retrieve kernel basis by direct computation
    kertop.SetDims(n,n);
    for (long i = 0; i < n; ++i)
        kertop[i].swap(cmat[ell+n+i]);
    kerbot.SetDims(ell,n);
    for (long i = 0; i < ell; ++i)
        kerbot[i].swap(cmat[n+i]);

    Mat<zz_p> imat;
    cmat.SetDims(n,n);
    NTL::negate(cmat,cmat);
    inv(imat,cmat);
    mul(kertop, kertop, imat);
    mul(kerbot, kerbot, imat);
#ifdef PROFILE_KERNEL_STEP1
    t_kernel = GetWallTime() - t_kernel;
    t_mul = GetWallTime();
#endif
    mul(kertop, cmat, kertop);
#ifdef PROFILE_KERNEL_STEP1
    t_mul = GetWallTime() - t_mul;
    t_total = GetWallTime() - t_total;
    std::cout << "\tKernel step1 profile:" << std::endl;
    std::cout << "\ttotal time:\t" << t_total << std::endl;
    std::cout << "\tKernel time:\t" << t_kernel << std::endl;
    std::cout << "\tMatMul time:\t" << t_mul << std::endl;
#endif
}

// ONE STEP in degree 2
// ---> redundant with below, degree d
// kernel basis of a degree 2  (ell+3n) x n matrix of the form
//        [       F_0      ]
//        [   -----------  ]
//  F =   [       F_1      ]
//        [   -----------  ]
//        [  x^2 I_n + F_2 ]
// where F_0 is  2n x n, degree <= 1
//       F_1 is ell x n, degree <= 1
//       F_2 is   n x n, degree <= 1
// Assumption: ell >= 0, i.e. nrows >= 3 ncols
// Assumption: [F_00 | F_01] is invertible (not checked)
//       where F_00 is the constant coefficient of F_0
//         and F_01 is the coefficient of degree 1 of F_0
// Consequence 1: left nullspace of
//   [ F_00 | F_01 ]
//   [ -----  ---- ]
//   [ F_10 | F_11 ]
//   [ -----  ---- ]
//   [ F_20 | F_21 ]
// has the form
//   [  K_0  |  I_ell  |  0  ]
//   [  ----   -------   --- ]
//   [  K_1  |    0    | I_n ]
// where
//   [ K_0 ]       [ F_10 | F_11 ]
//   [ --- ]  =    [ ----   ---- ] * iF_0
//   [ K_1 ]       [ F_20 | F_21 ]
// where iF_0 = inverse([-F_00 | -F_01])
// Note: all rows in the top part (involving K_0) are
//       obviously also in the left kernel of F; in fact
//       they form the constant part of Ker(F)
// Note: we have the property that
//                       [  I_n  ]
//      iF_0 * F_0  =  - [ ----- ]
//                       [ x I_n ]
//   which gives trivial syzygies E = [ x I_n | -I_n ]
// Consequence 2: one reduced left kernel basis of F is
//   [  -iF_0b + x iF_0t |    0    |   0  ]
//   [  ---------------    -------   ---- ]
//   [   K_1 + x iF_0b   |    0    |  I_n ]
//   [  --------------     -------   ---- ]
//   [        K_0        |  I_ell  |   0  ]
// where iF_0b and iF_0t are bottom and top of iF_0
//   -- note that E iF_0 = -iF_0b + x iF_0t
//   -- note that this matrix is also in HNF up to row permutation
//   -- note that top left 2n x 2n block has leading matrix
//      equal to iF_0, which is invertible
// and the Popov left kernel basis of F is
//   [  R + x I_2n  |     0   | -F_01 ]
//   [  -----------   -------   ----- ]
//   [     K_0      |  I_ell  |   0   ]
//   where
//                            [ iF_0b ]
//      R = [ -F_00 | -F_01 ] [ ----- ] 
//                            [  K_1  ]
// Input : (ell+3n) x 2n constant matrix
//   [ F_00 | F_01 ]
//   [ ----- ----- ]
//   [ F_10 | F_11 ]
//   [ ----- ----- ]
//   [ F_20 | F_21 ]
// Output : 
//   - kertop: 2n x 2n constant matrix (R above)
//   - kerbot: ell x n constant matrix kerbot (K_0 above)
//   - cmat: replaced by only the degree 1 part of its n top rows (-F_01 above)
void kernel_step2(Mat<zz_p> & kertop, Mat<zz_p> & kerbot, Mat<zz_p> & cmat)
{
#ifdef PROFILE_KERNEL_STEP2
    double t_total,t_kernel,t_mul;
    t_total = GetWallTime();
#endif
    // get dimensions, check numrows >= 2 numcols
    const long n = cmat.NumCols()>>1;
    const long ell = cmat.NumRows() - 3*n;
    if (ell < 0)
    {
        std::cout << "~~ERROR~~ numrows >= 3numcols required in kernel_step2; see docstring" << std::endl;
        return;
    }

#ifdef PROFILE_KERNEL_STEP2
    t_kernel = GetWallTime();
#endif
    // retrieve kernel basis by direct computation
    kertop.SetDims(n,2*n);
    for (long i = 0; i < n; ++i)
        kertop[i].swap(cmat[ell+2*n+i]);
    kerbot.SetDims(ell,2*n);
    for (long i = 0; i < ell; ++i)
        kerbot[i].swap(cmat[2*n+i]);

    cmat.SetDims(2*n,2*n);
    NTL::negate(cmat,cmat);
    Mat<zz_p> imat;
    inv(imat,cmat); // imat == iF_0
    mul(kertop, kertop, imat); // kertop = K_1
    mul(kerbot, kerbot, imat); // kerbot = K_0
#ifdef PROFILE_KERNEL_STEP2
    t_kernel = GetWallTime() - t_kernel;
    t_mul = GetWallTime();
#endif
    kertop.SetDims(2*n,2*n);
    for (long i = 0; i < n; ++i)
    {
        // row n+i <- row i of K1
        kertop[i].swap(kertop[n+i]);
        // row i <- - row i of iF_0b = - row n+i of iF_0
        kertop[i].swap(imat[n+i]);
        NTL::negate(kertop[i],kertop[i]);
    }
    mul(kertop, cmat, kertop);
#ifdef PROFILE_KERNEL_STEP2
    t_mul = GetWallTime() - t_mul;
#endif
    // remove left columns of cmat to keep only -F_01
    Mat<zz_p> tmp;
    tmp.SetDims(2*n,n);
    for (long i = 0; i < 2*n; ++i)
        for (long j = 0; j < n; ++j)
            tmp[i][j] = cmat[i][n+j];
    cmat.swap(tmp);
#ifdef PROFILE_KERNEL_STEP2
    t_total = GetWallTime() - t_total;
    std::cout << "\tKernel step2 profile:" << std::endl;
    std::cout << "\ttotal time:\t" << t_total << std::endl;
    std::cout << "\tKernel time:\t" << t_kernel << std::endl;
    std::cout << "\tMatMul time:\t" << t_mul << std::endl;
#endif
}

// Compute kernel basis of degree 1:
// kernel basis of a degree d  (ell+(d+1)n) x n matrix of the form
//        [       F_0      ]
//        [   -----------  ]
//  F =   [       F_1      ]
//        [   -----------  ]
//        [  x^d I_n + F_2 ]
// where F_0 is  dn x n, degree < d
//       F_1 is ell x n, degree < d
//       F_2 is   n x n, degree < d
// Assumption: ell >= 0, i.e. nrows >= (d+1) ncols
// Assumption: [F_00 | F_01 | .. | F_0l] is invertible (not checked)
//       where F_0k is the coefficient of degree k of F_0
//       where F_0l is the leading coefficient (coeff of degree d-1) of F_0
// Consequence 1: left nullspace of
//   [ F_00 | F_01 | ... | F_0l ]
//   [ ----   ----   ---   ---- ]
//   [ F_10 | F_11 | ... | F_1l ]
//   [ ----   ----   ---   ---- ]
//   [ F_20 | F_21 | ... | F_2l ]
// has the form
//   [  K_0  |  I_ell  |  0  ]
//   [  ----   -------   --- ]
//   [  K_1  |    0    | I_n ]
// where
//   [ K_0 ]      [  F_10 |  F_11 |  ...  |  F_1l ] 
//   [ --- ]  =   [ -----   -----   -----   ----- ] * iF_0
//   [ K_1 ]      [  F_20 |  F_21 |  ...  |  F_2l ] 
// where iF_0 = inverse([-F_00 | -F_01 | ... | -F_0l])
// Note: all rows in the top part (involving K_0) are
//       obviously also in the left kernel of F; in fact
//       they form the constant part of Ker(F)
// Note: we have the property that
//                     [     I_n     ]
//                     [  ---------  ]
//                     [    x I_n    ]
//   iF_0 * F_0 =  -   [  ---------  ]
//                     [      :      ]
//                     [  ---------  ]
//                     [ x^(d-1) I_n ]
//   which gives (d-1)*n trivial syzygies, via the (d-1)n x dn matrix
//          [ x I_n |  -I_n |      |     |       |      ]
//          [       | x I_n | -I_n |     |       |      ]
//    E =   [       |       |  ... | ... |       |      ]
//          [       |       |      | ... |  ...  |      ]
//          [       |       |      |     | x I_n | -I_n ]
// Consequence 2: one reduced left kernel basis of F is
//   [      E iF_0       |    0    |   0  ]
//   [  ---------------    -------   ---- ]
//   [  K_1 + x iF_0bb   |    0    |  I_n ]
//   [  --------------     -------   ---- ]
//   [        K_0        |  I_ell  |   0  ]
// where iF_0bb is the bottom n x dn part of iF_0
// -- note that this matrix is also in HNF up to row permutation
// -- note that top left dn x dn block has leading matrix
//                 equal to iF_0, which is invertible
// and the Popov left kernel basis of F is
//   [  R + x I_dn  |     0   | -F_0l ]
//   [  -----------   -------   ----- ]
//   [     K_0      |  I_ell  |   0   ]
//   where
//                                         [ E(0) iF_0b ]
//      R = [ -F_00 | -F_01 | .. | -F_0l ] [  --------  ]
//                                         [     K_1    ]
//   where
//       iF_0b = (d-1)n x dn bottom part of iF_0 (everything except n first rows)
// Input : degree d and  (ell+(d+1)n) x dn constant matrix
//           [ F_00 | F_01 | ... | F_0l ]
//           [ ----   ----   ---   ---- ]
//   cmat =  [ F_10 | F_11 | ... | F_1l ]       (where l = d-1)
//           [ ----   ----   ---   ---- ]
//           [ F_20 | F_21 | ... | F_2l ]
// Output/Side-effect :
//   - kertop: dn x dn constant matrix (R above)
//   - kerbot: ell x dn constant matrix kerbot (K_0 above)
//   - cmat: replaced by only the degree d-1 part of its n top rows (-F_0l above)
void kernel_degree1(Mat<zz_p> & kertop, Mat<zz_p> & kerbot, Mat<zz_p> & cmat, long d)
{
#ifdef PROFILE_KERNEL_DEGREE1
    double t_total,t_kernel,t_mul;
    t_total = GetWallTime();
#endif
    // get dimensions, check numrows >= (degree+1) * numcols
    const long n = cmat.NumCols()/d;
    const long ell = cmat.NumRows() - (d+1)*n;
    if (ell < 0)
    {
        std::cout << "~~ERROR~~ numrows >= (degree+1)*numcols required in kernel_degree1; see docstring" << std::endl;
        return;
    }

#ifdef PROFILE_KERNEL_DEGREE1
    t_kernel = GetWallTime();
#endif
    // retrieve kernel basis by direct computation
    kertop.SetDims(n,d*n);
    for (long i = 0; i < n; ++i)
        kertop[i].swap(cmat[ell+d*n+i]);
    kerbot.SetDims(ell,d*n);
    for (long i = 0; i < ell; ++i)
        kerbot[i].swap(cmat[d*n+i]);

    cmat.SetDims(d*n,d*n); // forget bottom rows, already copied
    NTL::negate(cmat,cmat);
    Mat<zz_p> imat;
    inv(imat,cmat); // imat == iF_0
    mul(kertop, kertop, imat); // kertop = K_1
    mul(kerbot, kerbot, imat); // kerbot = K_0
#ifdef PROFILE_KERNEL_DEGREE1
    t_kernel = GetWallTime() - t_kernel;
#endif
    kertop.SetDims(d*n,d*n); // augment with iF_0b
    for (long i = 0; i < n; ++i)
    {
        // row (d-1)*n+i <- row i of K1
        kertop[i].swap(kertop[(d-1)*n+i]); 
    }
    for (long i = 0; i < (d-1)*n; ++i)
    {
        // row i <- -row i of iF_0b = -row n+i of iF_0
        kertop[i].swap(imat[n+i]);
        NTL::negate(kertop[i], kertop[i]);
    }
#ifdef PROFILE_KERNEL_DEGREE1
    t_mul = GetWallTime();
#endif
    mul(kertop, cmat, kertop);
#ifdef PROFILE_KERNEL_DEGREE1
    t_mul = GetWallTime() - t_mul;
#endif
    // remove left columns of cmat to keep only -F_01
    Mat<zz_p> tmp;
    tmp.SetDims(d*n,n);
    for (long i = 0; i < d*n; ++i)
        for (long j = 0; j < n; ++j)
            tmp[i][j] = cmat[i][(d-1)*n+j];
    cmat.swap(tmp);
#ifdef PROFILE_KERNEL_DEGREE1
    t_total = GetWallTime() - t_total;
    std::cout << "\tKernel-degree1 profile:" << std::endl;
    std::cout << "\ttotal time:\t" << t_total << std::endl;
    std::cout << "\tKernel time:\t" << t_kernel << std::endl;
    std::cout << "\tMatMul time:\t" << t_mul << std::endl;
#endif
}

// Compute kernel basis of degree 1:
// Goal:
// """""
// compute kernel basis of a (D+ell+n) x n matrix F of column degree d=cdeg,
// where D = sum of column degrees = d_1 + ... + d_n,
// and F has the form
//        [      F0      ]
//        [  ----------  ]
//  F =   [      F1      ]
//        [  ----------  ]
//        [ x^d I_n + F2 ]
// where F0 is   D x n, cdeg(F0) < d
//       F1 is ell x n, cdeg(F1) < d
//       F2 is   n x n, cdeg(F2) < d.
//
// Notation:
// """""""""
//   * eF the (D+ell+n) x D matrix obtained by fully linearizing F
//         eF = [ eF_0 | eF_1 | ... | eF_{mu-1} ]
//         where mu = max(cdeg)
//         and eF_k has dimensions m x n_k
//         (i.e. n_k is the max integer such that d[n_k-1] > k)
//   * eF0 (resp. eF1, eF2) the D x D (resp ell x D, n x D) matrix
//       obtained by fully linearizing F0 (resp. F1, F2)
//   * eFbot = [[eF1],[eF2]] the (ell+n) x D matrix obtained by fully linearizing [[F1],[F2]]
//   * iF0 the inverse of -eF0  (note the minus sign)
//
// Assumption:
// """""""""""
//   * ell >= 0, i.e. nrows >= n+D
//   * eF0 is invertible (not checked)
//
// Consequence:
// """"""""""""
// left nullspace of eF has a basis of the form
//   [  K0  |  I_ell  |  0  ]
//   [ ----   -------   --- ]
//   [  K1  |    0    | I_n ]
// where
//   [ K0 ]    [ eF1 ]   
//   [ -- ]  = [ --- ] * iF0
//   [ K1 ]    [ eF2 ]   
//
// Note: all rows in the top part (involving K_0) are
//       obviously also in the left kernel of F; in fact
//       they form the constant part of Ker(F)
//
// Note: unlike in the above implementations, here we will
// not use K0 and K1 as such, but a more direct computation
// of the Popov kernel basis, see below
//
// Consequence:
// """"""""""""
// the Popov left kernel basis of F is
//   [  xI_D + K0 |    0    | -lF0  ]
//   [  ---------   -------   ----- ]
//   [     K_1    |  I_ell  |   0   ]
// where
//  [ K0 ]       [ lF0 eF2 - sF0 ]
//  [ -- ]   =   [ ------------- ] * iF0
//  [ K1 ]       [     -eF1      ]
// where
//   lF0 = column-leading-matrix(F0, [d[0]-1, .., d[n-1]-1])
//   sF0 = [ 0 | eF0_0 | ... | eF0_{mu-2} ]
//   mu = max(cdeg)
//   0 has size D x n
//   eF0_k is the D x n_{k+1} left submatrix of the coefficient of degree k of F0
//
// Input:
// """"""
//   - matp: polynomial matrix (D+ell+n) x n  (corresponding to F without x^d Id part)
//   - cdeg: column degrees of F (d above)
//   - if starting_index > 0, then it will be considered that only the columns/entries
//     of matp and cdeg on and after starting_index should be considered
//
// Output/Side-effect:
// """""""""""""""""""
//   - K0:  D x D constant matrix
//   - K1: ell x D constant matrix kerbot
//   - lF0: D x n constant matrix ( -lF0 above, note the minus sign)
//
// Required:
// """""""""
//   - eF0 invertible, as mentioned above
//   - cdeg is nonincreasing and cdeg > 0 (so that n = matp[0].NumCols())
//   - nrows >= n+D, where D == sum(cdeg)
void kernel_degree1_cdeg(Mat<zz_p> & K0, Mat<zz_p> & K1, Mat<zz_p> & lF0, const Vec<Mat<zz_p>> & matp, const VecLong & cdeg, long starting_index=0)
{
#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    double t_total,t_mulinv,t_main;
    std::cout << "\tKernel-degree1-cdeg starting, m | n | D | ell | mu = ";
    t_main=0;
    t_total = GetWallTime();
#endif
    // get dimensions, check numrows >= D+n
    const long D = std::accumulate(cdeg.begin()+starting_index, cdeg.end(), 0);
    const long n = matp[0].NumCols()-starting_index;
    const long ell = matp[0].NumRows() - D - n;
    const long mu = cdeg[starting_index];
#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    std::cout << matp[0].NumRows() << " | " << n << " | " << D << " | " << ell << " | " << mu << std::endl;
#endif
    if (ell < 0)
    {
        std::cout << "~~ERROR~~ m >= D+n required in kernel_degree1_cdeg; see docstring" << std::endl;
        return;
    }

    // build eF0 and sF0
    // store -eF1 in K1, eF2 in K0, lF0 in lF0
    Mat<zz_p> eF0, sF0;
    eF0.SetDims(D,D);
    sF0.SetDims(D,D);
    K0.SetDims(n,D);
    K1.SetDims(ell,D);
    lF0.SetDims(D,n);

    long d = 0;
    for (long k = 0; k < mu; ++k)
    {
        const long n_k = matp[k].NumCols()-starting_index;
        const long n_kpp = (k < mu-1) ? (matp[k+1].NumCols()-starting_index) : 0;

        // eF0 <- eF0 and sF0 <- sF0 and lF0 <- lF0
        for (long i = 0; i < D; ++i)
        {
            for (long j = 0; j < n_kpp; ++j)
            {
                eF0[i][d+j] = matp[k][i][starting_index+j];
                sF0[i][n_k+d+j] = matp[k][i][starting_index+j];
            }
            for (long j = n_kpp; j < n_k; ++j)
            {
                eF0[i][d+j] = matp[k][i][starting_index+j];
                lF0[i][j] = matp[k][i][starting_index+j];
            }
        }

        // K1 <- -eF1
        for (long i = 0; i < ell; ++i)
            for (long j = 0; j < n_k; ++j)
                NTL::negate(K1[i][d+j], matp[k][D+i][starting_index+j]);

        // K0 <- eF2
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n_k; ++j)
                K0[i][d+j] = matp[k][D+ell+i][starting_index+j];

        d += n_k;
    }

#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_mulinv = GetWallTime();
#endif
    // K0 = lF0 * eF2
    mul(K0, lF0, K0);
#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_mulinv = GetWallTime() - t_mulinv;
    t_main += t_mulinv;
    std::cout << "\t\tMatMul1 --> dimensions\t" << D << " x " << n << " x " << D << " ||  time " << t_mulinv << std::endl;
#endif

    // lF0 = -lF0
    NTL::negate(lF0, lF0);
    // K0 = lF0 * eF2 - sF0
    sub(K0, K0, sF0); 
    //       [ lF0 * eF2 - sF0 ]
    // K0 <- [ --------------- ]
    //       [      -eF1       ]
    K0.SetDims(D+ell, D);
    for (long i = 0; i < ell; ++i)
        K0[i+D].swap(K1[i]);

#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_mulinv = GetWallTime();
#endif
    // iF0 <- inverse of eF0
    Mat<zz_p> iF0;
    inv(iF0,eF0);
#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_mulinv = GetWallTime() - t_mulinv;
    t_main += t_mulinv;
    std::cout << "\t\tInversion --> dimensions\t" << D << " x " << D << " ||  time " << t_mulinv << std::endl;
#endif

#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_mulinv = GetWallTime();
#endif
    //       [ lF0 * eF2 - sF0 ]
    // K0 <- [ --------------- ]  iF0
    //       [      -eF1       ]
    mul(K0, K0, iF0);
#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_mulinv = GetWallTime() - t_mulinv;
    t_main += t_mulinv;
    std::cout << "\t\tMatMul2 --> dimensions\t" << D+ell << " x " << D << " x " << D << " ||  time " << t_mulinv << std::endl;
#endif

    // split back into K0 and K1
    for (long i = 0; i < ell; ++i)
        K0[i+D].swap(K1[i]);
    K0.SetDims(D,D);

#ifdef PROFILE_KERNEL_DEGREE1_CDEG
    t_total = GetWallTime() - t_total;
    std::cout << "\ttotal time " << t_total << ", non-dominant ";
    std::cout << std::setprecision(2) << 100*(t_total-t_main)/t_total << "%" << std::setprecision(8) << std::endl;
#endif
}

// Given matp represented by a matrix, which is gluing together the 
// coefficients matrices of sizes respectively sz[0], sz[1], ...
// (first entries are high degree terms)
// and given a constant matrix mat,
// compute fmatp which is
//     ( [diag([x]*k + [1]*(m-k) | 0]  +  mat) * matp
// represented the same way
// 
// Warning: 
// sizes is updated in place
// matp is modified
void fused_mul_add(Mat<zz_p> & fmatp, Mat<zz_p> & matp, long rdim, Mat<zz_p> & mat, long k)
{
    const long n = matp.NumCols();
    const long mm = mat.NumRows(); // new row dim
    const long nn = n+rdim; // new column dim

    double t = GetWallTime();
    mul(mat, mat, matp);
    cout << "MUL --> " << GetWallTime() - t << endl;

    fmatp.SetDims(mm,nn);
    for (long i = 0; i < k; ++i)
    {
        for (long j = 0; j < rdim; ++j)
            fmatp[i][j].LoopHole() = mat[i][j]._zz_p__rep;
        for (long j = rdim; j < n; ++j)
            add(fmatp[i][j], matp[i][j-rdim], mat[i][j]);
        for (long j = n; j < nn; ++j)
            fmatp[i][j].LoopHole() = matp[i][j-rdim]._zz_p__rep;
    }
    for (long i = k; i < mm; ++i)
        for (long j = 0; j < n; ++j)
            fmatp[i][j], mat[i][j]._zz_p__rep;
}

// Convert a matrix polynomial stored as below in determinant_shifted_form_smartkernel_updateall
// towards a polynomial matrix
// Actually this supports the rectangular m x n case with n < m; then the x id part is at the bottom
void conv_shifted_form(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp, const VecLong & cdeg)
{
    const long len = matp.length();
    const long n = cdeg.size();
    if (len == 0)
    {
        // convention: matp is the n x n identity matrix
        // (number of rows is a bit undefined...
        std::cout << "WARNING: conv_shifted_form with zero length matp" << std::endl;
        ident(pmat, n);
        return;
    }

    const long m = matp[0].NumRows();
    if (m<n)
    {
        std::cout << "~~ERROR~~ conv_shifted_form requires m >= n" << std::endl;
        return;
    }

    pmat.SetDims(m,n);

    for (long i = 0; i < m; ++i)
    {
        for (long j = 0; j < n; j+=CACHE_FRIENDLY_SIZE)
        {
            const long j_bound = std::min((long)CACHE_FRIENDLY_SIZE, n-j);
            // initialize vectors
            for (long jj = 0; jj < j_bound; ++jj)
                pmat[i][j+jj].SetLength(cdeg[j+jj]);

            // fill data (columns j to j+CACHE_FRIENDLY_SIZE-1
            // note: cdeg must be nonincreasing!
            for (long k = 0; k < cdeg[j]; k+=MATRIX_BLOCK_SIZE)
            {
                for (long jj = 0; jj < j_bound; ++jj)
                {
                    const long k_bound = std::min((long)MATRIX_BLOCK_SIZE, cdeg[j+jj]-k);
                    for (long kk = 0; kk < k_bound; ++kk)
                        pmat[i][j+jj][k+kk] = matp[k+kk][i][j+jj];
                }
            }

            // strip away leading zeroes
            for (long jj = 0; jj < j_bound; ++jj)
                pmat[i][j+jj].normalize();
        }
    }

    // add x^cdeg Id
    for (long j = m-n; j < n; ++j)
        SetCoeff(pmat[m-n+j][j],cdeg[j]);

    return;
}

// same as above but applying J on right and left
// TODO not optimal
void conv_shifted_form_mirrored(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp, const VecLong & cdeg)
{
    Mat<zz_pX> buf;
    conv_shifted_form(buf, matp, cdeg);
    const long m = buf.NumRows();
    const long n = buf.NumCols();
    pmat.SetDims(m,n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            pmat[m-1-i][n-1-j] = buf[i][j];
}

Mat<ZZ> degree_matrix(const Vec<Mat<zz_p>> & matp)
{
    const long m = matp[0].NumRows();
    const long n = matp[0].NumCols();
    Mat<ZZ> degmat;
    degmat.SetDims(m,n);
    for (long j = 0; j < n; ++j)
    {
        long k = 0;
        while (k < matp.length() && matp[k].NumCols() > j)
            ++k;
        for (long i = 0; i < m; ++i)
            degmat[i][j] = k-1;
        degmat[j][j] += 1;
    }
    return degmat;
}


// det: output determinant
// matp: input square m x m matrix
// cdeg: column degrees of the input matp (len(cdeg) = m), assumed nondecreasing
// threshold: when "??? TODO", switch to a classical determinant algorithm
//      where d is the smallest entry (i.e. last entry)
// target_degdet: degree of determinant of matp, used for certification
// REQUIREMENTS:
//   -- matp is in "shiftedform": matp = X^{cdeg} Id + R, where cdeg(R) < cdeg
//   -- cdeg is nonincreasing   (cdeg[j+1] <= cdeg[j])
//   -- storage: matp is a list of coefficient matrices; matp[i] has dimensions m x n_i
//      where n_i is the largest integer j such that cdeg[j] >= i
//      NOTE: the X^cdeg Id part is not stored
//   -- target_degdet: must be equal to deg(det(matp)), otherwise the
//   behaviour of certification is not guaranteed
// GUARANTEES IN OUTPUT:
//   -- det = det(matp)
//   -- if target_degdet in input, then the return value is true iff det = det(matp)
// PROPERTIES:
//   -- degree of determinant is sum(cdeg) = sum of diagonal degrees
bool determinant_shifted_form_smartkernel_updateall(zz_pX & det, Vec<Mat<zz_p>> & matp, VecLong & cdeg, long threshold, long target_degdet)
{
#ifdef DEBUGGING_NOW
    std::cout << target_degdet << std::endl;
    std::cout << cdeg << std::endl;
#endif // DEBUGGING_NOW
    
#ifdef PROFILE_SMART_UPMAT
    double t_total, t_kernel, t_update, t;
    t_update = 0;
    t_total = GetWallTime();
#endif // PROFILE_SMART_UPMAT
    if (matp.length() == 0) // identity matrix
    {
        det = 1;
#ifdef PROFILE_SMART_UPMAT
    std::cout << "\tbase case identity matrix" << std::endl;
#endif // PROFILE_SMART_UPMAT
        return (target_degdet == 0);
    }

    const long m = matp[0].NumRows();
#ifdef PROFILE_SMART_UPMAT
    std::cout << "Entering iteration with m , mindeg = " << m << " , " << cdeg[m-1] << std::endl;
#endif // PROFILE_SMART_UPMAT
    if (m != (long)cdeg.size())
    {
        std::cout << "~~ERROR~~ wrong length of input cdeg in determinant_shifted_form_smartkernel_updateall" << std::endl;
        return false;
    }

    // if sum(cdeg) is not target_degree, then something went wrong in earlier steps
    if (std::accumulate(cdeg.begin(), cdeg.end(), 0) != target_degdet)
    {
        std::cout << "~~ERROR~~  in determinant_shifted_form_smartkernel_updateall" << std::endl;
        std::cout << "         --> target_degdet is not equal to sum(cdeg)" << std::endl;
        std::cout << "         --> either cdeg is wrong" << std::endl;
        std::cout << "           (which is not supposed to happen if it was correct in input)" << std::endl;
        std::cout << "         --> or some kernel computation went wrong" << std::endl;
        std::cout << "           (which can happen due to genericity assumption not satisfied)" << std::endl;
        return false;
    }

    // if small dim, just run expansion by minors
    if (m<=4)
    {
        Mat<zz_pX> pmat;
        conv_shifted_form(pmat, matp, cdeg);
#ifdef PROFILE_SMART_UPMAT
        t=GetWallTime();
#endif // PROFILE_SMART_UPMAT
        determinant_expansion_by_minors(det, pmat);
#ifdef PROFILE_SMART_UPMAT
        std::cout << "\tbase case (expansion) --> " << GetWallTime()-t << std::endl;
        std::cout << "\ttotal of iteration    --> " << GetWallTime()-t_total << std::endl;
#endif // PROFILE_SMART_UPMAT
        return true;
    }

    const long d = cdeg[m-1];

    // above some threshold,
    // run the usual algo splitting column dimension in two equal parts
    if (d >= threshold || m <= d+1) // TODO
    {
#ifdef PROFILE_SMART_UPMAT
        std::cout << "\t-->Entering halving stage" << std::endl;
#endif // PROFILE_SMART_UPMAT
        Mat<zz_pX> pmat;
        conv_shifted_form_mirrored(pmat, matp, cdeg);

#ifdef PROFILE_SMART_UPMAT
        t=GetWallTime();
#endif // PROFILE_SMART_UPMAT
        const bool ok = determinant_generic_knowing_degree(det, pmat, target_degdet);
#ifdef PROFILE_SMART_UPMAT
        std::cout << "\ttotal of halving    --> " << GetWallTime()-t << std::endl;
        std::cout << "\ttotal of iteration  --> " << GetWallTime()-t_total << std::endl;
#endif // PROFILE_SMART_UPMAT
        return ok;
    }

    // find number of columns with degree d
    long block_size = 1;
    while (block_size < m && cdeg[m-block_size-1] == d)
        ++block_size;
    // ensure kernel requirement: m >= (d+1) * block_size
    if (m < (d+1)*block_size)
        block_size = m / (d+1);
    // note that m > d+1 here, so block_size >= 1

    // copy linearized version of block F of degree d
    // cmat = [ F_00 | F_01 | .. | F_0l ]  where l=d-1
    const long dn = d*block_size;
    Mat<zz_p> cmat;
    cmat.SetDims(m, dn);
    for (long k = 0; k < d; ++k)
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < block_size; ++j)
                cmat[i][k*block_size+j] = matp[k][i][m-block_size+j];

    // compute kernel
#ifdef PROFILE_SMART_UPMAT
    std::cout << "\tkernel ||  m, d, n, m-(d+1)n --> " << m << " x " << d << " x " << block_size << " x " << m-dn-block_size << std::endl;
    t_kernel = GetWallTime();
#endif // PROFILE_SMART_UPMAT
    Mat<zz_p> kertop, kerbot;
    kernel_degree1(kertop, kerbot, cmat, d);
#ifdef PROFILE_SMART_UPMAT
    t_kernel = GetWallTime() - t_kernel;
#endif // PROFILE_SMART_UPMAT
    // Recall this gives us constant matrices
    //     kertop = R, dimensions dn x dn 
    //     kerbot = K_0, dimensions ell x dn
    //     and now cmat = -F_0l, dimensions dn x n
    // where ell = m - (d+1)*n,  n = block_size
    // such that the Popov left kernel of the block F is
    //   [  R + x I_dn  |     0   | -F_0l ]
    //   [  -----------   -------   ----- ]
    //   [     K_0      |  I_ell  |   0   ]
    // Then left-multiplying matp by this kernel gives
    //                  [ M_0 + x^{cdeg1+1} I_dn |   M_1 + x^{cdeg2} L   |  0  ]
    //   ker * matp =   [ ---------------------- | --------------------- | --- ]
    //                  [          M_2           | M_3 + x^{cdeg2} I_ell |  0  ]
    // where
    //            cdeg = cdeg1 + cdeg2 + cdeg3
    //      with lengths -dn-    -ell-    -n-
    // and
    //         cdeg(M_0) < cdeg1+1, cdeg(M_2) < cdeg1+1
    //         cdeg(M_1) < cdeg2, cdeg(M_3) < cdeg2
    // and
    //         L = matp[ :dn, dn:dn+ell ]
    // ==> as a result, ignoring the zero columns, this has the shifted form
    // with the new cdeg = cdeg1+1  +  cdeg2,
    // UP TO modifying the top rows by removing the L part, which is simply done
    // by left-multiplication by the constant (ell+D)  x (ell+D)  matrix
    //       [ I_D  |   -L   ]
    //       [   0  |  I_ell ]
    // -> incorporating this into the kernel basis directly, we will left-multiply
    // by the matrix
    //   [  S + x I_dn  |   -L    | -F_0l ]
    //   [  -----------   -------   ----- ]
    //   [     K_0      |  I_ell  |   0   ]
    // where S = R - L K_0

    // discard columns of block F in matp
    Mat<zz_p> buf;
    long ell = m - block_size; // not yet the ell in comments above
    for (long k = 0; k < d; ++k)
    {
        buf.SetDims(m,ell);
        for (long i = 0; i < m; ++i)
            VectorCopy(buf[i], matp[k][i], ell);
        buf.swap(matp[k]);
    }

    ell = ell - dn; // now ell = m - (d+1)*block_size
    // compute S = R - L K_0
    Mat<zz_p> lmat;
    lmat.SetDims(dn, ell); // matrix -L
    for (long j = 0; j < ell; ++j)
        for (long i = 0; i < dn; ++i)
            lmat[i][j] = -matp[cdeg[dn+j]-1][i][dn+j];
#ifdef PROFILE_SMART_UPMAT
    t = GetWallTime();
#endif // PROFILE_SMART_UPMAT
    mul(buf, lmat, kerbot);  // now buf = -L K_0
    add(kertop, kertop, buf); // now kertop = S
#ifdef PROFILE_SMART_UPMAT
    t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT
    // form the row  [  S  |  -L  | -F_0l ] as defined above
    Mat<zz_p> trans_top;
    trans_top.SetDims(dn, m);
    for (long i = 0; i < dn; ++i)
    {
        for (long j = 0; j < dn; ++j)
            trans_top[i][j] = kertop[i][j];
        for (long j = dn; j < dn+ell; ++j)
            trans_top[i][j] = lmat[i][j-dn];
        for (long j = dn+ell; j < m; ++j)
            trans_top[i][j] = cmat[i][j-dn-ell];
    }

    // left multiply by kernel
    // update cdeg : discard last n entries ; add 1 to first dn entries
    VecLong(cdeg.begin(), cdeg.begin()+dn+ell).swap(cdeg);
    for (long i = 0; i < dn; ++i)
        cdeg[i] = cdeg[i] + 1;
    // update matp storage: add one coefficient, will be initialized in first iteration of loop
    matp.SetLength(cdeg[0]);
    // going through the coefficient of degree k from max to 0
    for (long k = cdeg[0]-2; k >= 0; --k)
    {
        const long cdim = matp[k].NumCols();
        // 1. multiply by X I_dn: add first rows to matp[k+1]
        // And update size of matp[k+1] if needed
        //  -- note that cdim2 <= cdim is always true
        //  -- and there is no need to resize if cdim2 == cdim || cdim2 >= dn
        const long cdim2 = matp[k+1].NumCols();
        if (cdim2 == 0)
        {
            // first iteration, k+1 == cdeg[0]-1
            matp[k+1].SetDims(dn+ell, std::min<long>(cdim,dn));
            for (long i = 0; i < dn; ++i)
                VectorCopy(matp[k+1][i], matp[k][i], std::min<long>(cdim,dn));
        }
        else if (cdim2 < std::min<long>(cdim,dn))
        {
            // matp[k+1] needs to be resized
            buf.swap(matp[k+1]);
            matp[k+1].SetDims(dn+ell,std::min<long>(cdim,dn)); // reallocates
            for (long i = 0; i < dn; ++i)
            {
#ifdef PROFILE_SMART_UPMAT
                t = GetWallTime();
#endif // PROFILE_SMART_UPMAT
                for (long j = 0; j < cdim2; ++j)
                    add(matp[k+1][i][j], buf[i][j], matp[k][i][j]);
#ifdef PROFILE_SMART_UPMAT
                t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT
                for (long j = cdim2; j < matp[k+1].NumCols(); ++j)
                    matp[k+1][i][j] = matp[k][i][j];
            }
            for (long i = dn; i < dn+ell; ++i)
                VectorCopy(matp[k+1][i], buf[i], matp[k+1].NumCols());
        }
        else // matp[k+1] does not need to be resized
        {
#ifdef PROFILE_SMART_UPMAT
            t = GetWallTime();
#endif // PROFILE_SMART_UPMAT
            for (long i = 0; i < dn; ++i)
                for (long j = 0; j < cdim2; ++j)
                    add(matp[k+1][i][j], matp[k+1][i][j], matp[k][i][j]);
#ifdef PROFILE_SMART_UPMAT
            t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT
        }

        // 2. multiply by constant kernel    [ K_0 | I_ell | 0 ]
        // and store in buf
        buf.kill(); // TODO remove?
        buf.SetDims(dn, cdim);
        for (long i = 0; i < dn; ++i)
            VectorCopy(buf[i], matp[k][i], cdim);
#ifdef PROFILE_SMART_UPMAT
        t = GetWallTime();
#endif // PROFILE_SMART_UPMAT
        mul(buf, kerbot, buf);
        for (long i = 0; i < ell; ++i)
            add(buf[i], buf[i], matp[k][i+dn]);
        // 3. multiply by degree1 kernel   trans_top = [  S | -L | -F_0l ]
        // and store in matp[k]
        mul(matp[k], trans_top, matp[k]);
#ifdef PROFILE_SMART_UPMAT
        t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT
        // 4. stack rows of buf below this product
        matp[k].SetDims(dn+ell,cdim);
        for (long i = dn; i < dn+ell; ++i)
            matp[k][i].swap(buf[i-dn]);
    }
    // finally don't forget to add 
    //   [  S  ]
    //   [ K_0 ]
    // in leading terms of leftmost rows since it is multiplied by the
    // (non-stored) X^cdeg1 I_dn
    // FIXME if needed, there surely is a more cache friendly way
    // to do this, by incorporating it into the above loop
#ifdef PROFILE_SMART_UPMAT
        t = GetWallTime();
#endif // PROFILE_SMART_UPMAT
    for (long j = 0; j < dn; ++j)
    {
        for (long i = 0; i < dn; ++i)
            matp[cdeg[j]-1][i][j] += trans_top[i][j];
        for (long i = dn; i < dn+ell; ++i)
            matp[cdeg[j]-1][i][j] += kerbot[i-dn][j];
    }
#ifdef PROFILE_SMART_UPMAT
    t_update += GetWallTime() - t;
    t_total = GetWallTime() - t_total;
    std::cout << "End of iteration; times before recursion:" << std::endl;
    std::cout << "            kernel --> " << t_kernel << std::endl;
    std::cout << "            update --> " << t_update << std::endl;
    std::cout << "            total  --> " << t_total << std::endl;
    std::cout << "  100*(ker+up)/tot --> " << 100*(t_kernel+t_update)/t_total << std::endl;
#endif // PROFILE_SMART_UPMAT

    return determinant_shifted_form_smartkernel_updateall(det, matp, cdeg, threshold, target_degdet);
}

// det: output determinant
// matp: input square m x m matrix
// cdeg: column degrees of the input matp (len(cdeg) = m), assumed nondecreasing
// threshold: when "??? TODO", switch to a classical determinant algorithm
//      where d is the smallest entry (i.e. last entry)
// target_degdet: degree of determinant of matp, used for certification
// REQUIREMENTS:
//   -- matp is in "shiftedform": matp = X^{cdeg} Id + R, where cdeg(R) < cdeg
//   -- cdeg is nonincreasing   (cdeg[j+1] <= cdeg[j])
//   -- storage: matp is a list of coefficient matrices; matp[i] has dimensions m x n_i
//      where n_i is the largest integer j such that cdeg[j] >= i
//      NOTE: the X^cdeg Id part is not stored
//   -- target_degdet: must be equal to deg(det(matp)), otherwise the
//   behaviour of certification is not guaranteed
// GUARANTEES IN OUTPUT:
//   -- det = det(matp)
//   -- if target_degdet in input, then the return value is true iff det = det(matp)
// PROPERTIES:
//   -- degree of determinant is sum(cdeg) = sum of diagonal degrees
bool determinant_shifted_form_smartkernel_updateall_v2(zz_pX & det, Vec<Mat<zz_p>> & matp, VecLong & cdeg, long threshold, long target_degdet)
{
#ifdef DEBUGGING_NOW
    std::cout << target_degdet << std::endl;
    std::cout << cdeg << std::endl;
#endif // DEBUGGING_NOW
    
#ifdef PROFILE_SMART_UPMAT_V2
    double t_total, t_kernel, t_update, t;
    t_update = 0;
    t_total = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
    if (matp.length() == 0) // identity matrix
    {
        det = 1;
#ifdef PROFILE_SMART_UPMAT_V2
    std::cout << "\tbase case identity matrix" << std::endl;
#endif // PROFILE_SMART_UPMAT_V2
        return (target_degdet == 0);
    }

    const long m = matp[0].NumRows();
#ifdef PROFILE_SMART_UPMAT_V2
    std::cout << "Entering iteration with m , mindeg = " << m << " , " << cdeg[m-1] << std::endl;
#endif // PROFILE_SMART_UPMAT_V2
    if (m != (long)cdeg.size())
    {
        std::cout << "~~ERROR~~ wrong length of input cdeg in determinant_shifted_form_smartkernel_updateall" << std::endl;
        return false;
    }

    // if sum(cdeg) is not target_degree, then something went wrong in earlier steps
    if (std::accumulate(cdeg.begin(), cdeg.end(), 0) != target_degdet)
    {
        std::cout << "~~ERROR~~  in determinant_shifted_form_smartkernel_updateall" << std::endl;
        std::cout << "         --> target_degdet is not equal to sum(cdeg)" << std::endl;
        std::cout << "         --> either cdeg is wrong" << std::endl;
        std::cout << "           (which is not supposed to happen if it was correct in input)" << std::endl;
        std::cout << "         --> or some kernel computation went wrong" << std::endl;
        std::cout << "           (which can happen due to genericity assumption not satisfied)" << std::endl;
        return false;
    }

    // if small dim, just run expansion by minors
    if (m<=4)
    {
        Mat<zz_pX> pmat;
        conv_shifted_form(pmat, matp, cdeg);
#ifdef PROFILE_SMART_UPMAT_V2
        t=GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
        determinant_expansion_by_minors(det, pmat);
#ifdef PROFILE_SMART_UPMAT_V2
        std::cout << "\tbase case (expansion) --> " << GetWallTime()-t << std::endl;
        std::cout << "\ttotal of iteration    --> " << GetWallTime()-t_total << std::endl;
#endif // PROFILE_SMART_UPMAT_V2
        return true;
    }

    // above some threshold,
    // run the usual algo splitting column dimension in two equal parts
    if (cdeg[m-1] >= threshold || m <= cdeg[m-1]+1) // TODO
    {
#ifdef PROFILE_SMART_UPMAT_V2
        std::cout << "\t-->Entering halving stage" << std::endl;
#endif // PROFILE_SMART_UPMAT_V2
        Mat<zz_pX> pmat;
        conv_shifted_form_mirrored(pmat, matp, cdeg);

#ifdef PROFILE_SMART_UPMAT_V2
        t=GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
        const bool ok = determinant_generic_knowing_degree(det, pmat, target_degdet);
#ifdef PROFILE_SMART_UPMAT_V2
        std::cout << "\ttotal of halving    --> " << GetWallTime()-t << std::endl;
        std::cout << "\ttotal of iteration  --> " << GetWallTime()-t_total << std::endl;
#endif // PROFILE_SMART_UPMAT_V2
        return ok;
    }

    // take nn = largest possible number of columns handled in this step
    // the requirement is  m >= D + nn
    long nn = 1;          // size of the block
    long D = cdeg[m-nn];  // sum of cdeg of the block
    long mu = cdeg[m-nn]; // degree of the block
    while (nn < m && D+cdeg[m-nn-1]+nn < m)
    {
        ++nn;
        D += cdeg[m-nn];
        mu = cdeg[m-nn];
    }

    // compute kernel
#ifdef PROFILE_SMART_UPMAT_V2
    std::cout << "\tkernel ||  m, D, nn, ell, mu --> " << m << " x " << D << " x " << nn << " x " << m-D-nn << " x " << mu << std::endl;
    t_kernel = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
    Mat<zz_p> K0, K1, lF0;
    kernel_degree1_cdeg(K0, K1, lF0, matp, cdeg, m-nn);
#ifdef PROFILE_SMART_UPMAT_V2
    t_kernel = GetWallTime() - t_kernel;
#endif // PROFILE_SMART_UPMAT_V2
    // Recall this gives us constant matrices
    //     K0, dimensions D x D 
    //     K1, dimensions ell x D
    //     lF0, dimensions D x n
    // where ell = m - D - nn
    // such that the Popov left kernel of the considered block is
    //   [ x I_D + K0 |    0    | lF0 ]
    //   [ ----------   -------   ----- ]
    //   [     K1     |  I_ell  |   0   ]
    // Then left-multiplying matp by this kernel gives
    //                  [ M_0 + x^{cdeg1+1} I_D |   M_1 + x^{cdeg2} L   |  0  ]
    //   ker * matp =   [ ---------------------- | --------------------- | --- ]
    //                  [          M_2           | M_3 + x^{cdeg2} I_ell |  0  ]
    // where
    //            cdeg = cdeg1 + cdeg2 + cdeg3
    //      with lengths -D-    -ell-    -n-
    // and
    //         cdeg(M_0) < cdeg1+1, cdeg(M_2) < cdeg1+1
    //         cdeg(M_1) < cdeg2, cdeg(M_3) < cdeg2
    // and
    //         L = matp[ :D, D:D+ell ]
    // ==> as a result, ignoring the zero columns, this has the shifted form
    // with the new cdeg = cdeg1+1  +  cdeg2,
    // UP TO modifying the top rows by removing the L part, which is simply done
    // by left-multiplication by the constant (ell+D)  x (ell+D)  matrix
    //       [ I_D  |   -L   ]
    //       [   0  |  I_ell ]
    // -> incorporating this into the kernel basis directly, we will left-multiply
    // by the matrix
    //   [  S + x I_D  |   -L    | lF0 ]
    //   [  ----------   -------   --- ]
    //   [      K1     |  I_ell  |  0  ]
    // where S = K0 - L K1

    // discard columns of block F in matp
    Mat<zz_p> buf;
    long ell = m - nn; // warning: not yet the ell in comments above
    for (long k = 0; k < mu; ++k)
    {
        buf.SetDims(m,ell);
        for (long i = 0; i < m; ++i)
            VectorCopy(buf[i], matp[k][i], ell);
        buf.swap(matp[k]);
    }

    ell = ell - D; // now ell is indeed ell = m - D - nn
    // compute S = K0 - L K1 and store it in K0
    Mat<zz_p> lmat;
    lmat.SetDims(D, ell); // matrix -L
    for (long j = 0; j < ell; ++j)
        for (long i = 0; i < D; ++i)
            lmat[i][j] = -matp[cdeg[D+j]-1][i][D+j];
#ifdef PROFILE_SMART_UPMAT_V2
    t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
    mul(buf, lmat, K1);  // now buf = -L K1
    add(K0, K0, buf); // now K0 = S
#ifdef PROFILE_SMART_UPMAT_V2
    t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V2
    // form the row  [  S  |  -L  | lF0 ] as defined above
    Mat<zz_p> trans_top;
    trans_top.SetDims(D, m);
    for (long i = 0; i < D; ++i)
    {
        for (long j = 0; j < D; ++j)
            trans_top[i][j] = K0[i][j];
        for (long j = D; j < D+ell; ++j)
            trans_top[i][j] = lmat[i][j-D];
        for (long j = D+ell; j < m; ++j)
            trans_top[i][j] = lF0[i][j-D-ell];
    }

    // left multiply by kernel
    // update cdeg : discard last nn entries ; add 1 to first D entries
    VecLong(cdeg.begin(), cdeg.begin()+D+ell).swap(cdeg);
    for (long i = 0; i < D; ++i)
        cdeg[i] = cdeg[i] + 1;
    // update matp storage: add one coefficient, will be initialized in first iteration of loop
    matp.SetLength(cdeg[0]);
    // going through the coefficient of degree k from max to 0
    for (long k = cdeg[0]-2; k >= 0; --k)
    {
        const long cdim = matp[k].NumCols();
        // 1. multiply by X I_D: add first rows to matp[k+1]
        // And update size of matp[k+1] if needed
        //  -- note that cdim2 <= cdim is always true
        //  -- and there is no need to resize if cdim2 == cdim || cdim2 >= D
        const long cdim2 = matp[k+1].NumCols();
        if (cdim2 == 0)
        {
            // first iteration, k+1 == cdeg[0]-1
            matp[k+1].SetDims(D+ell, std::min<long>(cdim,D));
            for (long i = 0; i < D; ++i)
                VectorCopy(matp[k+1][i], matp[k][i], std::min<long>(cdim,D));
        }
        else if (cdim2 < std::min<long>(cdim,D))
        {
            // matp[k+1] needs to be resized
            buf.swap(matp[k+1]);
            matp[k+1].SetDims(D+ell,std::min<long>(cdim,D)); // reallocates
            for (long i = 0; i < D; ++i)
            {
#ifdef PROFILE_SMART_UPMAT_V2
                t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
                for (long j = 0; j < cdim2; ++j)
                    add(matp[k+1][i][j], buf[i][j], matp[k][i][j]);
#ifdef PROFILE_SMART_UPMAT_V2
                t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V2
                for (long j = cdim2; j < matp[k+1].NumCols(); ++j)
                    matp[k+1][i][j] = matp[k][i][j];
            }
            for (long i = D; i < D+ell; ++i)
                VectorCopy(matp[k+1][i], buf[i], matp[k+1].NumCols());
        }
        else // matp[k+1] does not need to be resized
        {
#ifdef PROFILE_SMART_UPMAT_V2
            t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
            for (long i = 0; i < D; ++i)
                for (long j = 0; j < cdim2; ++j)
                    add(matp[k+1][i][j], matp[k+1][i][j], matp[k][i][j]);
#ifdef PROFILE_SMART_UPMAT_V2
            t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V2
        }

        // 2. multiply by constant kernel    [ K1 | I_ell | 0 ]
        // and store in buf
        buf.kill(); // TODO remove?
        buf.SetDims(D, cdim);
        for (long i = 0; i < D; ++i)
            VectorCopy(buf[i], matp[k][i], cdim);
#ifdef PROFILE_SMART_UPMAT_V2
        t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
        mul(buf, K1, buf);
        for (long i = 0; i < ell; ++i)
            add(buf[i], buf[i], matp[k][i+D]);
        // 3. multiply by degree1 kernel   trans_top = [  S | -L | lF0 ]
        // and store in matp[k]
        mul(matp[k], trans_top, matp[k]);
#ifdef PROFILE_SMART_UPMAT_V2
        t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V2
        // 4. stack rows of buf below this product
        matp[k].SetDims(D+ell,cdim);
        for (long i = D; i < D+ell; ++i)
            matp[k][i].swap(buf[i-D]);
    }
    // finally don't forget to add 
    //   [ S  ]
    //   [ K1 ]
    // in leading terms of leftmost rows since it is multiplied by the
    // (non-stored) X^cdeg1 I_D
    // FIXME if needed, there surely is a more cache friendly way to do this...
#ifdef PROFILE_SMART_UPMAT_V2
        t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V2
    for (long j = 0; j < D; ++j)
    {
        for (long i = 0; i < D; ++i)
            matp[cdeg[j]-1][i][j] += trans_top[i][j];
        for (long i = D; i < D+ell; ++i)
            matp[cdeg[j]-1][i][j] += K1[i-D][j];
    }
#ifdef PROFILE_SMART_UPMAT_V2
    t_update += GetWallTime() - t;
    t_total = GetWallTime() - t_total;
    std::cout << "End of iteration; times before recursion:" << std::endl;
    std::cout << "          kernel --> " << t_kernel << std::endl;
    std::cout << "          update --> " << t_update << std::endl;
    std::cout << "          total  --> " << t_total << std::endl;
    std::cout << "    non-dominant --> " << std::setprecision(2) << 100*(t_kernel+t_update)/t_total << "%" << std::setprecision(8) << std::endl;
#endif // PROFILE_SMART_UPMAT_V2

    return determinant_shifted_form_smartkernel_updateall_v2(det, matp, cdeg, threshold, target_degdet);
}

// det: output determinant
// matp: input square m x m matrix
// cdeg: column degrees of the input matp (len(cdeg) = m), assumed nondecreasing
// threshold: when "??? TODO", switch to a classical determinant algorithm
//      where d is the smallest entry (i.e. last entry)
// target_degdet: degree of determinant of matp, used for certification
// REQUIREMENTS:
//   -- matp is in "shiftedform": matp = X^{cdeg} Id + R, where cdeg(R) < cdeg
//   -- cdeg is nonincreasing   (cdeg[j+1] <= cdeg[j])
//   -- storage: matp is a list of coefficient matrices; matp[i] has dimensions m x n_i
//      where n_i is the largest integer j such that cdeg[j] >= i
//      NOTE: the X^cdeg Id part is not stored
//   -- target_degdet: must be equal to deg(det(matp)), otherwise the
//   behaviour of certification is not guaranteed
// GUARANTEES IN OUTPUT:
//   -- det = det(matp)
//   -- if target_degdet in input, then the return value is true iff det = det(matp)
// PROPERTIES:
//   -- degree of determinant is sum(cdeg) = sum of diagonal degrees
bool determinant_shifted_form_smartkernel_updateall_v3(zz_pX & det, Vec<Mat<zz_p>> & matp, VecLong & cdeg, long threshold, long target_degdet)
{
#ifdef DEBUGGING_NOW
    std::cout << target_degdet << std::endl;
    std::cout << cdeg << std::endl;
#endif // DEBUGGING_NOW
    
#ifdef PROFILE_SMART_UPMAT_V3
    double t_total, t_kernel, t_update, t;
    t_update = 0;
    t_total = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
    if (matp.length() == 0) // identity matrix
    {
        det = 1;
#ifdef PROFILE_SMART_UPMAT_V3
    std::cout << "\tbase case identity matrix" << std::endl;
#endif // PROFILE_SMART_UPMAT_V3
        return (target_degdet == 0);
    }

    const long m = matp[0].NumRows();
#ifdef PROFILE_SMART_UPMAT_V3
    std::cout << "Entering iteration with m , mindeg = " << m << " , " << cdeg[m-1] << std::endl;
#endif // PROFILE_SMART_UPMAT_V3
    if (m != (long)cdeg.size())
    {
        std::cout << "~~ERROR~~ wrong length of input cdeg in determinant_shifted_form_smartkernel_updateall" << std::endl;
        return false;
    }

    // if sum(cdeg) is not target_degree, then something went wrong in earlier steps
    if (std::accumulate(cdeg.begin(), cdeg.end(), 0) != target_degdet)
    {
        std::cout << "~~ERROR~~  in determinant_shifted_form_smartkernel_updateall" << std::endl;
        std::cout << "         --> target_degdet is not equal to sum(cdeg)" << std::endl;
        std::cout << "         --> either cdeg is wrong" << std::endl;
        std::cout << "           (which is not supposed to happen if it was correct in input)" << std::endl;
        std::cout << "         --> or some kernel computation went wrong" << std::endl;
        std::cout << "           (which can happen due to genericity assumption not satisfied)" << std::endl;
        return false;
    }

    // if small dim, just run expansion by minors
    if (m<=4)
    {
        Mat<zz_pX> pmat;
        conv_shifted_form(pmat, matp, cdeg);
#ifdef PROFILE_SMART_UPMAT_V3
        t=GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
        determinant_expansion_by_minors(det, pmat);
#ifdef PROFILE_SMART_UPMAT_V3
        std::cout << "\tbase case (expansion) --> " << GetWallTime()-t << std::endl;
        std::cout << "\ttotal of iteration    --> " << GetWallTime()-t_total << std::endl;
#endif // PROFILE_SMART_UPMAT_V3
        return true;
    }

    // above some threshold,
    // run the usual algo splitting column dimension in two equal parts
    if (cdeg[m-1] >= threshold || m <= cdeg[m-1]+1) // TODO
    {
#ifdef PROFILE_SMART_UPMAT_V3
        std::cout << "\t-->Entering halving stage" << std::endl;
#endif // PROFILE_SMART_UPMAT_V3
        Mat<zz_pX> pmat;
        conv_shifted_form_mirrored(pmat, matp, cdeg);

#ifdef PROFILE_SMART_UPMAT_V3
        t=GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
        const bool ok = determinant_generic_knowing_degree(det, pmat, target_degdet);
#ifdef PROFILE_SMART_UPMAT_V3
        std::cout << "\ttotal of halving    --> " << GetWallTime()-t << std::endl;
        std::cout << "\ttotal of iteration  --> " << GetWallTime()-t_total << std::endl;
#endif // PROFILE_SMART_UPMAT_V3
        return ok;
    }

    // take nn = largest possible number of columns handled in this step
    // the requirement is  m >= D + nn
    long nn = 1;          // size of the block
    long D = cdeg[m-nn];  // sum of cdeg of the block
    long mu = cdeg[m-nn]; // degree of the block
    while (nn < m && D+cdeg[m-nn-1]+nn < m)
    {
        ++nn;
        D += cdeg[m-nn];
        mu = cdeg[m-nn];
    }

    // compute kernel
#ifdef PROFILE_SMART_UPMAT_V3
    std::cout << "\tkernel ||  m, D, nn, ell, mu --> " << m << " x " << D << " x " << nn << " x " << m-D-nn << " x " << mu << std::endl;
    t_kernel = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
    Mat<zz_p> K0, K1, lF0;
    kernel_degree1_cdeg(K0, K1, lF0, matp, cdeg, m-nn);
#ifdef PROFILE_SMART_UPMAT_V3
    t_kernel = GetWallTime() - t_kernel;
#endif // PROFILE_SMART_UPMAT_V3
    // Recall this gives us constant matrices
    //     K0, dimensions D x D 
    //     K1, dimensions ell x D
    //     lF0, dimensions D x n
    // where ell = m - D - nn
    // such that the Popov left kernel of the considered block is
    //   [ x I_D + K0 |    0    | lF0 ]
    //   [ ----------   -------   ----- ]
    //   [     K1     |  I_ell  |   0   ]
    // Then left-multiplying matp by this kernel gives
    //                  [ M_0 + x^{cdeg1+1} I_D |   M_1 + x^{cdeg2} L   |  0  ]
    //   ker * matp =   [ ---------------------- | --------------------- | --- ]
    //                  [          M_2           | M_3 + x^{cdeg2} I_ell |  0  ]
    // where
    //            cdeg = cdeg1 + cdeg2 + cdeg3
    //      with lengths -D-    -ell-    -n-
    // and
    //         cdeg(M_0) < cdeg1+1, cdeg(M_2) < cdeg1+1
    //         cdeg(M_1) < cdeg2, cdeg(M_3) < cdeg2
    // and
    //         L = matp[ :D, D:D+ell ]
    // ==> as a result, ignoring the zero columns, this has the shifted form
    // with the new cdeg = cdeg1+1  +  cdeg2,
    // UP TO modifying the top rows by removing the L part, which is simply done
    // by left-multiplication by the constant (ell+D)  x (ell+D)  matrix
    //       [ I_D  |   -L   ]
    //       [   0  |  I_ell ]
    // -> incorporating this into the kernel basis directly, we will left-multiply
    // by the matrix
    //   [  S + x I_D  |   -L    | lF0 ]
    //   [  ----------   -------   --- ]
    //   [      K1     |  I_ell  |  0  ]
    // where S = K0 - L K1

    // discard columns of block F in matp
    Mat<zz_p> buf;
    long ell = m - nn; // warning: not yet the ell in comments above
    for (long k = 0; k < mu; ++k)
    {
        buf.SetDims(m,ell);
        for (long i = 0; i < m; ++i)
            VectorCopy(buf[i], matp[k][i], ell);
        buf.swap(matp[k]);
    }

    ell = ell - D; // now ell is indeed ell = m - D - nn
    // compute S = K0 - L K1 and store it in K0
    Mat<zz_p> lmat;
    lmat.SetDims(D, ell); // matrix -L
    for (long j = 0; j < ell; ++j)
        for (long i = 0; i < D; ++i)
            lmat[i][j] = -matp[cdeg[D+j]-1][i][D+j];
#ifdef PROFILE_SMART_UPMAT_V3
    t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
    mul(buf, lmat, K1);  // now buf = -L K1
    add(K0, K0, buf); // now K0 = S
#ifdef PROFILE_SMART_UPMAT_V3
    t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V3
    // form the row  [  S  |  -L  | lF0 ] as defined above
    Mat<zz_p> trans_top;
    trans_top.SetDims(D, m);
    for (long i = 0; i < D; ++i)
    {
        for (long j = 0; j < D; ++j)
            trans_top[i][j] = K0[i][j];
        for (long j = D; j < D+ell; ++j)
            trans_top[i][j] = lmat[i][j-D];
        for (long j = D+ell; j < m; ++j)
            trans_top[i][j] = lF0[i][j-D-ell];
    }

    // left multiply by kernel
    // update cdeg : discard last nn entries ; add 1 to first D entries
    VecLong(cdeg.begin(), cdeg.begin()+D+ell).swap(cdeg);
    for (long i = 0; i < D; ++i)
        cdeg[i] = cdeg[i] + 1;
    // update matp storage: add one coefficient, will be initialized in first iteration of loop
    matp.SetLength(cdeg[0]);
    // going through the coefficient of degree k from max to 0
    for (long k = cdeg[0]-2; k >= 0; --k)
    {
        const long cdim = matp[k].NumCols();
        // 1. multiply by X I_D: add first rows to matp[k+1]
        // And update size of matp[k+1] if needed
        //  -- note that cdim2 <= cdim is always true
        //  -- and there is no need to resize if cdim2 == cdim || cdim2 >= D
        const long cdim2 = matp[k+1].NumCols();
        if (cdim2 == 0)
        {
            // first iteration, k+1 == cdeg[0]-1
            matp[k+1].SetDims(D+ell, std::min<long>(cdim,D));
            for (long i = 0; i < D; ++i)
                VectorCopy(matp[k+1][i], matp[k][i], std::min<long>(cdim,D));
        }
        else if (cdim2 < std::min<long>(cdim,D))
        {
            // matp[k+1] needs to be resized
            buf.swap(matp[k+1]);
            matp[k+1].SetDims(D+ell,std::min<long>(cdim,D)); // reallocates
            for (long i = 0; i < D; ++i)
            {
#ifdef PROFILE_SMART_UPMAT_V3
                t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
                for (long j = 0; j < cdim2; ++j)
                    add(matp[k+1][i][j], buf[i][j], matp[k][i][j]);
#ifdef PROFILE_SMART_UPMAT_V3
                t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V3
                for (long j = cdim2; j < matp[k+1].NumCols(); ++j)
                    matp[k+1][i][j] = matp[k][i][j];
            }
            for (long i = D; i < D+ell; ++i)
                VectorCopy(matp[k+1][i], buf[i], matp[k+1].NumCols());
        }
        else // matp[k+1] does not need to be resized
        {
#ifdef PROFILE_SMART_UPMAT_V3
            t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
            for (long i = 0; i < D; ++i)
                for (long j = 0; j < cdim2; ++j)
                    add(matp[k+1][i][j], matp[k+1][i][j], matp[k][i][j]);
#ifdef PROFILE_SMART_UPMAT_V3
            t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V3
        }

        // 2. multiply by constant kernel    [ K1 | I_ell | 0 ]
        // and store in buf
        buf.kill(); // TODO remove?
        buf.SetDims(D, cdim);
        for (long i = 0; i < D; ++i)
            VectorCopy(buf[i], matp[k][i], cdim);
#ifdef PROFILE_SMART_UPMAT_V3
        t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
        mul(buf, K1, buf);
        for (long i = 0; i < ell; ++i)
            add(buf[i], buf[i], matp[k][i+D]);
        // 3. multiply by degree1 kernel   trans_top = [  S | -L | lF0 ]
        // and store in matp[k]
        mul(matp[k], trans_top, matp[k]);
#ifdef PROFILE_SMART_UPMAT_V3
        t_update += GetWallTime() - t;
#endif // PROFILE_SMART_UPMAT_V3
        // 4. stack rows of buf below this product
        matp[k].SetDims(D+ell,cdim);
        for (long i = D; i < D+ell; ++i)
            matp[k][i].swap(buf[i-D]);
    }
    // finally don't forget to add 
    //   [ S  ]
    //   [ K1 ]
    // in leading terms of leftmost rows since it is multiplied by the
    // (non-stored) X^cdeg1 I_D
    // FIXME if needed, there surely is a more cache friendly way to do this...
#ifdef PROFILE_SMART_UPMAT_V3
        t = GetWallTime();
#endif // PROFILE_SMART_UPMAT_V3
    for (long j = 0; j < D; ++j)
    {
        for (long i = 0; i < D; ++i)
            matp[cdeg[j]-1][i][j] += trans_top[i][j];
        for (long i = D; i < D+ell; ++i)
            matp[cdeg[j]-1][i][j] += K1[i-D][j];
    }
#ifdef PROFILE_SMART_UPMAT_V3
    t_update += GetWallTime() - t;
    t_total = GetWallTime() - t_total;
    std::cout << "End of iteration; times before recursion:" << std::endl;
    std::cout << "          kernel --> " << t_kernel << std::endl;
    std::cout << "          update --> " << t_update << std::endl;
    std::cout << "          total  --> " << t_total << std::endl;
    std::cout << "    non-dominant --> " << std::setprecision(2) << 100*(t_kernel+t_update)/t_total << "%" << std::setprecision(8) << std::endl;
#endif // PROFILE_SMART_UPMAT_V3

    return determinant_shifted_form_smartkernel_updateall_v2(det, matp, cdeg, threshold, target_degdet);
}



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
        //zz_p::init(29);
        //cout << zz_p::modulus() << ", bit length = " << nbits << endl;
    }

    // retrieve degree matrix from file (this one has decreasing diag degrees)
    Mat<long> dmat;
    Vec<long> cdeg;
    retrieve_degree_matrix(dmat,cdeg,filename);
    const long dim = dmat.NumRows();
    const long degdet = std::accumulate(cdeg.begin(),cdeg.end(),0);
    //// small example
    //dmat.SetDims(8,8);
    //cdeg.SetLength(8);
    //cdeg[0] = 5; cdeg[1] = 3; cdeg[2] = 3; cdeg[3] = 3; cdeg[4] = 3; cdeg[5] = 3; cdeg[6] = 1; cdeg[7] = 1;
    //for (long i = 0; i < 8; ++i)
    //    for (long j = 0; j < 8; ++j)
    //        dmat[i][j] = cdeg[j] - 1;
    //for (long i = 0; i < 8; ++i)
    //    dmat[i][i] = cdeg[i];
    //dim = dmat.NumRows();
    //degdet = std::accumulate(cdeg.begin(),cdeg.end(),0);

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
#if true
    std::cout << std::endl;
    std::cout << "sequence of degrees:\t";
    std::copy(unique_cdeg.begin(), unique_cdeg.end(), std::ostream_iterator<long>(std::cout, "\t"));
    std::cout << std::endl;
    std::cout << "ncols for each degree:\t";
    std::copy(mult_cdeg.begin(), mult_cdeg.end(), std::ostream_iterator<long>(std::cout, "\t"));
    std::cout << std::endl << "proportion of total:";
    for (size_t k = 0; k < mult_cdeg.size(); ++k)
        std::cout << "\t" << std::setprecision(3) << 100 * mult_cdeg[k] / (double)dim << std::setprecision(8);
    std::cout << std::endl;
#endif

    std::vector<double> timings;
    double t,tt;
    long nb_iter;

#ifdef TIME_LINSOLVE
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
            std::cout << "~~~Warning~~~ verification of determinant failed in LINSOLVE approach" << std::endl;
    }
#endif

#ifdef TIME_TRI_DECR
    { // generic determinant, decreasing diagonal degrees
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
            std::cout << "~~~Warning~~~ verification of determinant failed in TRI_DECR approach" << std::endl;
    }
#endif

#ifdef TIME_TRI_INCR
    { // generic determinant, increasing diagonal degrees
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
            std::cout << "~~~Warning~~~ verification of determinant failed in TRI_INCR approach" << std::endl;
    }
#endif

#ifdef TIME_DEG_UPMAT
    for (long thres=4; thres < 6; ++thres)
    { // shifted form specific, degree-aware, update matrix
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
            std::cout << "~~~Warning~~~ verification of determinant failed in DEG_UPMAT approach" << std::endl;
    }
#endif

#ifdef TIME_DEG_UPKER
    for (long thres=3; thres < 4; ++thres)
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
            std::cout << "~~~Warning~~~ verification of determinant failed in DEG_UPKER approach" << std::endl;
    }
#endif

#ifdef TIME_SMART_UPMAT
    for (long thres=20; thres < 101; thres+=10)
    { // shifted form specific, degree-aware, update matrix, smart kernel
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            //std::cout << degdet << std::endl;
            VecLong copy_cdeg(dim);
            for (long j = 0; j < dim; ++j)
                copy_cdeg[j] = cdeg[j];
            Vec<Mat<zz_p>> matp;
            matp.SetLength(cdeg[0]);
            for (long k = 0; k < cdeg[0]; ++k)
            {
                // find size
                long n_k = 0;
                while (n_k < cdeg.length() && cdeg[n_k] > k)
                    ++n_k;
                // random dim x n_k constant term (that of degree k)
                random(matp[k], dim, n_k);
            }
            Mat<zz_pX> pmat;
            conv_shifted_form(pmat,matp,copy_cdeg);
            //for (long i = 0; i < pmat.NumRows(); ++i)
            //{
            //    for (long j = 0; j < pmat.NumCols(); ++j)
            //    {
            //        std::cout << "pR([" ;
            //        for (long k = 0; k <= deg(pmat[i][j]); ++k)
            //        {
            //            std::cout << pmat[i][j][k] << ", ";
            //        }
            //        std::cout << "])," << std::endl;
            //    }
            //}
            //std::cout << matp << std::endl;
            //std::cout << pmat << std::endl;
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_shifted_form_smartkernel_updateall(det, matp, copy_cdeg, thres, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in KER_UPMAT approach" << std::endl;
    }
#endif

#ifdef TIME_SMART_UPMAT_V2
    std::cout << std::endl;
    for (long thres=30; thres < 71; thres+=10)
    { // shifted form specific, degree-aware, update matrix, smart kernel, version 2
        t=0.0; nb_iter=0;
        bool ok = true;
        while (ok && t<1)
        {
            VecLong copy_cdeg(dim);
            for (long j = 0; j < dim; ++j)
                copy_cdeg[j] = cdeg[j];
            Vec<Mat<zz_p>> matp;
            matp.SetLength(cdeg[0]);
            for (long k = 0; k < cdeg[0]; ++k)
            {
                // find size
                long n_k = 0;
                while (n_k < cdeg.length() && cdeg[n_k] > k)
                    ++n_k;
                // random dim x n_k constant term (that of degree k)
                random(matp[k], dim, n_k);
            }
            Mat<zz_pX> pmat;
            conv_shifted_form(pmat,matp,copy_cdeg);
            //for (long i = 0; i < pmat.NumRows(); ++i)
            //{
            //    for (long j = 0; j < pmat.NumCols(); ++j)
            //    {
            //        std::cout << "pR([" ;
            //        for (long k = 0; k <= deg(pmat[i][j]); ++k)
            //        {
            //            std::cout << pmat[i][j][k] << ", ";
            //        }
            //        std::cout << "])," << std::endl;
            //    }
            //}
            //std::cout << matp << std::endl;
            //std::cout << pmat << std::endl;
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_shifted_form_smartkernel_updateall_v2(det, matp, copy_cdeg, thres, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        timings.push_back(t/nb_iter);
        std::cout << "thresh: " << thres << "\t timing: " << t/nb_iter << std::endl;
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in KER_UPMAT approach" << std::endl;
    }
#endif


#ifdef TIME_ZLS_AVG
    { // shifted form specific, degree-aware, update one
        t=0.0; nb_iter=0;
        bool ok = true;
#ifdef DEBUGGING_NOW
        while (ok && t<1 && nb_iter<1)
#else
        while (ok && t<1)
#endif
        {
            Mat<zz_pX> pmat,kerbas;
            random(pmat, dmat2);
            VecLong diag_deg = col_degree(pmat);
            VecLong shifted_rdeg(dim);
            VecLong split_sizes(mult_cdeg2);
            tt = GetWallTime();
            zz_pX det;
            ok = ok && determinant_shifted_form_zls_1(det, pmat, diag_deg, shifted_rdeg, split_sizes, degdet);
            t += GetWallTime()-tt;
            ++nb_iter;
            ok = ok && verify_determinant(det, pmat, true, true);
        }
        timings.push_back(t/nb_iter);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in ZLS_AVG approach" << std::endl;
    }
#endif

#ifdef ESTIMATE_WIEDEMANN
    { // estimate Wiedemann
        t=0.0; nb_iter=0;
        while (t<1)
        {
            Mat<zz_p> mat;
            Vec<zz_p> vec,buf;
            random(mat, dim, degdet);
            random(vec, degdet);
            tt = GetWallTime();
            for (long d = 0; d < 2*degdet; ++d)
            {
                mul(buf, mat, vec);
                for (long i = degdet-dim; i < degdet; ++i)
                    vec[i] = buf[i];
            }
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        timings.push_back(t/nb_iter);
    }
#endif // ESTIMATE_WIEDEMANN

#ifdef ESTIMATE_KELLERGEHRIG
    { // estimate Keller-Gehrig
        t=0.0; nb_iter=0;
        while (t<1)
        {
            // working in columns to avoid introducing artificial
            // slowliness due to data movements
            Mat<zz_p> mat;
            Mat<zz_p> projs,buf;
            ident(mat, degdet);
            for (long i = degdet-dim; i < degdet; ++i)
                random(mat[i], degdet);
            // transpose, work in columns
            transpose(mat, mat);
            random(projs, 1, degdet);
            tt = GetWallTime();
            for (long d = 1; d <= degdet; d=2*d)
            {
                // mat = power(M , d) where M = input matrix
                // projs = [row i is V * M**i, i=0..d-1] where V = input vector
                mul(buf, projs, mat);
                projs.SetDims(2*d, degdet);
                for (long i = d; i < 2*d; ++i)
                    projs[i].swap(buf[i-d]);
                // now proj = [row i is V * M**i, i=0..2*d-1]
                if (2*d <= degdet)
                    sqr(mat, mat);
                // now mat = power(M , 2*d)
            }
            t += GetWallTime()-tt;
            ++nb_iter;
        }
        timings.push_back(t/nb_iter);
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

#if defined PROFILE_KERNEL_STEP1
// to compare variants for kernel step1
// note that the "naive" kernel variant seems faster when very thin matrix
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long m = atoi(argv[1]);
    const long n = atoi(argv[2]);
    zz_p::init(NTL::GenPrime_long(atoi(argv[3])));

    Mat<zz_p> cmat,cmat2,cmat3;
    random(cmat, m, n);
    cmat2 = cmat;
    cmat3 = cmat;
    //std::cout << cmat << std::endl;

    Mat<zz_pX> pmat;
    pmat.SetDims(m,n);
    for (long i = 0; i < n; ++i)
        SetX(pmat[m-n+i][i]);
    pmat = pmat + cmat;
    //std::cout << pmat << std::endl;

    double t = GetWallTime();
    //Mat<zz_pX> kerbas;
    //VecLong shift(m);
    //pmbasis(kerbas, pmat, 2, shift);
    //t = GetWallTime() - t;
    //std::cout << "Via approx:\t" << t << std::endl;

    t = GetWallTime();
    Mat<zz_p> ckertop,ckerbot;
    kernel_step1(ckertop, ckerbot, cmat);
    t = GetWallTime() - t;
    std::cout << "Smart Gauss:\t" << t << std::endl;
    //std::cout << ckertop << std::endl<< std::endl << std::endl;
    //std::cout << ckerbot << std::endl<< std::endl<< std::endl;

    t = GetWallTime();
    Mat<zz_p> ckertop2,ckerbot2;
    kernel_step1_direct(ckertop2, ckerbot2, cmat2);
    t = GetWallTime() - t;
    std::cout << "Smart kernel:\t" << t << std::endl;
    //std::cout << ckertop2  << std::endl<< std::endl << std::endl;
    //std::cout << ckerbot2   << std::endl<< std::endl<< std::endl;

    return 0;
}

#elif defined PROFILE_KERNEL_STEP2
// to compare variants for kernel step2
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long m = atoi(argv[1]);
    const long n = atoi(argv[2]);
    zz_p::init(NTL::GenPrime_long(atoi(argv[3])));

    Mat<zz_p> cmat0,cmat1,cmat;
    random(cmat0, m, n);
    random(cmat1, m, n);
    //std::cout << cmat << std::endl;

    Mat<zz_pX> pmat;
    pmat.SetDims(m,n);
    SetCoeff(pmat, 0, cmat0);
    SetCoeff(pmat, 1, cmat1);
    for (long i = 0; i < n; ++i)
        SetCoeff(pmat[m-n+i][i],2);
    //std::cout << pmat << std::endl;

    double t = GetWallTime();
    Mat<zz_pX> kerbas;
    VecLong shift(m);
    pmbasis(kerbas, pmat, 4, shift);
    t = GetWallTime() - t;
    std::cout << "Via approx:\t" << t << std::endl;

    cmat.SetDims(m,2*n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
        {
            cmat[i][j] = cmat0[i][j];
            cmat[i][n+j] = cmat1[i][j];
        }
    t = GetWallTime();
    Mat<zz_p> ckertop,ckerbot;
    kernel_step2(ckertop, ckerbot, cmat);
    t = GetWallTime() - t;
    std::cout << "Smart kernel:\t" << t << std::endl;
    //std::cout << ckertop << std::endl<< std::endl << std::endl;
    //std::cout << ckerbot << std::endl<< std::endl<< std::endl;

    return 0;
}

#elif defined PROFILE_KERNEL_DEGREE1
// to compare variants for kernel degree1
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long m = atoi(argv[1]);
    const long n = atoi(argv[2]);
    const long d = atoi(argv[3]);
    long p = NTL::GenPrime_long(atoi(argv[4]));
    zz_p::init(p);
#if defined DEBUGGING_NOW
    std::cout << "PRIME " << p << std::endl;
#endif // DEBUGGING_NOW

    Vec<Mat<zz_p>> cmats;
    cmats.SetLength(d);
    for (long i = 0; i < d; ++i)
        random(cmats[i],m,n);
    // // concrete example in degree 2:
    // // here: long p = 594397;
    //for (long i = 0; i < 2; ++i)
    //    cmats[i].SetDims(m,n);
    //cmats[0][0][0] = 540635 ;
    //cmats[0][1][0] = 393278 ;
    //cmats[0][2][0] = 25615  ;
    //cmats[0][3][0] = 361430 ;
    //cmats[0][4][0] = 426228 ;
    //cmats[0][5][0] = 533219 ;
    //cmats[0][6][0] = 12266  ;
    //cmats[0][7][0] = 182116 ;
    //cmats[0][8][0] = 355140 ;
    //cmats[0][9][0] = 506910 ;
    //cmats[0][0][1] = 185850;
    //cmats[0][1][1] = 319804;
    //cmats[0][2][1] = 456968;
    //cmats[0][3][1] = 538120;
    //cmats[0][4][1] = 471472;
    //cmats[0][5][1] = 353017;
    //cmats[0][6][1] = 585546;
    //cmats[0][7][1] = 499752;
    //cmats[0][8][1] = 482345;
    //cmats[0][9][1] = 57564 ;
    //cmats[0][0][2] = 184671;
    //cmats[0][1][2] = 83659 ;
    //cmats[0][2][2] = 12273 ;
    //cmats[0][3][2] = 466317;
    //cmats[0][4][2] = 233868;
    //cmats[0][5][2] = 247343;
    //cmats[0][6][2] = 403402;
    //cmats[0][7][2] = 52394 ;
    //cmats[0][8][2] = 552391;
    //cmats[0][9][2] = 7779  ;
    //cmats[1][0][0] = 353628 ;
    //cmats[1][1][0] = 410695 ;
    //cmats[1][2][0] = 334048 ;
    //cmats[1][3][0] = 251812 ;
    //cmats[1][4][0] = 461843 ;
    //cmats[1][5][0] = 248587 ;
    //cmats[1][6][0] = 125230 ;
    //cmats[1][7][0] = 120852 ;
    //cmats[1][8][0] = 363768 ;
    //cmats[1][9][0] = 579117 ;
    //cmats[1][0][1] = 178165;
    //cmats[1][1][1] = 500174;
    //cmats[1][2][1] = 342199;
    //cmats[1][3][1] = 416403;
    //cmats[1][4][1] = 437184;
    //cmats[1][5][1] = 301989;
    //cmats[1][6][1] = 547901;
    //cmats[1][7][1] = 313877;
    //cmats[1][8][1] = 173123;
    //cmats[1][9][1] = 305000;
    //cmats[1][0][2] = 22627;
    //cmats[1][1][2] = 106362;
    //cmats[1][2][2] = 126392;
    //cmats[1][3][2] = 493266;
    //cmats[1][4][2] = 541673;
    //cmats[1][5][2] = 273329;
    //cmats[1][6][2] = 5298;
    //cmats[1][7][2] = 96693;
    //cmats[1][8][2] = 116158;
    //cmats[1][9][2] = 389773;
 
#ifdef DEBUGGING_NOW
    std::cout << cmats << std::endl;
#endif // DEBUGGING_NOW

    Mat<zz_pX> pmat;
    pmat.SetDims(m,n);
    for (long k = 0; k < d; ++k)
        SetCoeff(pmat, k, cmats[k]);
    for (long i = 0; i < n; ++i)
        SetCoeff(pmat[m-n+i][i],3);
    //std::cout << pmat << std::endl;

    double t = GetWallTime();
    Mat<zz_pX> kerbas;
    VecLong shift(m);
    pmbasis(kerbas, pmat, 2*d, shift);
    t = GetWallTime() - t;
    std::cout << "Via approx:\t" << t << std::endl;

    // general degree d smart kernel
    Mat<zz_p> cmat;
    cmat.SetDims(m,d*n);
    for (long i = 0; i < m; ++i)
        for (long k = 0; k < d; ++k)
            for (long j = 0; j < n; ++j)
                cmat[i][k*n+j] = cmats[k][i][j];

    t = GetWallTime();
    Mat<zz_p> ckertop,ckerbot;
    kernel_degree1(ckertop, ckerbot, cmat, d);
    t = GetWallTime() - t;
    std::cout << "Smart kernel:\t" << t << std::endl;
#ifdef DEBUGGING_NOW
    std::cout << ckertop << std::endl<< std::endl << std::endl;
    std::cout << ckerbot << std::endl<< std::endl<< std::endl;
    std::cout << cmat << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW

    // degree 1 smart kernel
    if (d>=1)
    {
        Mat<zz_p> cmat1 = cmats[0];
        t = GetWallTime();
        Mat<zz_p> ckertop1,ckerbot1;
        kernel_step1_direct(ckertop1, ckerbot1, cmat1);
        t = GetWallTime() - t;
        std::cout << "Smart kernel(1):\t" << t << std::endl;
#ifdef DEBUGGING_NOW
        std::cout << ckertop1 << std::endl<< std::endl << std::endl;
        std::cout << ckerbot1 << std::endl<< std::endl<< std::endl;
        std::cout << cmat1 << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW
    }

    // degree 2 smart kernel
    if (d>=2)
    {
        Mat<zz_p> cmat2;
        cmat2.SetDims(m,2*n);
        for (long i = 0; i < m; ++i)
            for (long k = 0; k < 2; ++k)
                for (long j = 0; j < n; ++j)
                    cmat2[i][k*n+j] = cmats[k][i][j];
        t = GetWallTime();
        Mat<zz_p> ckertop2,ckerbot2;
        kernel_step2(ckertop2, ckerbot2, cmat2);
        t = GetWallTime() - t;
        std::cout << "Smart kernel(2):\t" << t << std::endl;
#ifdef DEBUGGING_NOW
        std::cout << ckertop2 << std::endl<< std::endl << std::endl;
        std::cout << ckerbot2 << std::endl<< std::endl<< std::endl;
        std::cout << cmat2 << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW
    }

    return 0;
}

#elif defined PROFILE_KERNEL_DEGREE1_CDEG
// to compare variants for kernel degree1 cdeg
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long m = atoi(argv[1]);
    const long n = atoi(argv[2]);
    const long d = atoi(argv[3]);
    long p = NTL::GenPrime_long(atoi(argv[4]));
    zz_p::init(p);
    std::cout << "PRIME " << p << std::endl;

    Vec<Mat<zz_p>> cmats;
    cmats.SetLength(d);
    for (long i = 0; i < d; ++i)
        random(cmats[i],m,n);
    // // concrete example in degree 2:
    // // here: long p = 594397;
    //for (long i = 0; i < 2; ++i)
    //    cmats[i].SetDims(m,n);
    //cmats[0][0][0] = 540635 ;
    //cmats[0][1][0] = 393278 ;
    //cmats[0][2][0] = 25615  ;
    //cmats[0][3][0] = 361430 ;
    //cmats[0][4][0] = 426228 ;
    //cmats[0][5][0] = 533219 ;
    //cmats[0][6][0] = 12266  ;
    //cmats[0][7][0] = 182116 ;
    //cmats[0][8][0] = 355140 ;
    //cmats[0][9][0] = 506910 ;
    //cmats[0][0][1] = 185850;
    //cmats[0][1][1] = 319804;
    //cmats[0][2][1] = 456968;
    //cmats[0][3][1] = 538120;
    //cmats[0][4][1] = 471472;
    //cmats[0][5][1] = 353017;
    //cmats[0][6][1] = 585546;
    //cmats[0][7][1] = 499752;
    //cmats[0][8][1] = 482345;
    //cmats[0][9][1] = 57564 ;
    //cmats[0][0][2] = 184671;
    //cmats[0][1][2] = 83659 ;
    //cmats[0][2][2] = 12273 ;
    //cmats[0][3][2] = 466317;
    //cmats[0][4][2] = 233868;
    //cmats[0][5][2] = 247343;
    //cmats[0][6][2] = 403402;
    //cmats[0][7][2] = 52394 ;
    //cmats[0][8][2] = 552391;
    //cmats[0][9][2] = 7779  ;
    //cmats[1][0][0] = 353628 ;
    //cmats[1][1][0] = 410695 ;
    //cmats[1][2][0] = 334048 ;
    //cmats[1][3][0] = 251812 ;
    //cmats[1][4][0] = 461843 ;
    //cmats[1][5][0] = 248587 ;
    //cmats[1][6][0] = 125230 ;
    //cmats[1][7][0] = 120852 ;
    //cmats[1][8][0] = 363768 ;
    //cmats[1][9][0] = 579117 ;
    //cmats[1][0][1] = 178165;
    //cmats[1][1][1] = 500174;
    //cmats[1][2][1] = 342199;
    //cmats[1][3][1] = 416403;
    //cmats[1][4][1] = 437184;
    //cmats[1][5][1] = 301989;
    //cmats[1][6][1] = 547901;
    //cmats[1][7][1] = 313877;
    //cmats[1][8][1] = 173123;
    //cmats[1][9][1] = 305000;
    //cmats[1][0][2] = 22627;
    //cmats[1][1][2] = 106362;
    //cmats[1][2][2] = 126392;
    //cmats[1][3][2] = 493266;
    //cmats[1][4][2] = 541673;
    //cmats[1][5][2] = 273329;
    //cmats[1][6][2] = 5298;
    //cmats[1][7][2] = 96693;
    //cmats[1][8][2] = 116158;
    //cmats[1][9][2] = 389773;
 
#ifdef DEBUGGING_NOW
    std::cout << cmats << std::endl;
#endif // DEBUGGING_NOW

    Mat<zz_pX> pmat;
    pmat.SetDims(m,n);
    for (long k = 0; k < d; ++k)
        SetCoeff(pmat, k, cmats[k]);
    for (long i = 0; i < n; ++i)
        SetCoeff(pmat[m-n+i][i],3);
    //std::cout << pmat << std::endl;

    double t = GetWallTime();
    //Mat<zz_pX> kerbas;
    //VecLong shift(m);
    //pmbasis(kerbas, pmat, 2*d, shift);
    //t = GetWallTime() - t;
    //std::cout << "Via approx:\t" << t << std::endl;

    // general degree d smart kernel
    Mat<zz_p> cmat;
    cmat.SetDims(m,d*n);
    for (long i = 0; i < m; ++i)
        for (long k = 0; k < d; ++k)
            for (long j = 0; j < n; ++j)
                cmat[i][k*n+j] = cmats[k][i][j];

    VecLong cdeg(n,d);
    t = GetWallTime();
    Mat<zz_p> K0, K1, lF0;
    kernel_degree1_cdeg(K0,K1,lF0,cmats,cdeg);
    t = GetWallTime() - t;
    std::cout << "Smart kernel, general:\t" << t << std::endl;
#ifdef DEBUGGING_NOW
    std::cout << K0 << std::endl<< std::endl << std::endl;
    std::cout << K1 << std::endl<< std::endl<< std::endl;
    std::cout << lF0 << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW

    t = GetWallTime();
    Mat<zz_p> ckertop,ckerbot;
    kernel_degree1(ckertop, ckerbot, cmat, d);
    t = GetWallTime() - t;
    std::cout << "Smart kernel:\t" << t << std::endl;
#ifdef DEBUGGING_NOW
    std::cout << ckertop << std::endl<< std::endl << std::endl;
    std::cout << ckerbot << std::endl<< std::endl<< std::endl;
    std::cout << cmat << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW

    // degree 1 smart kernel
    if (d>=1)
    {
        Mat<zz_p> cmat1 = cmats[0];
        t = GetWallTime();
        Mat<zz_p> ckertop1,ckerbot1;
        kernel_step1_direct(ckertop1, ckerbot1, cmat1);
        t = GetWallTime() - t;
        std::cout << "Smart kernel(1):\t" << t << std::endl;
#ifdef DEBUGGING_NOW
        std::cout << ckertop1 << std::endl<< std::endl << std::endl;
        std::cout << ckerbot1 << std::endl<< std::endl<< std::endl;
        std::cout << cmat1 << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW
    }

    // degree 2 smart kernel
    if (d>=2)
    {
        Mat<zz_p> cmat2;
        cmat2.SetDims(m,2*n);
        for (long i = 0; i < m; ++i)
            for (long k = 0; k < 2; ++k)
                for (long j = 0; j < n; ++j)
                    cmat2[i][k*n+j] = cmats[k][i][j];
        t = GetWallTime();
        Mat<zz_p> ckertop2,ckerbot2;
        kernel_step2(ckertop2, ckerbot2, cmat2);
        t = GetWallTime() - t;
        std::cout << "Smart kernel(2):\t" << t << std::endl;
#ifdef DEBUGGING_NOW
        std::cout << ckertop2 << std::endl<< std::endl << std::endl;
        std::cout << ckerbot2 << std::endl<< std::endl<< std::endl;
        std::cout << cmat2 << std::endl<< std::endl<< std::endl;
#endif // DEBUGGING_NOW
    }

    return 0;
}

#elif defined TIME_FMA
// to compare variants for kernel degree1 cdeg
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long m = atoi(argv[1]);
    const long D = atoi(argv[2]);
    const long k = atoi(argv[3]);
    long p = NTL::GenPrime_long(atoi(argv[4]));
    zz_p::init(p);
    std::cout << "PRIME " << p << std::endl;
    std::cout << "Time FMA" << std::endl;

    Mat<zz_p> fmatp, matp, mat;
    random(mat, k+5, m);
    random(matp, m, D);

    double t=GetWallTime();
    fused_mul_add(fmatp, matp, m, mat, k);
    std::cout << GetWallTime() - t << std::endl;
}


#elif defined UNIFORM_KSHIFTED
// to check time for charpoly of d-shifted
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    const long m = atoi(argv[1]);
    const long d = atoi(argv[2]);
    long p = NTL::GenPrime_long(atoi(argv[3]));
    zz_p::init(p);
    std::cout << "PRIME " << p << std::endl;

    VecLong cdeg(m,d);
    Vec<Mat<zz_p>> matp;
    for (long thres = 10; thres < 51; thres+=10)
    {
        matp.kill();
        matp.SetLength(cdeg[0]);
        for (long k = 0; k < cdeg[0]; ++k)
            random(matp[k], m, m);
        Mat<zz_pX> pmat;
        conv(pmat,matp);
        for (long i = 0; i < m; ++i)
            SetCoeff(pmat[i][i], d);
        double t = GetWallTime();
        VecLong copy_cdeg = cdeg;
        zz_pX det;
        bool ok = determinant_shifted_form_smartkernel_updateall(det, matp, copy_cdeg, thres, m*d);
        t = GetWallTime()-t;
        ok = ok && verify_determinant(det, pmat, true, true);
        if (not ok)
            std::cout << "~~~Warning~~~ verification of determinant failed in CHARPOLY" << std::endl;
        std::cout << "Time: " << t << std::endl;
    }
}

#else
int main(int argc, char ** argv)
{
    NTL::SetSeed(ZZ(2));
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    std::vector<string> labels;
#ifdef TIME_LINSOLVE
    labels.push_back("rand-linsolve");
#endif
#ifdef TIME_TRI_DECR
    labels.push_back("triangular-decr");
#endif
#ifdef TIME_TRI_INCR
    labels.push_back("triangular-incr");
#endif
#ifdef TIME_DEG_UPMAT
    labels.push_back("deg-upmat1");
    labels.push_back("deg-upmat2");
    labels.push_back("deg-upmat3");
    labels.push_back("deg-upmat4");
    labels.push_back("deg-upmat5");
    labels.push_back("deg-upmat6");
    labels.push_back("deg-upmat7");
#endif
#ifdef TIME_DEG_UPKER
    labels.push_back("deg-upker1");
    labels.push_back("deg-upker2");
    labels.push_back("deg-upker3");
    labels.push_back("deg-upker4");
    labels.push_back("deg-upker5");
    labels.push_back("deg-upker6");
    labels.push_back("deg-upker7");
#endif
#ifdef TIME_SMART_UPMAT
    labels.push_back("ker-umat20");
    labels.push_back("ker-umat25");
    labels.push_back("ker-umat30");
    labels.push_back("ker-umat35");
    labels.push_back("ker-umat40");
    labels.push_back("ker-umat45");
#endif
#ifdef TIME_ZLS_AVG
    labels.push_back("zls-avgdeg");
#endif
#ifdef ESTIMATE_WIEDEMANN
    labels.push_back("est-Wied.");
#endif
#ifdef ESTIMATE_KELLERGEHRIG
    labels.push_back("est-K.-G.");
#endif

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
        cout << "nbits\tdim\tdegdet\tdensity\t";
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
#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
