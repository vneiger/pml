// timing charpoly of matrix in "shifted forms" [Pernet-Storjohann]
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

#define PROFILE_DEG_UPMAT
#define PROFILE_DEG_UPKER
//#define PROFILE_ZLS_AVG
//#define PROFILE_KER_AVG
//#define PROFILE_KERNEL_STEP1
//#define DEBUGGING_NOW

//#define TIME_LINSOLVE // via linear system solving with random rhs
//#define TIME_TRI_DECR // generic determinant on matrix with decreasing diagonal degrees
#define TIME_TRI_INCR // generic determinant on matrix with increasing diagonal degrees
#define TIME_DEG_UPMAT // going degree by degree, updating remaining matrix columns at each iteration
#define TIME_DEG_UPKER // going degree by degree, updating kernel basis at each iteration
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
//   [ --- ]  =  - [ --- ] * iF_0
//   [ K_1 ]       [ F_2 ]
// where iF_0 = inverse(F_0)
// Consequence 2: one reduced left kernel basis of F is
//   [  K_1 - x iF_0  |    0    |  I_n ]
//   [  -------------   -------   ---- ]
//   [       K_0      |  I_ell  |   0  ]
// and the Popov left kernel basis of F is
//   [  F_0 K_1 + x I_n  |    0    |  F_0 ]
//   [  ---------------    -------   ---- ]
//   [        K_0        |  I_ell  |   0  ]
// Input : (ell+2n) x n constant matrix
//        [  F_0  ]
//        [  ---- ]
//  F =   [  F_1  ]
//        [  ---- ]
//        [  F_2  ]
// Output : 
//   - kertop: n x n constant matrix kertop = F_0 K_1 
//   - kerbot: ell x n constant matrix kerbot = K_0
//   - cmat: replaced by only its n top rows F_0
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
    Mat<zz_p> imattop;
    kertop.SetDims(n,n);
    for (long i = 0; i < n; ++i)
        kertop[i].swap(cmat[ell+n+i]);
    kerbot.SetDims(ell,n);
    for (long i = 0; i < ell; ++i)
        kerbot[i].swap(cmat[n+i]);

    cmat.SetDims(n,n);
    inv(imattop,cmat);
    mul(kertop, kertop, imattop);
    NTL::negate(kertop, kertop);
    mul(kerbot, kerbot, imattop);
    NTL::negate(kerbot, kerbot);
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
// kernel basis of a degree-1  (ell+2n) x n matrix of the form
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
//   [ ----- ----- ]
//   [ F_10 | F_11 ]
//   [ ----- ----- ]
//   [ F_20 | F_21 ]
// has the form
//   [  K_0  |  I_ell  |  0  ]
//   [  ----   -------   --- ]
//   [  K_1  |    0    | I_n ]
// where
//   [ K_0 ]       [ F_10 | F_11 ]
//   [ --- ]  =  - [ ----   ---- ] * iF_0
//   [ K_1 ]       [ F_20 | F_11 ]
// where iF_0 = inverse([F_00 | F_01])
// Note: all rows in the top part (involving K_0) are
//       obviously also in the left kernel of F; in fact
//       they form the constant part of Ker(F)
// Consequence 2: one reduced left kernel basis of F is
//   [  K_1 - x iF_0  |    0    |  I_n ]
//   [  -------------   -------   ---- ]
//   [       K_0      |  I_ell  |   0  ]
// and the Popov left kernel basis of F is
//   [  F_0 K_1 + x I_n  |    0    |  F_0 ]
//   [  ---------------    -------   ---- ]
//   [        K_0        |  I_ell  |   0  ]
// Input : (ell+2n) x n constant matrix
//        [  F_0  ]
//        [  ---- ]
//  F =   [  F_1  ]
//        [  ---- ]
//        [  F_2  ]
// Output : 
//   - kertop: n x n constant matrix kertop = F_0 K_1 
//   - kerbot: ell x n constant matrix kerbot = K_0
//   - cmat: replaced by only its n top rows F_0
void kernel_step2(Mat<zz_p> & kertop, Mat<zz_p> & kerbot, Mat<zz_p> & cmat)
{
    return;
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
    for (long thres=1; thres < 8; ++thres)
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
            std::cout << "~~~Warning~~~ verification of determinant failed in DEG_UPKER approach" << std::endl;
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

#ifdef PROFILE_KERNEL_STEP1
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
#else
int main(int argc, char ** argv)
{
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
