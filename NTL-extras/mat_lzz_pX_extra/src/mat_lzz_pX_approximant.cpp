#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

//#define MBASIS1_PROFILE // FIXME
//#define PMBASIS_PROFILE // FIXME

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ALGORITHMS FOR APPROXIMANT BASES                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* VERIFICATION                                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* TODO doc currently in .h --> should be moved here??        */
/*------------------------------------------------------------*/
// follows ideas from Algorithm 1 in Giorgi-Neiger, ISSAC 2018
bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const Order & order,
                          const Shift & shift,
                          const PolMatForm & form,
                          const bool row_wise,
                          const bool randomized
                         )
{
    // TODO far from optimal except in the balanced case
    // (e.g. could be improved when deg(appbas)<<deg(pmat) like in Hermite-Pade,
    // or when appbas has strange column degrees or strange row degrees)
    if (randomized)
        throw std::logic_error("==is_approximant_basis== Fast randomized approximant basis verification not implemented yet");
    else
        std::cout << "==is_approximant_basis== Warning: using *randomized* algorithm for testing generation, in the verification that det(appbas) = c X^d" << std::endl;

    // test whether appbas has the right dimensions
    if (appbas.NumRows() != appbas.NumCols()
        || (row_wise && appbas.NumCols() != pmat.NumRows())
        || ((not row_wise) && appbas.NumRows() != pmat.NumCols()))
        return false;

    // test whether appbas is shift-reduced with form at least 'form'
    if (not is_polmatform(appbas,form,shift,row_wise))
        return false;

    // test whether appbas consists of approximants (if row-wise: appbas * pmat = 0 mod X^order)
    // and retrieve the constant coefficient "cmat" of degree order (if row-wise: cmat = (appbas * pmat * X^{-order})  mod X)
    // (we reserve some additional space in cmat because later it will store the constant coefficient of appbas)
    Mat<zz_pX> residual;
    if (row_wise)
        multiply_evaluate(residual,appbas,pmat);
    else
        multiply_evaluate(residual,pmat,appbas);
    // TODO this multiplication could be:
    //   - truncated mod X^{order+1}
    //   - improved by taking degree profile into account

    Mat<zz_p> cmat;
    if (row_wise)
        cmat.SetDims(residual.NumRows(),residual.NumCols()+appbas.NumRows());
    else
        cmat.SetDims(residual.NumRows()+appbas.NumRows(),residual.NumCols());

    for (long i = 0; i < residual.NumRows(); ++i)
    {
        for (long j = 0; j < residual.NumCols(); ++j)
        {
            long ord = row_wise ? order[j] : order[i];
            GetCoeff(cmat[i][j],residual[i][j],ord);
            trunc(residual[i][j],residual[i][j],ord);
            if (residual[i][j] != 0)
                return false;
        }
    }

    // test whether det(appbas) is a monomial (power of x):
    // det(appbas) = det(appbas(1)) * X^(deg(det(appbas)))
    // the degree of det(appbas) is known since appbas is reduced:
    DegVec degs = vector_degree(appbas,shift,row_wise);
    long sum_degs = std::accumulate(degs.begin(), degs.end(), (long)0);
    long sum_shift = std::accumulate(shift.begin(), shift.end(), (long)0);
    long degdet = sum_degs - sum_shift;
    // for a random point pt, verify that det(appbas(pt)) = det(appbas(1)) * pt^degdet
    zz_p pt = random_zz_p();
    zz_p det_pt = determinant(eval(appbas, pt));
    zz_p det_one = determinant(eval(appbas, zz_p(1)));
    if (det_pt != det_one * power(pt,degdet))
        return false;

    // generation test: verify that [ cmat  P(0) ] has full rank (see Giorgi-Neiger ISSAC 2018)
    if (row_wise)
        for (long i = 0; i < appbas.NumRows(); ++i)
            for (long j = 0; j < appbas.NumCols(); ++j)
                cmat[i][j+residual.NumCols()] = coeff(appbas[i][j],0);
    else
        for (long i = 0; i < appbas.NumRows(); ++i)
            for (long j = 0; j < appbas.NumCols(); ++j)
                cmat[i+residual.NumRows()][j] = coeff(appbas[i][j],0);
    long rank = gauss(cmat);
    if (rank != std::min(cmat.NumRows(),cmat.NumCols()))
        return false;

    return true;
}

bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order,
                          const Shift & shift,
                          const PolMatForm & form,
                          const bool row_wise,
                          const bool randomized
                         )
{
    Order orders(pmat.NumRows(),order);
    return is_approximant_basis(appbas,pmat,orders,shift,form,row_wise,randomized);
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ITERATIVE ALGORITHMS                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* for the general case                                       */
/*------------------------------------------------------------*/
DegVec appbas_iterative(
                        Mat<zz_pX> & appbas,
                        const Mat<zz_pX> & pmat,
                        const Order & order,
                        const Shift & shift,
                        bool order_wise
                       )
{
    /** Two possibilities (among others) for next coefficient to deal with:
     *   - process 'pmat' order-wise (choose column with largest order)
     *   - process 'pmat' column-wise (choose leftmost column not yet completed)
     **/

    long rdim = pmat.NumRows();
    long cdim = pmat.NumCols();

    // initial approximant basis: identity of dimensions 'rdim x rdim'
    set(appbas,rdim);

    // initial residual: the whole input matrix
    Mat<zz_pX> residual( pmat );

    // order that remains to be dealt with
    Order rem_order( order );

    // indices of columns/orders that remain to be dealt with
    std::vector<long> rem_index( cdim );
    std::iota(rem_index.begin(), rem_index.end(), 0);

    // shifted row degrees of approximant basis
    // (initially, of the identity matrix, i.e. rdeg == shift)
    DegVec rdeg(shift);

    while (not rem_order.empty())
    {
        /** Invariant:
         *  - appbas is a shift-ordered weak Popov approximant basis for
         *  (pmat,reached_order) where doneorder is the tuple such that
         *  -->reached_order[j] + rem_order[j] == order[j] for j appearing in rem_index
         *  -->reached_order[j] == order[j] for j not appearing in rem_index
         *  - rdeg == the shift-row degree of appbas
         *  - residual == submatrix of columns (appbas * pmat)[:,j] for all j such that reached_order[j] < order[j]
         **/

        long j=0; // value if columnwise (order_wise==False)
        if (order_wise)
        {
            // FIXME if seems to slow (e.g. compared to non-order-wise), could be
            // valuable to simply initially permute the columns of pmat and order:
            // then 'j' is obvious
            j = std::distance(rem_order.begin(), std::max_element(rem_order.begin(), rem_order.end()));
        }

        long deg = order[rem_index[j]] - rem_order[j];

        // record the coefficients of degree deg of the column j of residual
        // also keep track of which of these are nonzero,
        // and among the nonzero ones, which is the first with smallest shift
        Vec<zz_p> const_residual;
        const_residual.SetLength(rdim);
        std::vector<long> indices_nonzero;
        long piv = -1;
        for (long i = 0; i < rdim; ++i)
        {
            const_residual[i] = coeff(residual[i][j],deg);
            if (const_residual[i] != 0)
            {
                indices_nonzero.push_back(i);
                if (piv<0 || rdeg[i] < rdeg[piv])
                    piv = i;
            }
        }

        // if indices_nonzero is empty, const_residual is already zero, there is nothing to do
        if (not indices_nonzero.empty())
        {
            // update all rows of appbas and residual in indices_nonzero except row piv
            for (long row : indices_nonzero)
            {
                if (row!=piv)
                {
                    zz_p c = - const_residual[row] / const_residual[piv];
                    for (long k=0; k<rdim; ++k)
                        appbas[row][k] += c * appbas[piv][k];
                    for (long k=0; k<residual.NumCols(); ++k)
                        residual[row][k] += c * residual[piv][k];
                }
            }

            // update row piv
            ++rdeg[piv]; // shifted row degree of row piv increases
            for (long k=0; k<rdim; ++k) // TODO use shiftRow
                appbas[piv][k] <<= 1; // row piv multiplied by X
            for (long k=0; k<residual.NumCols(); ++k)
            {
                residual[piv][k] <<= 1; // row piv multiplied by X
                trunc(residual[piv][k],residual[piv][k],order[rem_index[k]]); // truncate
            }
        }

        // now column j (or rather rem_index[j] in original 'pmat') is zero mod X^(deg+1)
        if (rem_order[j] > 1) // one step of progress towards order[j]
            --rem_order[j];
        else // rem_order[j] == 1: work completed for column j
        {
            rem_order.erase(rem_order.begin() + j);
            rem_index.erase(rem_index.begin() + j);
            if (!rem_order.empty())
            {
                Mat<zz_pX> buffer( residual );
                residual.kill();
                residual.SetDims(rdim,rem_order.size());
                for (long i=0; i<rdim; ++i)
                {
                    for (long k=0; k<j; ++k)
                        residual[i][k] = buffer[i][k];
                    for (long k=j+1; k<buffer.NumCols(); ++k)
                        residual[i][k-1] = buffer[i][k];
                }
            }
        }
    }

    // make rdeg contain the pivot degree rather than shifted row degree, i.e.:
    // rdeg = rdeg-shift   , done entry-wise
    std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::minus<long>());
    return rdeg;
}

/*------------------------------------------------------------*/
/* general case, output in Popov form                         */
/*------------------------------------------------------------*/
DegVec popov_appbas_iterative(
                              Mat<zz_pX> & appbas,
                              const Mat<zz_pX> & pmat,
                              const Order & order,
                              const Shift & shift,
                              bool order_wise
                             )
{
    // TODO: first call can be very slow if strange degrees --> rather implement
    // BecLab00's "continuous" normalization?
    DegVec pivdeg = appbas_iterative(appbas,pmat,order,shift,order_wise);
    Shift new_shift( pivdeg );
    std::transform(new_shift.begin(), new_shift.end(), new_shift.begin(), std::negate<long>());
    clear(appbas);
    appbas_iterative(appbas,pmat,order,new_shift,order_wise);
    Mat<zz_p> lmat;
    leading_matrix(lmat, appbas, new_shift, true);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas);
    return pivdeg;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ALGORITHMS FOR UNIFORM ORDER                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// These algorithms are defined for arbitrary shifts, but
// work best for shifts close to uniform
// (except the next one, popov_mbasis1, where the shift has roughly no
// influence)

DegVec popov_mbasis1(
                     Mat<zz_p> & kerbas,
                     const Mat<zz_p> & pmat,
                     const Shift & shift
                    )
{
    // compute permutation which realizes stable sort of the shift
    // (i.e. sorts (shift[0],0)....(shift[len],len) lexicographically increasingly)
#ifdef MBASIS1_PROFILE
    double t_perm1,t_perm2,t_pivind,t_ker,t_now;
    t_now = GetWallTime();
#endif
    std::vector<long> perm_shift(pmat.NumRows());
    std::iota(perm_shift.begin(), perm_shift.end(), 0);
    stable_sort(perm_shift.begin(), perm_shift.end(),
                [&](const long& a, const long& b)->bool
                {
                return (shift[a] < shift[b]);
                } );

    // permute rows of pmat accordingly
    Mat<zz_p> mat;
    mat.SetDims(pmat.NumRows(),pmat.NumCols());
    for (long i = 0; i < pmat.NumRows(); ++i)
        //for (long j = 0; j < mat.NumRows(); ++j)
        mat[i] = pmat[perm_shift[i]];
#ifdef MBASIS1_PROFILE
    t_perm1 = GetWallTime() - t_now;
#endif

    // find the permuted kernel basis in row echelon form
#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    Mat<zz_p> p_kerbas;
    kernel(p_kerbas,mat);
#ifdef MBASIS1_PROFILE
    t_ker = GetWallTime() - t_now;
#endif
    if (p_kerbas.NumRows()==0)
        return DegVec(pmat.NumRows(),1);
    if (p_kerbas.NumRows()==pmat.NumRows())
        return DegVec(pmat.NumRows(),0);

    // compute the (permuted) pivot indices
    // (NTL doesn't return the pivot indices in Gaussian elimination, we might
    // hack the NTL code to retrieve them directly but it seems that the next
    // lines have negligible time compared to the kernel computation)
#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    std::vector<long> p_pivind(p_kerbas.NumRows(),p_kerbas.NumCols()-1); // pivot indices in permuted kernel basis
    //p_pivind.back() = p_kerbas.NumCols()-1;
    for (long i = 0; i<p_kerbas.NumRows(); ++i)
    {
        while (p_pivind[i]>=0 && p_kerbas[i][p_pivind[i]]==0)
            --p_pivind[i];
    }
#ifdef MBASIS1_PROFILE
    t_pivind = GetWallTime() - t_now;
#endif

    // permute everything back to original order:
    // prepare kernel permutation by permuting kernel pivot indices
    // pivot degrees corresponding to kernel pivot indices are 0, others are 1
#ifdef MBASIS1_PROFILE
    t_now = GetWallTime();
#endif
    DegVec pivdeg(pmat.NumRows(),1);
    std::vector<long> pivind(p_kerbas.NumRows());
    for (long i = 0; i < p_kerbas.NumRows(); ++i)
    {
        pivind[i] = perm_shift[p_pivind[i]];
        pivdeg[pivind[i]] = 0;
    }

    std::vector<long> perm_rows_ker(p_kerbas.NumRows());
    std::iota(perm_rows_ker.begin(), perm_rows_ker.end(), 0);
    sort(perm_rows_ker.begin(), perm_rows_ker.end(),
         [&](const long& a, const long& b)->bool
         {
         return (pivind[a] < pivind[b]);
         } );


    kerbas.SetDims(p_kerbas.NumRows(),p_kerbas.NumCols());
    for (long i = 0; i < kerbas.NumRows(); ++i)
        for (long j = 0; j < kerbas.NumCols(); ++j)
            kerbas[i][perm_shift[j]] = p_kerbas[perm_rows_ker[i]][j];
#ifdef MBASIS1_PROFILE
    t_perm2 = GetWallTime() - t_now;
#endif

#ifdef MBASIS1_PROFILE
    std::cout << "~~popov_mbasis1~~ input dimensions " << pmat.NumRows() << " x " << pmat.NumCols() << std::endl;
    std::cout << "  Initial permutations: " << t_perm1 << std::endl;
    std::cout << "  Computing kernel: " << t_ker << std::endl;
    std::cout << "  Computing pivot indices: " << t_pivind << std::endl;
    std::cout << "  Final permutations: " << t_perm2 << std::endl;
#endif

    return pivdeg;
}

/*------------------------------------------------------------*/
/* plain mbasis with polynomial matrices                      */
/*------------------------------------------------------------*/
DegVec mbasis(
              Mat<zz_pX> & appbas,
              const Mat<zz_pX> & pmat,
              const long order,
              const Shift & shift
             )
{
    // initially, appbas is the identity matrix
    set(appbas,pmat.NumRows());

    // holds the current shifted row degree of appbas
    // initially, this is exactly shift
    DegVec rdeg(shift);

    // TODO should we keep this?
    // (is the code below really doing something if zero matrix?)
    if ( IsZero(pmat) )
        return rdeg;

    long deg_pmat = deg(pmat);

    // holds the current pivot degree of appbas
    // initially tuple of zeroes
    // (note that at all times pivdeg+shift = rdeg entrywise)
    DegVec pivdeg(pmat.NumRows());

    // will store the pivot degree at each call of mbasis1
    DegVec diff_pivdeg;

    // matrix to store the kernels in mbasis1 calls
    Mat<zz_p> kerbas;
    // matrix to store residuals, initially constant coeff of pmat
    Mat<zz_p> residual( coeff(pmat,0) );

    // declare matrices
    Mat<zz_p> res_coeff,res_coeff1,res_coeff2; // will store coefficient matrices used to compute the residual
    Mat<zz_pX> kerapp; // will store constant-kernel * appbas

    for (long ord = 1; ord <= order; ++ord)
    {
        diff_pivdeg = popov_mbasis1(kerbas,residual,rdeg);

        if (kerbas.NumRows()==0)
        {
            // computation is already finished: the final basis is X^(order-ord+1)*appbas
            appbas <<= (order-ord+1);
            // update pivdeg accordingly, and return
            std::for_each(pivdeg.begin(), pivdeg.end(), [&order,&ord](long& a) { a+=order-ord+1; });
            return pivdeg;
        }

        // kerbas.NumRows()==residual.NumRows() --> approximant basis is already
        // correct for this order, just go to the next

        if (kerbas.NumRows()<residual.NumRows())
        {
            // I/ Update degrees:
            // new shifted row degree = old rdeg + diff_pivdeg
            std::transform(rdeg.begin(), rdeg.end(), diff_pivdeg.begin(), rdeg.begin(), std::plus<long>());
            // new pivot degree = old pivot_degree + diff_pivdeg
            std::transform(pivdeg.begin(), pivdeg.end(), diff_pivdeg.begin(), pivdeg.begin(), std::plus<long>());
            // deduce degree of appbas; note that it is a property of this algorithm
            // that deg(appbas) = max(pivot degree) (i.e. max(degree of diagonal
            // entries); this does not hold in general for ordered weak Popov forms
            long deg_appbas = *std::max_element(pivdeg.begin(), pivdeg.end());

            // II/ update approximant basis

            // submatrix of rows with diff_pivdeg==0 is replaced by kerbas*appbas
            mul(kerapp,kerbas,appbas);
            long row=0;
            // TODO have function to copy into submatrix??
            for (long i = 0; i < appbas.NumRows(); ++i)
            {
                if (diff_pivdeg[i]==0)
                {
                    appbas[i] = kerapp[row];
                    ++row;
                }
            }
            // rows with diff_pivdeg=1 are simply multiplied by X
            for (long i = 0; i < appbas.NumRows(); ++i)
                if (diff_pivdeg[i]==1)
                    LeftShiftRow(appbas,appbas,i,1);

            // III/ compute next residual, if needed
            // this is coefficient of degree ord in appbas * pmat
            if (ord<order)
            {
                clear(residual);
                for (long d = std::max<long>(0,ord-deg_pmat); d <= deg_appbas; ++d) // note that deg_appbas <= ord holds
                {
                    res_coeff1 = coeff(appbas,d);
                    res_coeff2 = coeff(pmat,ord-d);
                    mul(res_coeff, res_coeff1, res_coeff2);
                    add(residual, residual, res_coeff);
                }
            }
        }
    }

    return pivdeg;
}

/*------------------------------------------------------------*/
/* mbasis, using vectors of matrices                          */
/*------------------------------------------------------------*/
// FIXME representing output as vector of matrices seems wrong if max degree is large (Hermite-Pade with unbalanced degrees) --> should be avoided in this case???)
DegVec mbasis_vector(
                     Mat<zz_pX> & appbas,
                     const Mat<zz_pX> & pmat,
                     const long order,
                     const Shift & shift
                    )
{
    Vec<Mat<zz_p>> coeffs_pmat = conv(pmat,order);
    long nrows = coeffs_pmat[0].NumRows();
    Vec<Mat<zz_p>> coeffs_appbas;

    // initially, coeffs_appbas is the identity matrix
    coeffs_appbas.SetLength(1);
    coeffs_appbas[0] = ident_mat_zz_p(nrows);

    // holds the current shifted row degree of coeffs_appbas
    // initially, this is exactly shift
    DegVec rdeg( shift );

    // holds the current pivot degree of coeffs_appbas
    // initially tuple of zeroes
    // (note that at all times pivdeg+shift = rdeg entrywise)
    DegVec pivdeg(nrows);

    // will store the pivot degree at each call of mbasis1
    DegVec diff_pivdeg;

    // matrix to store the kernels in mbasis1 calls
    Mat<zz_p> kerbas;
    // matrix to store residuals, initially constant coeff of coeffs_pmat
    Mat<zz_p> residual( coeffs_pmat[0] );

    // declare matrices
    Mat<zz_p> res_coeff,res_coeff1,res_coeff2; // will store coefficient matrices used to compute the residual
    Mat<zz_p> kerapp; // will store constant-kernel * coeffs_appbas[d]

    for (long ord = 1; ord <= order; ++ord)
    {
        // call MBasis1 to retrieve kernel and pivdeg
        diff_pivdeg = popov_mbasis1(kerbas,residual,rdeg);

        if (kerbas.NumRows()==0)
        {
            // computation is already finished: the final basis is X^(order-ord+1)*coeffs_appbas
            appbas = conv(coeffs_appbas);
            appbas <<= (order-ord+1);
            // update pivdeg accordingly, and return
            std::for_each(pivdeg.begin(), pivdeg.end(), [&order,&ord](long& a) { a+=order-ord+1; });
            return pivdeg;
        }

        // kerbas.NumRows()==residual.NumRows() --> approximant basis is already
        // correct for this order, just go to the next

        if (kerbas.NumRows()<residual.NumRows())
        {
            // I/ Update degrees:
            // new shifted row degree = old rdeg + diff_pivdeg
            std::transform(rdeg.begin(), rdeg.end(), diff_pivdeg.begin(), rdeg.begin(), std::plus<long>());
            // new pivot degree = old pivot_degree + diff_pivdeg
            std::transform(pivdeg.begin(), pivdeg.end(), diff_pivdeg.begin(), pivdeg.begin(), std::plus<long>());
            // deduce degree of coeffs_appbas; note that it is a property of this algorithm
            // that deg(coeffs_appbas) = max(pivot degree) (i.e. max(degree of diagonal
            // entries); this does not hold in general for ordered weak Popov forms
            long deg_appbas = *std::max_element(pivdeg.begin(), pivdeg.end());
            coeffs_appbas.SetLength(deg_appbas+1);
            coeffs_appbas[deg_appbas].SetDims(nrows, nrows);

            // II/ update approximant basis

            // submatrix of rows with diff_pivdeg==0 is replaced by kerbas*coeffs_appbas
            // while rows with diff_pivdeg=1 are simply multiplied by X
            // --> the loop goes downwards, so that we can do both in the same iteration
            for (long d = deg_appbas-1; d >= 0; --d)
            {
                kerapp = kerbas * coeffs_appbas[d];
                long row=0;
                for (long i = 0; i < nrows; ++i)
                {
                    if (diff_pivdeg[i]==0)
                    {
                        coeffs_appbas[d][i] = kerapp[row];
                        ++row;
                    }
                    else  // diff_pivdeg[i]==1 --> multiply by X
                    {
                        coeffs_appbas[d+1][i] = coeffs_appbas[d][i];
                        if (d==0) // put zero row
                            clear(coeffs_appbas[0][i]);
                    }
                }
            }

            // III/ compute next residual, if needed
            // this is coefficient of degree ord in appbas * pmat
            if (ord<order)
            {
                clear(residual);
                for (long d = std::max<long>(0,ord-coeffs_pmat.length()+1); d <= deg_appbas; ++d) // we have deg_appbas <= ord
                {
                    mul(res_coeff, coeffs_appbas[d], coeffs_pmat[ord-d]);
                    add(residual, residual, res_coeff);
                }
            }
        }
    }

    appbas = conv(coeffs_appbas);
    return pivdeg;
}

/*------------------------------------------------------------*/
/* variant with continuous update of the residual             */
/*------------------------------------------------------------*/
DegVec mbasis_resupdate(
                        Mat<zz_pX> & appbas,
                        const Mat<zz_pX> & pmat,
                        const long order,
                        const Shift & shift
                       )
{
    // initially, appbas is the identity matrix
    set(appbas,pmat.NumRows());

    // holds the current shifted row degree of appbas
    // initially, this is exactly shift
    DegVec rdeg( shift );

    // buffer for temporary pivdegs at each mbasis1 call
    DegVec pivdeg(pmat.NumRows());
    // here, unlike in mbasis(), we do not hold the current pivot degree of
    // appbas (there the main motivation for having it was that we needed
    // to maintain the actual degree of the basis)

    // matrix to store the kernels in mbasis1 calls
    Mat<zz_p> kerbas;

    // matrix to temporarily store constant coeff of residual
    Mat<zz_p> res_const;

    // matrices which will be constant-kernel*appbas and constant-kernel*residual
    Mat<zz_pX> kerapp,kerres;

    // matrix to store the residual, initially equal to pmat mod X^order
    Mat<zz_pX> residual;
    trunc(residual,pmat,order);

    for (long ord = 1; ord <= order; ++ord)
    {
        // call MBasis1 to retrieve kernel and pivdeg
        pivdeg = popov_mbasis1(kerbas,coeff(residual,0),rdeg);

        if (kerbas.NumRows()==0)
        {
            // computation is already finished: the final basis is X^(order-ord+1)*appbas
            appbas <<= (order-ord+1);
            // update rdeg accordingly, then transform it to pivdeg, and return
            std::for_each(rdeg.begin(), rdeg.end(), [&order,&ord](long& a) { a+=order-ord+1; });
            std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::minus<long>());
            return rdeg;
        }
        else if (kerbas.NumRows()==residual.NumRows())
        {
            // approximant basis is already correct for this order
            // just update the residual (discard zero constant coefficient), if needed
            if (ord<order)
                residual >>= 1;
        }
        else
        {
            // update approximant basis
            // submatrix of rows with pivdeg=0 is replaced by kerbas*appbas
            mul(kerapp,kerbas,appbas);
            long row=0;
            // TODO have function to copy into submatrix??
            for (long i = 0; i < appbas.NumRows(); ++i)
            {
                if (pivdeg[i]==0)
                {
                    appbas[i] = kerapp[row];
                    ++row;
                }
            }
            // rows with pivdeg=1 are simply multiplied by X
            for (long i = 0; i < appbas.NumRows(); ++i)
                if (pivdeg[i]==1)
                    LeftShiftRow(appbas,appbas,i,1);

            // if we are not finished, update residual so that it remains equal to
            // X^(-ord) appbas*pmat mod X^(order-ord)
            if (ord<order)
            {
                // keep the constant matrix as a temp, and shift residual = X^(-1) residual
                GetCoeff(res_const, residual, 0);
                residual >>= 1;
                // submatrix of rows with pivdeg=0 is replaced by kerbas*residual
                mul(kerres,kerbas,residual);
                row=0;
                // TODO have function to copy into submatrix??
                for (long i = 0; i < residual.NumRows(); ++i)
                {
                    if (pivdeg[i]==0)
                    {
                        residual[i] = kerres[row];
                        ++row;
                    }
                }
                // rows with pivdeg=1 are multiplied by X
                for (long i = 0; i < residual.NumRows(); ++i)
                {
                    if (pivdeg[i]==1)
                    {
                        // multiply by X and truncate mod X^(order-ord)
                        LeftShiftRow(residual,residual,i,1);
                        truncRow(residual,residual,i,order-ord);
                        for (long j = 0; j < residual.NumCols(); ++j)
                            SetCoeff(residual[i][j],0,res_const[i][j]); 
                    }
                }
            }

            // new shifted row degree = old rdeg + pivdeg  (entrywise)
            std::transform(rdeg.begin(), rdeg.end(), pivdeg.begin(), rdeg.begin(), std::plus<long>());
        }
    }

    std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::minus<long>());
    return rdeg;
}

/*------------------------------------------------------------*/
/* M-Basis returning Popov basis                              */
/*------------------------------------------------------------*/
DegVec popov_mbasis(
                    Mat<zz_pX> &appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    const Shift & shift
                   )
{
    DegVec pivdeg = mbasis(appbas,pmat,order,shift);
    Shift new_shift( pivdeg );
    std::transform(new_shift.begin(), new_shift.end(), new_shift.begin(), std::negate<long>());
    clear(appbas);
    mbasis(appbas,pmat,order,new_shift);
    Mat<zz_p> lmat;
    leading_matrix(lmat, appbas, new_shift, true);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas);
    return pivdeg;
}

/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis                                */
/*------------------------------------------------------------*/
DegVec pmbasis(
               Mat<zz_pX> & appbas,
               const Mat<zz_pX> & pmat,
               const long order,
               const Shift & shift
              )
{
#ifdef PMBASIS_PROFILE
    std::cout << "Order: " << order << std::endl;
    double t1,t2;
#endif
    // TODO thresholds to be determined:
    //  --> from mbasis (only linalg) to pmbasis with low-degree polmatmul (Karatsuba...)
    //  --> from this pmbasis to pmbasis with eval-based polmatmul (FFT, geometric..)
#ifdef PMBASIS_PROFILE
    if (order <= 32)
    {
        t1 = GetWallTime();
        DegVec rdeg = mbasis_vector(appbas,pmat,order,shift);
        t2 = GetWallTime();
        std::cout << "\tTime(base-case): " << (t2-t1) << "s" << std::endl;
        return rdeg;
    }
#else
    if (order <= 32)
        return mbasis_vector(appbas,pmat,order,shift);
#endif
    DegVec pivdeg; // pivot degree, first call
    DegVec pivdeg2; // pivot degree, second call
    DegVec rdeg(pmat.NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> trunc_pmat; // truncated pmat for first call
    Mat<zz_pX> appbas2; // basis for second call
    Mat<zz_pX> residual; // for the residual

    // first recursive call, with 'pmat' and 'shift'
#ifdef PMBASIS_PROFILE
    t1 = GetWallTime();
#endif
    trunc(trunc_pmat,pmat,order1);
    pivdeg = pmbasis(appbas,trunc_pmat,order1,shift);

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(first-call): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif
    // residual = (appbas * pmat * X^-order1) mod X^order2
    middle_product(residual, appbas, pmat, order1, order2-1);
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(middle-prod): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif

    // second recursive call, with 'residual' and 'rdeg'
    pivdeg2 = pmbasis(appbas2,residual,order2,rdeg);

    // final basis = appbas2 * appbas
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(second-call): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif
    multiply(appbas,appbas2,appbas);
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(basis-mul): " << (t2-t1) << "s" << std::endl;
#endif

    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;
}

/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis returning Popov                */
/*------------------------------------------------------------*/
DegVec popov_pmbasis(
                     Mat<zz_pX> &appbas,
                     const Mat<zz_pX> & pmat,
                     const long order,
                     const Shift & shift
                    )
{
    DegVec pivdeg = pmbasis(appbas,pmat,order,shift);
    Shift new_shift( pivdeg );
    std::transform(new_shift.begin(), new_shift.end(), new_shift.begin(), std::negate<long>());
    clear(appbas);
    pmbasis(appbas,pmat,order,new_shift);
    Mat<zz_p> lmat;
    leading_matrix(lmat, appbas, new_shift, true);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas);
    return pivdeg;
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PMBASIS-GENERIC                                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

DegVec popov_mbasis1_generic(
                             Mat<zz_pX> & appbas,
                             const Mat<zz_pX> & pmat,
                             const long order,
                             const Shift & shift
                             //Mat<zz_p> & kerbas,
                             //const Mat<zz_p> & pmat,
                             //const Shift & shift
                            )
{
    Mat<zz_pX> L;
    column_partial_linearization(L, pmat, order-1, 0);

    long m = L.NumRows();
    long n = L.NumCols();

    // copy the constant coeffs
    Mat<zz_p> L_cp;
    L_cp.SetDims(m,n);
    for (long r = 0; r < m; r++)
        for (long c = 0; c < n; c++)
            L_cp[r][c] = ConstTerm(L[r][c]);


    Mat<zz_p> inv_L;
    inv(inv_L, L_cp);

    Mat<zz_p> Ls;
    Ls.SetDims(L.NumRows(), L.NumCols());
    for (long c = 1; c < L.NumCols(); c++)
        for (long r = 0; r < L.NumRows(); r++)
            Ls[r][c] = L_cp[r][c-1];

    L_cp = Ls*inv_L;
    zz_pX x;
    SetCoeff(x,1,1);

    appbas.SetDims(L.NumRows(), L.NumCols());
    for (long r = 0; r < L.NumRows(); r++)
    {
        for (long c = 0; c < L.NumCols(); c++)
        {
            zz_pX tmp;
            SetCoeff(tmp, 0, -L_cp[r][c]);
            if (r == c) tmp += x;
            appbas[r][c] = tmp;
        }
    }
    DegVec dv (L.NumRows());
    for (long i = 0; i < L.NumRows(); i++)
        dv[i] = 1;
    return dv;
}


DegVec popov_mbasis1_generic2(
                             Mat<zz_pX> & appbas,
                             const Mat<zz_pX> & pmat,
                             const long order,
                             const Shift & shift
                             //Mat<zz_p> & kerbas,
                             //const Mat<zz_p> & pmat,
                             //const Shift & shift
                            )
{
    Mat<zz_pX> L;
    column_partial_linearization(L, pmat, order-1, 0);

    long m = L.NumRows();
    long n = L.NumCols();
            
    Mat<zz_p> LLs; // contains entries of L over Ls
    LLs.SetDims(2*L.NumRows(), L.NumCols());
    // copy L
    for (long r = 0; r < L.NumRows(); r++)
        for (long c = 0; c < L.NumCols(); c++)
            LLs[r][c] = ConstTerm(L[r][c]);
    // copy Ls
    for (long c = 1; c < L.NumCols(); c++)
        for (long r = 0; r < L.NumRows(); r++)
            LLs[r+L.NumRows()][c] = LLs[r][c-1];

    Mat<zz_p> kern;
    kernel(kern,LLs);
    
    zz_pX x;
    SetCoeff(x,1,1);
            
    appbas.SetDims(L.NumRows(), L.NumCols());
    for (long r = 0; r < L.NumRows(); r++)
    {
        for (long c = 0; c < L.NumCols(); c++)
        {
            zz_pX tmp;
            SetCoeff(tmp, 0, kern[r][c]);
            if (r == c) tmp += x;
            appbas[r][c] = tmp;
        }
    }
    
    DegVec dv (L.NumRows());
    for (long i = 0; i < L.NumRows(); i++)
        dv[i] = 1;
    return dv;
}


DegVec mbasis_generic(
                      Mat<zz_pX> & appbas,
                      const Mat<zz_pX> & pmat,
                      const long order,
                      const Shift & shift
                     )
{
    // TODO
    DegVec dv (pmat.NumRows());
    return dv;
}



/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis                                */
/*------------------------------------------------------------*/
DegVec pmbasis_generic(
                       Mat<zz_pX> & appbas,
                       const Mat<zz_pX> & pmat,
                       const long order,
                       const Shift & shift
                      )
{
#ifdef PMBASIS_PROFILE
    std::cout << "Order: " << order << std::endl;
    double t1,t2;
#endif
    // TODO thresholds to be determined:
    //  --> from mbasis (only linalg) to pmbasis with low-degree polmatmul (Karatsuba...)
    //  --> from this pmbasis to pmbasis with eval-based polmatmul (FFT, geometric..)
#ifdef PMBASIS_PROFILE
    if (order <= 32)
    {
        t1 = GetWallTime();
        DegVec rdeg = mbasis_vector(appbas,pmat,order,shift);
        t2 = GetWallTime();
        std::cout << "\tTime(base-case): " << (t2-t1) << "s" << std::endl;
        return rdeg;
    }
#else
    if (order == pmat.NumRows()){
        //cout << "pmat: " << pmat << endl;
        //cout << "blah blah" << endl;
        //popov_mbasis1_generic2(appbas,pmat,order,shift);
        return popov_mbasis1_generic(appbas,pmat,order,shift);
        //cout << "appbas1: " << appbas << endl;
        //auto t = mbasis_vector(appbas,pmat,order,shift);
        //cout << "appbas2: " << appbas << endl;
        //return t;
    }
#endif

    DegVec pivdeg; // pivot degree, first call
    DegVec pivdeg2; // pivot degree, second call
    DegVec rdeg(pmat.NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> trunc_pmat; // truncated pmat for first call
    Mat<zz_pX> appbas2; // basis for second call
    Mat<zz_pX> residual; // for the residual

    // first recursive call, with 'pmat' and 'shift'
#ifdef PMBASIS_PROFILE
    t1 = GetWallTime();
#endif
    trunc(trunc_pmat,pmat,order1);
    pivdeg = pmbasis_generic(appbas,trunc_pmat,order1,shift);

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(first-call): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif
    // residual = (appbas * pmat * X^-order1) mod X^order2
    long deg_sp = (pmat.NumCols() * order)/ (2*pmat.NumRows());
    right_parlin_middle_product(residual, appbas, pmat, order-1, deg_sp, order1, order2);
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(middle-prod): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif

    // second recursive call, with 'residual' and 'rdeg'
    pivdeg2 = pmbasis_generic(appbas2,residual,order2,rdeg);

    // final basis = appbas2 * appbas
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(second-call): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif
    multiply(appbas,appbas2,appbas);
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(basis-mul): " << (t2-t1) << "s" << std::endl;
#endif

    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;
}

DegVec pmbasis_generic2(
                       Mat<zz_pX> & appbas,
                       const Mat<zz_pX> & pmat,
                       const long order,
                       const Shift & shift
                      )
{
#ifdef PMBASIS_PROFILE
    std::cout << "Order: " << order << std::endl;
    double t1,t2;
#endif
    // TODO thresholds to be determined:
    //  --> from mbasis (only linalg) to pmbasis with low-degree polmatmul (Karatsuba...)
    //  --> from this pmbasis to pmbasis with eval-based polmatmul (FFT, geometric..)
#ifdef PMBASIS_PROFILE
    if (order <= 32)
    {
        t1 = GetWallTime();
        DegVec rdeg = mbasis_vector(appbas,pmat,order,shift);
        t2 = GetWallTime();
        std::cout << "\tTime(base-case): " << (t2-t1) << "s" << std::endl;
        return rdeg;
    }
#else
    if (order == pmat.NumRows()){
        //cout << "pmat: " << pmat << endl;
        //cout << "blah blah" << endl;
        return popov_mbasis1_generic2(appbas,pmat,order,shift);
        //return popov_mbasis1_generic(appbas,pmat,order,shift);
        //cout << "appbas1: " << appbas << endl;
        //auto t = mbasis_vector(appbas,pmat,order,shift);
        //cout << "appbas2: " << appbas << endl;
        //return t;
    }
#endif

    DegVec pivdeg; // pivot degree, first call
    DegVec pivdeg2; // pivot degree, second call
    DegVec rdeg(pmat.NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> trunc_pmat; // truncated pmat for first call
    Mat<zz_pX> appbas2; // basis for second call
    Mat<zz_pX> residual; // for the residual

    // first recursive call, with 'pmat' and 'shift'
#ifdef PMBASIS_PROFILE
    t1 = GetWallTime();
#endif
    trunc(trunc_pmat,pmat,order1);
    pivdeg = pmbasis_generic(appbas,trunc_pmat,order1,shift);

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(first-call): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif
    // residual = (appbas * pmat * X^-order1) mod X^order2
    long deg_sp = (pmat.NumCols() * order)/ (2*pmat.NumRows());
    right_parlin_middle_product(residual, appbas, pmat, order-1, deg_sp, order1, order2);
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(middle-prod): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif

    // second recursive call, with 'residual' and 'rdeg'
    pivdeg2 = pmbasis_generic(appbas2,residual,order2,rdeg);

    // final basis = appbas2 * appbas
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(second-call): " << (t2-t1) << "s" << std::endl;
    t1 = GetWallTime();
#endif
    multiply(appbas,appbas2,appbas);
#ifdef PMBASIS_PROFILE
    t2 = GetWallTime();
    std::cout << "\tTime(basis-mul): " << (t2-t1) << "s" << std::endl;
#endif

    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
