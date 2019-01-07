#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "mat_lzz_pX_approximant.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ALGORITHMS FOR APPROXIMANT BASES                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* GENERAL USER INTERFACE                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
// TODO check shift and order length, pmat dims (positive), ...
void approximant_basis(
                         Mat<zz_pX> & appbas,
                         const Mat<zz_pX> & pmat,
                         const VecLong & order,
                         const VecLong & shift,
                         const PolMatForm form,
                         const bool row_wise,
                         const bool generic
                        )
{
    std::cout << "==approximant_basis== WARNING: SLOW ITERATIVE ALGO" << std::endl;
    std::cout << "==approximant_basis== NOT READY FOR USE YET" << std::endl;
    if (row_wise && not generic)
    {
        VecLong pivdeg;
        if (form<PolMatForm::POPOV)
            pivdeg = appbas_iterative(appbas,pmat,order,shift);
        else if (form==PolMatForm::POPOV)
            pivdeg = popov_appbas_iterative(appbas,pmat,order,shift);
    }
    else
        throw std::logic_error("==approximant_basis== Not implemented yet");
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* VERIFICATION                                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// follows ideas from Algorithm 1 in Giorgi-Neiger, ISSAC 2018
bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const VecLong & order,
                          const VecLong & shift,
                          const PolMatForm & form,
                          const bool randomized
                         )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // TODO far from optimal except in the balanced case
    // (e.g. could be improved when deg(appbas)<<deg(pmat) like in Hermite-Pade,
    // or when appbas has strange column degrees or strange row degrees)

    if (not randomized)
        throw std::logic_error("==is_approximant_basis== Deterministic approximant basis verification not implemented.");

    // test whether appbas has the right dimensions
    if (appbas.NumRows() != appbas.NumCols() || appbas.NumCols() != m)
    {
        std::cout << "~~is_approx~~ wrong dimensions" << std::endl;
        return false;
    }

    // test whether appbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,appbas,shift))
    {
        std::cout << "~~is_approx~~ not in the required shifted form" << std::endl;
        std::cout << degree_matrix(appbas) << std::endl;
        return false;
    }

    // test whether appbas consists of approximants (if row-wise: appbas * pmat = 0 mod X^order)
    // and retrieve the constant coefficient "cmat" of degree order (if row-wise: cmat = (appbas * pmat * X^{-order})  mod X)
    // (we reserve some additional space in cmat because later it will store the constant coefficient of appbas)
    Mat<zz_pX> residual;
    multiply(residual,appbas,pmat);
    // TODO this multiplication could be:
    //   - truncated mod X^{order+1}
    //   - improved by taking degree profile into account

    Mat<zz_p> cmat;
    cmat.SetDims(m,m+n);

    for (long i = 0; i < m; ++i)
    {
        for (long j = 0; j < n; ++j)
        {
            long ord = order[j];
            GetCoeff(cmat[i][j],residual[i][j],ord);
            trunc(residual[i][j],residual[i][j],ord);
            if (not IsZero(residual[i][j]))
            {
                std::cout << "~~is_approx~~ not approx (" << i << "," << j << "): " << std::endl;
                return false;
            }
        }
    }

    // test whether det(appbas) is a monomial (power of x):
    // det(appbas) = det(appbas(1)) * X^(deg(det(appbas)))
    // the degree of det(appbas) is known since appbas is reduced:
    VecLong degs;
    row_degree(degs, appbas, shift);
    long sum_degs = std::accumulate(degs.begin(), degs.end(), (long)0);
    long sum_shift = std::accumulate(shift.begin(), shift.end(), (long)0);
    long degdet = sum_degs - sum_shift;
    // for a random point pt, verify that det(appbas(pt)) = det(appbas(1)) * pt^degdet
    zz_p pt = random_zz_p();
    zz_p det_pt = determinant(eval(appbas, pt));
    zz_p det_one = determinant(eval(appbas, zz_p(1)));
    if (det_pt != det_one * power(pt,degdet))
    {
        std::cout << "~~is_approx~~ wrong det" << std::endl;
        return false;
    }

    // generation test: verify that [ cmat  P(0) ] has full rank (see Giorgi-Neiger ISSAC 2018)
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < m; ++j)
            cmat[i][j+n] = coeff(appbas[i][j],0);
    long rank = gauss(cmat);
    if (rank != m)
    {
        std::cout << "~~is_approx~~ cmat not full rank" << std::endl;
        return false;
    }

    return true;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ITERATIVE ALGORITHMS                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* for the general case                                       */
/*------------------------------------------------------------*/
// TODO not optimized at all; is it ever faster than mbasis and others?
// (maybe with highly non-uniform orders/shifts?)
VecLong appbas_iterative(
                        Mat<zz_pX> & appbas,
                        const Mat<zz_pX> & pmat,
                        const VecLong & order,
                        const VecLong & shift,
                        bool order_wise
                       )
{
    long rdim = pmat.NumRows();
    long cdim = pmat.NumCols();

    // initial approximant basis: identity of dimensions 'rdim x rdim'
    ident(appbas,rdim);

    // initial residual: the whole input matrix
    Mat<zz_pX> residual(pmat);

    // order that remains to be dealt with
    VecLong rem_order(order);

    // indices of columns/orders that remain to be dealt with
    VecLong rem_index(cdim);
    std::iota(rem_index.begin(), rem_index.end(), 0);

    // shifted row degrees of approximant basis
    // (initially, of the identity matrix, i.e. rdeg == shift)
    VecLong rdeg(shift);

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
            j = std::distance(rem_order.begin(), std::max_element(rem_order.begin(), rem_order.end()));

        long deg = order[rem_index[j]] - rem_order[j];

        // record the coefficients of degree deg of the column j of residual
        // also keep track of which of these are nonzero,
        // and among the nonzero ones, which is the first with smallest shift
        Vec<zz_p> const_residual;
        const_residual.SetLength(rdim);
        VecLong indices_nonzero;
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
            LeftShiftRow(appbas, appbas, piv, 1); // row piv multiplied by X
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
                Mat<zz_pX> buffer(residual);
                residual.SetDims(rdim,rem_order.size()); // storage is freed+reallocated
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
VecLong popov_appbas_iterative(
                              Mat<zz_pX> & appbas,
                              const Mat<zz_pX> & pmat,
                              const VecLong & order,
                              const VecLong & shift,
                              bool order_wise
                             )
{
    // FIXME: first call can be very slow if strange degrees (can it? find example)
    // --> rather implement BecLab00's "continuous" normalization?
    VecLong pivdeg = appbas_iterative(appbas,pmat,order,shift,order_wise);
    VecLong new_shift(pivdeg);
    std::transform(new_shift.begin(), new_shift.end(), new_shift.begin(), std::negate<long>());
    appbas.kill();
    appbas_iterative(appbas,pmat,order,new_shift,order_wise);
    Mat<zz_p> lmat;
    row_leading_matrix(lmat, appbas, new_shift);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas); // FIXME special mult (expand columns of appbas)
    return pivdeg;
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
