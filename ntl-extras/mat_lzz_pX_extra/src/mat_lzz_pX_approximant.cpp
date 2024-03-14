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
                       VecLong & rdeg,
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
        rdeg = shift;
        if (form<PolMatForm::POPOV)
        {
            appbas_iterative(appbas,pmat,order,rdeg);
            return;
        }
        if (form==PolMatForm::POPOV)
        {
            popov_appbas_iterative(appbas,pmat,order,rdeg);
            return;
        }
        throw std::logic_error("==approximant_basis== PolMatForm > POPOV: Not implemented yet");
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
        return false;

    // test whether appbas is shift-reduced with form at least 'form'
    if (not is_row_polmatform(form,appbas,shift))
        return false;

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
                return false;
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
        return false;

    // generation test: verify that [ cmat  P(0) ] has full rank (see Giorgi-Neiger ISSAC 2018)
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < m; ++j)
            cmat[i][j+n] = coeff(appbas[i][j],0);
    long rank = gauss(cmat);
    if (rank != m)
        return false;

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
void appbas_iterative(
                      Mat<zz_pX> & appbas,
                      const Mat<zz_pX> & pmat,
                      const VecLong & order,
                      VecLong & shift,
                      bool order_wise
                     )
{
    const long rdim = pmat.NumRows();
    const long cdim = pmat.NumCols();

    // initial approximant basis: identity of dimensions 'rdim x rdim'
    ident(appbas,rdim);

    // initial residual: the whole input matrix
    Mat<zz_pX> residual(pmat);

    // order that remains to be dealt with
    VecLong rem_order(order);

    // indices of columns/orders that remain to be dealt with
    VecLong rem_index(cdim);
    std::iota(rem_index.begin(), rem_index.end(), 0);

    // all along the algorithm, shift = shifted row degrees of approximant basis
    // (initially, input shift = shifted row degree of the identity matrix)

    while (not rem_order.empty())
    {
        /** Invariant:
         *  - appbas is a shift-ordered weak Popov approximant basis for
         *  (pmat,reached_order) where doneorder is the tuple such that
         *  -->reached_order[j] + rem_order[j] == order[j] for j appearing in rem_index
         *  -->reached_order[j] == order[j] for j not appearing in rem_index
         *  - shift == the "input shift"-row degree of appbas
         *  - residual == submatrix of columns (appbas * pmat)[:,j] for all j such that reached_order[j] < order[j]
         */

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
                if (piv<0 || shift[i] < shift[piv])
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
            ++shift[piv]; // shifted row degree of row piv increases
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
}

/*------------------------------------------------------------*/
/* specific to 2 x 1 input                                    */
/*------------------------------------------------------------*/
void appbas_iterative_2x1(
                           zz_pX & p00,
                           zz_pX & p01,
                           zz_pX & p10,
                           zz_pX & p11,
                           const zz_pX & f0,
                           const zz_pX & f1,
                           long order,
                           long & s0,
                           long & s1
                          )
{
    // initial approximant basis: identity of dimensions '2 x 2'
    p00 = 1;
    p01 = 0;
    p10 = 0;
    p11 = 1;

    // initial residual: the whole input f,g
    zz_pX r0 = f0;
    zz_pX r1 = f1;
    // will store the considered low-degree coefficient of residual,
    // and their quotient when applicable
    zz_p c0,c1,c;

    // all along the algorithm, shift (s0,s1) = shifted row degrees of approximant basis
    // (initially, input shift = shifted row degree of the identity matrix)
    for (long ord=0; ord<order; ++ord)
    {
        /** Invariant:
         *  - appbas is a shift-ordered weak Popov approximant basis for (pmat,ord)
         *    at the beginning of the loop, and (pmat,ord+1) at the end of the loop
         *  - (s0,s1) == the "input shift"-row degree of appbas
         *  - (r0,r1) == (appbas * [f,g]^t)
         */
        c0 = coeff(r0,ord); c1 = coeff(r1,ord);
        if (not IsZero(c0) && not IsZero(c1))
        {
            // most expected case: both coeffs are nonzero
            if (s0 <= s1)
            {
                // s0 <= s1 ==> use row 1 as pivot
                c = -c1/c0;
                p10 = p10 + c*p00;
                p11 = p11 + c*p01;
                p00 <<= 1; p01 <<= 1;
                r1 = r1 + c*r0;
                r0 <<= 1;
                s0 += 1;
            }
            else
            {
                // s0 > s1 ==> use row 2 as pivot
                c = -c0/c1;
                p00 = p00 + c*p10;
                p01 = p01 + c*p11;
                p10 <<= 1; p11 <<= 1;
                r0 = r0 + c*r1;
                r1 <<= 1;
                s1 += 1;
            }
        }
        else if (not IsZero(c1))
        {
            // multiply row 2 by x, do nothing on row 1
            p10 <<= 1; p11 <<= 1;
            r1 <<= 1;
            s1 += 1;
        }
        else if (not IsZero(c0))
        {
            // multiply row 1 by x, do nothing on row 2
            p00 <<= 1; p01 <<= 1;
            r0 <<= 1;
            s0 += 1;
        }
        // else, both coeffs are zero, do nothing
    }
}

/*------------------------------------------------------------*/
/* general case, output in Popov form                         */
/*------------------------------------------------------------*/
void popov_appbas_iterative(
                            Mat<zz_pX> & appbas,
                            const Mat<zz_pX> & pmat,
                            const VecLong & order,
                            VecLong & shift,
                            bool order_wise
                           )
{
    // compute first basis, copy the shift because we need it afterwards
    VecLong rdeg(shift);
    appbas_iterative(appbas,pmat,order,rdeg,order_wise);

    // shift for second call: negated pivot degree
    VecLong popov_shift(pmat.NumRows());
    std::transform(shift.begin(), shift.end(), rdeg.begin(),
                   popov_shift.begin(), std::minus<long>());

    // output shifted row degree
    shift=rdeg;

    // save `popov_shift` using `rdeg` as a buffer
    rdeg=popov_shift;

    // second call, basis shifted Popov up to constant transformation
    appbas_iterative(appbas,pmat,order,popov_shift,order_wise);

    // perform the constant transformation
    Mat<zz_p> lmat;
    row_leading_matrix(lmat, appbas, rdeg);
    inv(lmat, lmat);
    mul(appbas,lmat,appbas);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
