#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota
#include <NTL/BasicThreadPool.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* interpolant basis verification                             */
/* TODO incomplete                                            */
/*------------------------------------------------------------*/
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & evals,
                          const Vec<zz_p> & pts,
                          const Shift & shift,
                          const PolMatForm & form,
                          const bool row_wise,
                          const bool randomized
                         )
{
    if (randomized)
        throw std::logic_error("==is_interpolant_basis== Fast randomized interpolant basis verification not implemented yet");

    std::cout << "==is_interpolant_basis== WARNING: not fully implemented: not checking generation" << std::endl;

    // test that appbas is shift-reduced with form at least 'form'
    if (not is_polmatform(intbas,form,shift,row_wise))
        return false;

    // test that the matrix consists of interpolants
    Vec<Mat<zz_p>> intbas_vec;
    conv(intbas_vec,intbas);
    long deg = intbas_vec.length()-1;
    for (long pt=0; pt<pts.length(); ++pt)
    {
        // evaluate with Horner
        Mat<zz_p> ev = intbas_vec[deg];
        for (long d = deg-1; d >= 0; --d)
        {
            ev = pts[pt] * ev;
            add(ev, ev, intbas_vec[d]);
        }
        mul(ev, ev, evals[pt]);
        if (not IsZero(ev))
            return false;
    }
    return true;
}



//DegVec intbas_iterative(
//		Mat<zz_pX> & intbas,
//		const Mat<zz_pX> & pmat,
//		const Points & pts,
//		const Shift & shift,
//		bool point_wise,
//		)
//{
//	/** Two possibilities (among others) for next coefficient to deal with:
//	 *   - process 'pmat' column-wise (choose leftmost column not yet completed)
//	 *   - process 'pmat' point-wise (choose points in order given in the list)
//	 **/
//
//	long rdim = pmat.NumRows();
//	long cdim = pmat.NumCols();
//
//	// initial approximant basis: identity of dimensions 'rdim x rdim'
//	ident(intbas,rdim);
//
//	// initial residual: the whole input matrix
//	Mat<zz_pX> residual( pmat );
//
//	// order that remains to be dealt with
//	Points rem_pts( pts );
//
//	// indices of columns that remain to be dealt with
//	std::vector<long> rem_index( cdim );
//	std::iota(rem_index.begin(), rem_index.end(), 0);
//
//	// shifted row degrees of approximant basis
//	// (initially, of the identity matrix, i.e. rdeg == shift)
//	DegVec rdeg(shift);
//
//	while (not rem_pts.empty())
//	{
//		/** Invariant:
//		 *  - intbas is a shift-ordered weak Popov approximant basis for
//		 *  (pmat,reached_order) where doneorder is the tuple such that
//		 *  -->reached_order[j] + rem_order[j] == order[j] for j appearing in rem_index
//		 *  -->reached_order[j] == order[j] for j not appearing in rem_index
//		 *  - rdeg == the shift-row degree of intbas
//		 *  - residual == submatrix of columns (intbas * pmat)[:,j] for all j such that reached_order[j] < order[j]
//		 **/
//
//		long j=0; // value if columnwise (order_wise==False)
//		if (point_wise)
//		{
//			// FIXME if seems to slow (e.g. compared to non-point-wise), could be
//			// valuable to simply initially permute the columns of pmat and order:
//			// then 'j' is obvious
//			j = // TODO
//		}
//
//		// TODO FIXME below here, not correct.. (copy paste from appbas)
//
//		long deg = order[rem_index[j]] - rem_order[j];
//
//		// record the coefficients of degree deg of the column j of residual
//		// also keep track of which of these are nonzero,
//		// and among the nonzero ones, which is the first with smallest shift
//		Vec<zz_p> const_residual;
//		const_residual.SetLength(rdim);
//		std::vector<long> indices_nonzero;
//		long piv = -1;
//		for (long i = 0; i < rdim; ++i)
//		{
//			const_residual[i] = coeff(residual[i][j],deg);
//			if (const_residual[i] != 0)
//			{
//				indices_nonzero.push_back(i);
//				if (piv<0 || rdeg[i] < rdeg[piv])
//					piv = i;
//			}
//		}
//
//		// if indices_nonzero is empty, const_residual is already zero, there is nothing to do
//		if (not indices_nonzero.empty())
//		{
//			// update all rows of intbas and residual in indices_nonzero except row piv
//			for (long row : indices_nonzero)
//			{
//				if (row!=piv)
//				{
//					zz_p c = - const_residual[row] / const_residual[piv];
//					for (long k=0; k<rdim; ++k)
//						intbas[row][k] += c * intbas[piv][k];
//					for (long k=0; k<residual.NumCols(); ++k)
//						residual[row][k] += c * residual[piv][k];
//				}
//			}
//
//			// update row piv
//			++rdeg[piv]; // shifted row degree of row piv increases
//			for (long k=0; k<rdim; ++k) // TODO use shiftRow
//				intbas[piv][k] <<= 1; // row piv multiplied by X
//			for (long k=0; k<residual.NumCols(); ++k)
//			{
//				residual[piv][k] <<= 1; // row piv multiplied by X
//				trunc(residual[piv][k],residual[piv][k],order[rem_index[k]]); // truncate
//			}
//		}
//
//		// now column j (or rather rem_index[j] in original 'pmat') is zero mod X^(deg+1)
//		if (rem_order[j] > 1) // one step of progress towards order[j]
//			--rem_order[j];
//		else // rem_order[j] == 1: work completed for column j
//		{
//			rem_order.erase(rem_order.begin() + j);
//			rem_index.erase(rem_index.begin() + j);
//			if (!rem_order.empty())
//			{
//				Mat<zz_pX> buffer( residual );
//				residual.kill();
//				residual.SetDims(rdim,rem_order.size());
//				for (long i=0; i<rdim; ++i)
//				{
//					for (long k=0; k<j; ++k)
//							residual[i][k] = buffer[i][k];
//					for (long k=j+1; k<buffer.NumCols(); ++k)
//						residual[i][k-1] = buffer[i][k];
//				}
//			}
//		}
//	}
//
//	// make rdeg contain the pivot degree rather than shifted row degree, i.e.:
//	// rdeg = rdeg-shift   , done entry-wise
//	std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::minus<long>());
//	return rdeg;
//}
//
//DegVec popov_intbas_iterative(
//		Mat<zz_pX> & intbas,
//		const Mat<zz_pX> & pmat,
//		const Points & pts,
//		const Shift & shift,
//		bool point_wise
//		)
//{
//	// TODO: implement BecLab00's "continuous" normalization and compare timings
//	// (for various shifts and point configurations)
//	DegVec pivdeg = intbas_iterative(intbas,pmat,pts,shift,point_wise);
//	Shift new_shift( pivdeg );
//	std::transform(new_shift.begin(), new_shift.end(), new_shift.begin(), std::negate<long>());
//	// TODO write zero method for polmats
//	for (long i = 0; i < intbas.NumCols(); ++i)
//		for (long j = 0; j < intbas.NumRows(); ++j)
//			intbas[i][j] = 0;
//	intbas_iterative(intbas,pmat,order,new_shift,point_wise);
//	Mat<zz_p> lmat;
//	leading_matrix(lmat, intbas, new_shift, true);
//	inv(lmat, lmat);
//	mul(intbas,lmat,intbas);
//	return pivdeg;
//}

DegVec mbasis(
              Mat<zz_pX> & intbas,
              const Vec<Mat<zz_p>> & evals,
              const Vec<zz_p> & pts,
              const Shift & shift
             )
{
    long nrows = evals[0].NumRows();

    // interpolation basis represented as a vector (polynomial) of constant matrices
    Vec<Mat<zz_p>> coeffs_intbas;

    // initially, intbas is the identity matrix
    coeffs_intbas.SetLength(1);
    ident(coeffs_intbas[0], nrows);

    // holds the current shifted row degree of coeffs_intbas
    // initially, this is exactly shift
    DegVec rdeg(shift);

    // holds the current pivot degree of coeffs_intbas
    // initially tuple of zeroes
    // (note that at all times pivdeg+shift = rdeg entrywise)
    DegVec pivdeg(nrows);

    // will store the pivot degree at each call of mbasis1
    DegVec diff_pivdeg;

    // matrix to store the kernels in mbasis1 calls
    Mat<zz_p> kerbas;
    // matrix to store residuals, initially constant coeff of evals
    Mat<zz_p> residual( evals[0] );
    Mat<zz_p> intbas_eval;
    intbas_eval.SetDims(nrows,nrows);

    // declare matrices
    Mat<zz_p> res_coeff; // will store coefficient matrices used to compute the residual
    Mat<zz_p> kerapp; // will store constant-kernel * coeffs_intbas[d]

    for (long pt = 1; pt <= pts.length(); ++pt)
    {
        // call MBasis1 to retrieve kernel and pivdeg
        diff_pivdeg = popov_mbasis1(kerbas,residual,rdeg);

        if (kerbas.NumRows()==0)
        {
            // the basis is simply multiplied by X-pt
            // new degree is old degree + 1, that is, the length
            long d = coeffs_intbas.length();
            coeffs_intbas.SetLength(d+1);
            if (d>0)
                coeffs_intbas[d] = coeffs_intbas[d-1];
            --d;
            for (; d > 0; --d)
            {
                coeffs_intbas[d] *= -pts[pt-1];
                coeffs_intbas[d] += coeffs_intbas[d-1];
            }
            coeffs_intbas[0] *= -pts[pt-1];
            // update pivdeg accordingly
            std::for_each(pivdeg.begin(), pivdeg.end(), [](long& a) { ++a; });

            // TODO what about residual??
        }

        // kerbas.NumRows()==residual.NumRows() --> interpolant basis is already
        // correct for this point, just go to the next

        else if (kerbas.NumRows()<residual.NumRows())
        {
            // I/ Update degrees:
            // new shifted row degree = old rdeg + diff_pivdeg
            std::transform(rdeg.begin(), rdeg.end(), diff_pivdeg.begin(), rdeg.begin(), std::plus<long>());
            // new pivot degree = old pivot_degree + diff_pivdeg
            std::transform(pivdeg.begin(), pivdeg.end(), diff_pivdeg.begin(), pivdeg.begin(), std::plus<long>());
            // deduce degree of coeffs_intbas; note that it is a property of this algorithm
            // that deg(coeffs_intbas) = max(pivot degree) (i.e. max(degree of diagonal
            // entries); this does not hold in general for ordered weak Popov forms
            long deg_coeffs_intbas = *std::max_element(pivdeg.begin(), pivdeg.end());

            // II/ update approximant basis

            // submatrix of rows with diff_pivdeg==0 is replaced by kerbas*coeffs_intbas
            // while rows with diff_pivdeg=1 are multiplied by X-pt
            // --> the loop goes downwards, so that we can do both in the same iteration
            long row;
            long d = deg_coeffs_intbas;
            // separate treatment of highest degree terms, if degree has increased
            if (d>=coeffs_intbas.length()) // then, it is equal
            {
                coeffs_intbas.SetLength(d+1); // increases length by 1
                coeffs_intbas[d].SetDims(nrows,nrows);
                for (long i = 0; i < nrows; ++i)
                    if (diff_pivdeg[i]==1)
                        coeffs_intbas[d][i] = coeffs_intbas[d-1][i];
                --d;
            }

            // normal treatment of other terms
            for (; d >= 0; --d)
            {
                mul(kerapp,kerbas,coeffs_intbas[d]);
                row=0;
                for (long i = 0; i < nrows; ++i)
                {
                    if (diff_pivdeg[i]==0)
                    {
                        coeffs_intbas[d][i] = kerapp[row];
                        ++row;
                    }
                    else  // diff_pivdeg[i]==1 --> multiply by X-pt
                    {
                        coeffs_intbas[d][i] *= -pts[pt-1];
                        if (d>0)
                            coeffs_intbas[d][i] += coeffs_intbas[d-1][i];
                    }
                }
            }

            // III/ compute next residual, if needed
            // it is coeffs_intbas(pts[pt]) * evals[pt]
            // evaluate with Horner method
            if (pt<pts.length())
            {
                intbas_eval = coeffs_intbas[deg_coeffs_intbas];
                d = deg_coeffs_intbas-1;
                for (; d >= 0; --d)
                {
                    intbas_eval = pts[pt] * intbas_eval;
                    add(intbas_eval, intbas_eval, coeffs_intbas[d]);
                }
                mul(residual,intbas_eval,evals[pt]);
            }
        }
    }

    conv(intbas,coeffs_intbas);
    return pivdeg;
}

DegVec pmbasis_geometric(
               Mat<zz_pX> & intbas,
               const Mat<zz_pX> & pmat,
               const zz_p & r,
               const long order,
               const Shift & shift,
               Vec<Mat<zz_p>> &evals,
               Vec<zz_p> &pts
              )
{
    DegVec rdeg;
    rdeg.resize(pmat.NumRows());
    row_degree(rdeg, pmat);
    long max_rowdeg = rdeg[0];
    for (auto i : rdeg)
        if (i > max_rowdeg) max_rowdeg = i;
    max_rowdeg = max(order, max_rowdeg);
    zz_pX_Multipoint_Geometric eval(r,max_rowdeg+1);
    
    // set up pts
    zz_pX x;
    SetCoeff(x,1,1);
    Vec<zz_p> pts2;
    eval.evaluate(pts2, x); // just gets powers of r
    
    // set up evaluations of pmat
    evals.SetLength(order);
    pts.SetLength(order);
    for (long d = 0; d < order; d++)
    {
        pts[d] = pts2[d];
        evals[d].SetDims(pmat.NumRows(), pmat.NumCols());
    }
    for (long r = 0; r < pmat.NumRows(); r++)
    {
        for (long c = 0; c < pmat.NumCols(); c++)
        {
            Vec<zz_p> vals;
            eval.evaluate(vals, pmat[r][c]);
            for (long d = 0; d < order; d++)
                evals[d][r][c] = vals[d];
        }
    }
    
    return pmbasis_geometric(intbas, evals, pts, r, shift);
}

DegVec pmbasis_geometric(
                         Mat<zz_pX> & intbas,
                         const Vec<Mat<zz_p>> & evals,
                         const Vec<zz_p> & pts,
                         const zz_p & r,
                         const Shift & shift
                        )
{
    const long order = pts.length();
    zz_pContext context;
    
    if (order <= 32)
        return mbasis(intbas, evals, pts, shift);

    DegVec pivdeg; // pivot degree, first call
    DegVec pivdeg2; // pivot degree, second call
    DegVec rdeg(evals[0].NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> intbas2; // basis for second call

    // first recursive call
    Vec<zz_p>  pts1;
    pts1.SetLength(order1);
    for (long i = 0; i < order1; i++)
        pts1[i] = pts[i];
    pivdeg = pmbasis_geometric(intbas, evals, pts1, r, shift);

    long max_pivdeg = pivdeg[0];
    for (auto i : pivdeg)
        if (max_pivdeg < i) max_pivdeg = i;
        
    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

    // get the product of evaluations intbas(x_i) * pmat(x_i)
    // for the second half of the points
    Vec<zz_p> pts2;
    pts2.SetLength(order2);
    for (long i=0; i<order2; ++i)
        pts2[i] = pts[order1+i];
        
    // geometric progression of r, starting at pts2[0]
    // first ensure we have enough points for the degree
    max_pivdeg = max(order2, max_pivdeg); 
    zz_pX_Multipoint_Geometric ev(r,pts2[0], max_pivdeg+1);
    
    Vec<Mat<zz_p>> evals2;
    evals2.SetLength(order2);
    for (long i = 0; i < order2; i++)
        evals2[i].SetDims(intbas.NumRows(),intbas.NumCols());
    
    context.save();  

    // evaluate and store
NTL_EXEC_RANGE(intbas.NumRows(),first,last)
    
    context.restore();    
    for (long r = first; r < last; r++)
        for (long c = 0; c < intbas.NumCols(); c++)
        {
            Vec<zz_p> val;
            ev.evaluate(val, intbas[r][c]);
            
            for (long i = 0; i < order2; i++)
                evals2[i][r][c] = val[i];
            
        }
NTL_EXEC_RANGE_END

    context.save();  

    // multiply and store    
NTL_EXEC_RANGE(order2,first,last)
    context.restore();    
    
    for (long i = first; i < last; i++)
    {
        evals2[i] = evals2[i] * evals[order1+i];
    }
NTL_EXEC_RANGE_END
    
    // second recursive call
    pivdeg2 = pmbasis_geometric(intbas2, evals2, pts2, r, rdeg);
    
    multiply(intbas,intbas2,intbas);
    
    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;    
}

DegVec pmbasis(
              Mat<zz_pX> & intbas,
              const Vec<Mat<zz_p>> & evals,
              const Vec<zz_p> & pts,
              const Shift & shift
             )
{

    const long order = pts.length();
    zz_pContext context;
    
    if (order <= 32)
        return mbasis(intbas, evals, pts, shift);

    DegVec pivdeg; // pivot degree, first call
    DegVec pivdeg2; // pivot degree, second call
    DegVec rdeg(evals[0].NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> intbas2; // basis for second call

    // first recursive call
    Vec<zz_p>  pts1;
    pts1.SetLength(order1);
    for (long i = 0; i < order1; i++)
        pts1[i] = pts[i];
    pivdeg = pmbasis(intbas, evals, pts1, shift);

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

    // get the product of evaluations intbas(x_i) * pmat(x_i)
    // for the second half of the points
    Vec<zz_p> pts2;
    pts2.SetLength(order2);
    for (long i=0; i<order2; ++i)
        pts2[i] = pts[order1+i];
        
    zz_pX_Multipoint_General ev(pts2);
    
    Vec<Mat<zz_p>> evals2;
    evals2.SetLength(order2);
    for (long i = 0; i < order2; i++)
        evals2[i].SetDims(intbas.NumRows(),intbas.NumCols());
    
    context.save();  

    // evaluate and store
NTL_EXEC_RANGE(intbas.NumRows(),first,last)
    
    context.restore();    
    for (long r = first; r < last; r++)
        for (long c = 0; c < intbas.NumCols(); c++)
        {
            Vec<zz_p> val;
            ev.evaluate(val, intbas[r][c]);
            for (long i = 0; i < order2; i++)
                evals2[i][r][c] = val[i];
            
        }
NTL_EXEC_RANGE_END

    context.save();  

    // multiply and store    
NTL_EXEC_RANGE(order2,first,last)
    context.restore();    
    
    for (long i = first; i < last; i++)
    {
        evals2[i] = evals2[i] * evals[order1+i];
    }
NTL_EXEC_RANGE_END
    
    // second recursive call
    pivdeg2 = pmbasis(intbas2, evals2, pts2, rdeg);
    
    multiply(intbas,intbas2,intbas);
    
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
