#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

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
//	intbas.SetDims(rdim,rdim);
//	for (long i = 0; i < rdim; ++i)
//		SetCoeff(intbas[i][i],0);
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
//	DegVec rdeg( shift );
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

// code inspired from mbasis_vector (input+output are Vec<Mat<zz_p>>)
DegVec mbasis(
		Vec<Mat<zz_p>> & intbas,
		const Vec<Mat<zz_p>> & evals,
		const Vec<zz_p> & pts,
		const Shift & shift
		)
{
	long nrows = evals[0].NumRows();
	// TODO FIXME expected degree in non-shifted case --> assuming generic and shift==0...
	long expected_degree = ceil(evals[0].NumCols()*pts.length()/nrows);
	intbas.SetLength(expected_degree+1);
	for (long d = 0; d < expected_degree+1; ++d)
		intbas[d].SetDims(nrows,nrows);
	// initially, intbas is the identity matrix
	long current_degree = 0;
	for (long i = 0; i < nrows; ++i)
		intbas[0].put(i,i,zz_p(1));

	// holds the current shifted row degree of intbas
	// initially, this is exactly shift
	DegVec rdeg( shift );

	// holds the current pivot degree of intbas
	// initially tuple of zeroes
	// (note that at all times pivdeg+shift = rdeg entrywise)
	DegVec pivdeg(nrows);

	// will store the pivot degree at each call of mbasis1
	DegVec diff_pivdeg;

	// matrix to store the kernels in mbasis1 calls
	Mat<zz_p> kerbas;
	// matrix to store residuals, initially constant coeff of evals
	Mat<zz_p> residual( evals[0] );

	// declare matrices
	Mat<zz_p> res_coeff; // will store coefficient matrices used to compute the residual
	Mat<zz_p> kerapp; // will store constant-kernel * intbas[d]

	for (long pt = 1; pt <= pts.length(); ++pt)
	{
		// call MBasis1 to retrieve kernel and pivdeg
		diff_pivdeg = popov_mbasis1(kerbas,residual,rdeg);

		if (kerbas.NumRows()==0)
		{
			// the basis is simply multiplied by X-pt
			++current_degree;
			intbas[current_degree] = intbas[current_degree-1];
			for (long d = current_degree-1; d > 0; --d)
			{
				intbas[d] *= -pts[pt-1];
				intbas[d] += intbas[d-1];
			}
			intbas[0] *= pts[pt-1];
			// update pivdeg accordingly
			std::for_each(pivdeg.begin(), pivdeg.end(), [](long& a) { ++a; });
		}

		// kerbas.NumRows()==residual.NumRows() --> interpolant basis is already
		// correct for this point, just go to the next

		if (kerbas.NumRows()<residual.NumRows())
		{
			// I/ Update degrees:
			// new shifted row degree = old rdeg + diff_pivdeg
			std::transform(rdeg.begin(), rdeg.end(), diff_pivdeg.begin(), rdeg.begin(), std::plus<long>());
			// new pivot degree = old pivot_degree + diff_pivdeg
			std::transform(pivdeg.begin(), pivdeg.end(), diff_pivdeg.begin(), pivdeg.begin(), std::plus<long>());
			// deduce degree of intbas; note that it is a property of this algorithm
			// that deg(intbas) = max(pivot degree) (i.e. max(degree of diagonal
			// entries); this does not hold in general for ordered weak Popov forms
			long deg_intbas = *std::max_element(pivdeg.begin(), pivdeg.end());

			if (deg_intbas > expected_degree) {
				std::cout << "~~mbasis-variant not implemented yet, currently assuming generic evals and uniform shift" << std::endl;
				throw;
			}

			// II/ update approximant basis

			// submatrix of rows with diff_pivdeg==0 is replaced by kerbas*intbas
			// while rows with diff_pivdeg=1 are multiplied by X-pt
			// --> the loop goes downwards, so that we can do both in the same iteration
			long row;
			// separate treatment of highest degree terms, if degree has increased
			if (deg_intbas>current_degree)
				for (long i = 0; i < nrows; ++i)
					if (diff_pivdeg[i]==1)
						intbas[deg_intbas][i] = intbas[deg_intbas-1][i];
			// normal treatment of other terms
			for (long d = current_degree; d >= 0; --d)
			{
				mul(kerapp,kerbas,intbas[d]);
				row=0;
				for (long i = 0; i < nrows; ++i)
				{
					if (diff_pivdeg[i]==0)
					{
						intbas[d][i] = kerapp[row];
						++row;
					}
					else  // diff_pivdeg[i]==1 --> multiply by X-pt
					{
						intbas[d][i] *= -pts[pt-1];
						if (d>0)
							intbas[d][i] = intbas[d-1][i];
					}
				}
			}

			// III/ compute next residual, if needed
			// it is intbas(pts[pt]) * evals[pt]
			// evaluate with Horner method
			if (pt<pts.length())
			{
				residual = intbas[deg_intbas];
				for (long d = deg_intbas-1; d >= 0; --d) // note that deg_intbas <= pt holds
				{
					residual = pts[pt] * residual;
					add(residual, residual, intbas[d]);
				}
				mul(residual,residual,evals[pt]);
			}
			current_degree = deg_intbas;
		}
	}

	return pivdeg;
}


