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


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MINIMAL APPROXIMANT BASES ALGORITHMS                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* VERIFICATION                                               */
/*------------------------------------------------------------*/

bool is_approximant_basis(
		const Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const std::vector<long> & order,
		const std::vector<long> & shift,
		const PolMatForm & form,
		const bool row_wise,
		const bool randomized
		)
{
	if (randomized)
		throw std::logic_error("==is_approximant_basis== Fast randomized approximant basis not implemented yet");

	std::cout << "==is_approximant_basis== WARNING: not fully implemented: not checking generation" << std::endl;

	// test that appbas is shift-reduced with form at least 'form'
	if (not is_polmatform(appbas,form,shift,row_wise))
		return false;

	// test that the matrix consists of approximants
	Mat<zz_pX> residual;
	if (row_wise)
		multiply_naive(residual,appbas,pmat); // TODO this mul should be (truncated and) non-naive
	else
		multiply_naive(residual,pmat,appbas); // TODO this mul should be (truncated and) non-naive
	for (long i = 0; i < residual.NumRows(); ++i)
	{
		for (long j = 0; j < residual.NumCols(); ++j)
		{
			long ord = row_wise ? order[j] : order[i];
			trunc(residual[i][j],residual[i][j],ord);
			if (residual[i][j] != 0)
				return false;
		}
	}

	return true;
}

bool is_approximant_basis(
		const Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift,
		const PolMatForm & form,
		const bool row_wise,
		const bool randomized
		)
{
	std::vector<long> orders(pmat.NumRows(),order);
	return is_approximant_basis(appbas,pmat,orders,shift,form,row_wise,randomized);
}

/*------------------------------------------------------------*/
/* ITERATIVE ALGORITHMS                                       */
/*------------------------------------------------------------*/


std::vector<long> appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const std::vector<long> order,
		const std::vector<long> & shift,
		bool order_wise
		)
{
	/** Three possibilities (among others) for next coefficient to deal with:
	 *   - process 'pmat' order-wise (choose column with largest order)
	 *   - process 'pmat' column-wise (choose leftmost column not yet completed)
	 **/

	long rdim = pmat.NumRows();
	long cdim = pmat.NumCols();

	// initial approximant basis: identity of dimensions 'rdim x rdim'
	appbas.SetDims(rdim,rdim);
	for (long i = 0; i < rdim; ++i)
		appbas[i][i] = 1;

	// initial residual: the whole input matrix
	Mat<zz_pX> residual( pmat );

	// order that remains to be dealt with
	std::vector<long> rem_order( order );

	// indices of columns/orders that remain to be dealt with
	std::vector<long> rem_index( cdim );
	std::iota(rem_index.begin(), rem_index.end(), 0);

	// shifted row degrees of approximant basis
	// (initially, of the identity matrix, i.e. rdeg == shift)
	std::vector<long> rdeg( shift );

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
			for (long k=0; k<rdim; ++k)
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

std::vector<long> popov_appbas_iterative(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const std::vector<long> order,
		const std::vector<long> & shift,
		bool order_wise
		)
{
	// TODO: first call can be very slow if strange degrees --> rather implement
	// BecLab00's "continuous" normalization
	std::vector<long> pivdeg = appbas_iterative(appbas,pmat,order,shift,order_wise);
	std::transform(pivdeg.begin(), pivdeg.end(), pivdeg.begin(), std::negate<long>());
	pivdeg = appbas_iterative(appbas,pmat,order,pivdeg,order_wise);
	// TODO multiply appbas by inverse of leading matrix
	return pivdeg;
}
