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
// for the general case


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

/*------------------------------------------------------------*/
/* ALGORITHMS FOR UNIFORM ORDER                               */
/*------------------------------------------------------------*/
// defined for arbitrary shifts, but
// works best for shifts close to uniform

std::vector<long> popov_mbasis1(
		Mat<zz_p> & kerbas,
		const Mat<zz_p> & pmat,
		const std::vector<long> & shift
		)
{
	// compute permutation which realizes stable sort of the shift
	// (i.e. sorts (shift[0],0)....(shift[len],len) lexicographically increasingly)
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

	// find the permuted kernel basis in row echelon form
	Mat<zz_p> p_kerbas;
	kernel(p_kerbas,mat);

	// compute the (permuted) pivot indices
	// FIXME unfortunately NTL doesn't return the pivot indices in Gaussian elimination
	// elimination... investigate if retrieving them as below actually takes time
	std::vector<long> p_pivind(p_kerbas.NumRows()); // pivot indices in permuted kernel basis
	p_pivind.back() = p_kerbas.NumCols()-1;
	for (long i = p_kerbas.NumRows()-1; i>=0; --i)
	{
		while (p_pivind[i]>=0 && p_kerbas[i][p_pivind[i]]==0)
			--p_pivind[i];
		// pivot of row i-1 is always less than pivot of row i
		if (i>0)
			p_pivind[i-1] = p_pivind[i]-1;
	}
	// note: "generically", the while loop will not perform any iteration

	// permute everything back to original order:
	// prepare kernel permutation by permuting kernel pivot indices
	// pivot degrees corresponding to kernel pivot indices are 0, others are 1
	std::vector<long> pivdeg(pmat.NumRows(),1);
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

	return pivdeg;
}

std::vector<long> mbasis(
		Mat<zz_pX> & appbas,
		const Mat<zz_pX> & pmat,
		const long order,
		const std::vector<long> & shift
		)
{
	// initially, appbas is the identity matrix
	appbas.SetDims(pmat.NumRows(),pmat.NumRows());
	for (long i = 0; i < appbas.NumRows(); ++i)
		SetCoeff(appbas[i][i],0);

	// holds the current shifted row degree of appbas
	// initially, this is exactly shift
	std::vector<long> rdeg( shift );

	// temporary buffer for local pivdegs at each mbasis1 call
	std::vector<long> pivdeg(pmat.NumRows());

	// matrix to store the kernels in mbasis1 calls
	Mat<zz_p> kerbas;
	// matrix to store residuals, initially constant coeff of pmat
	Mat<zz_p> residual( coeff(pmat,0) );

	for (long ord = 0; ord < order; ++ord)
	{
		// call MBasis1 to retrieve kernel and pivdeg
		pivdeg = popov_mbasis1(kerbas,residual,rdeg);

		// update approximant basis
		for (long i = 0; i < appbas.NumRows(); ++i) {
			if (pivdeg[i]==1) {
				LeftShiftRow(appbas,appbas,i,1);
			} else {
				
			}
		}

		// find new residual
		
		// TODO: other option with continuous full update of pmat

		// new shifted row degree = old rdeg + pivdeg  (entrywise)
		std::transform(rdeg.begin(), rdeg.end(), pivdeg.begin(), rdeg.begin(), std::plus<long>());
	}

	std::transform(rdeg.begin(),rdeg.end(),shift.begin(),rdeg.begin(),std::minus<long>());
	return rdeg;
}
