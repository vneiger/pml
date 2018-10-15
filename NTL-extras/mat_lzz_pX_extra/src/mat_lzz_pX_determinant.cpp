#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"

bool verify_determinant(const zz_pX & det, const Mat<zz_pX> & pmat, bool up_to_constant, bool randomized)
{
    if (not randomized)
        throw std::logic_error("==verify_determinant== *Deterministic* polynomial matrix determinant verification not implemented yet");


    // take a random field element and evaluate
    zz_p pt = random_zz_p();
    zz_p det_pt = determinant(eval(pmat, pt));
    // evaluate the provided determinant
    zz_p ev_det = eval(det,pt);

    // if allowing mult by a constant factor, do a first attempt
    // scale determinant to ensure its constant coefficient is correct
    if (not up_to_constant)
    {
        return (det_pt == ev_det);
    } 
    else
    {
        // take another random field element and evaluate
        pt = random_zz_p();
        zz_p det_pt2 = determinant(eval(pmat, pt));
        // evaluate the provided determinant
        zz_p ev_det2 = eval(det,pt);
        return (det_pt * ev_det2 == det_pt2 * ev_det);
    }
}

void determinant_via_diagonal_of_hermite(zz_pX & det, const Mat<zz_pX> & pmat)
{
    Vec<zz_pX> diag;
    diagonal_of_hermite(diag, pmat);
    set(det); // det = 1
    for (long i = 0; i < diag.length(); ++i)
        det *= diag[i];
}

bool determinant_generic_knowing_degree(zz_pX & det, const Mat<zz_pX> & pmat, long degree)
{
    long dim = pmat.NumRows();
    if (dim==1)
    {
        det = pmat[0][0];
        return (degree==deg(det));
    }
    else
    {
        long cdim1 = (dim>>1);  // cdim1 ~ dim/2
        long cdim2 = dim-cdim1;  // cdim2 ~ dim/2, cdim1+cdim2 = dim
        // TODO write and use "generic" kernel
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
        // TODO is computing the degree at each recursion level necessary?
        // (goes over the whole matrix each time... info could be transmitted
        // through the successive calls)
        Mat<zz_pX> appbas;
        // degree of kernel basis will be (generically)  D = dim * deg(pmat_l) / (dim - cdim1)
        // --> compute approximants at order deg(pmat_l) + D + 1
        // (cf for example Neiger-Rosenkilde-Solomatov, Lemma 4.3) <-- find better reference
        long deg_pmat_l = deg(pmat_l);
        long deg_ker = ceil( cdim1 * deg(pmat_l) / (double)(dim-cdim1) );
        long order = deg_pmat_l + deg_ker + 1;

        //std::cout << "computing approx basis of matrix with degrees:" << std::endl;
        //std::cout << degree_matrix(pmat_l) << std::endl;
        //std::cout << "and at order " << order << std::endl;
        Shift shift(dim,0);
        pmbasis(appbas, pmat_l, order, shift);

        // FIXME debug
        //std::cout << "appbas" << std::endl << degree_matrix(appbas) << std::endl;
        Mat<zz_pX> prod;
        multiply(prod,appbas,pmat);
        //std::cout << "prod" << std::endl << degree_matrix(prod) << std::endl;

        // minimal left kernel basis of pmat_r : last rows of app
        Mat<zz_pX> kerbas;
        kerbas.SetDims(cdim2,dim);
        for (long i = 0; i < cdim2; ++i)
            for (long j = 0; j < dim; ++j)
                kerbas[i][j] = appbas[i+cdim1][j];

        // then compute the product
        Mat<zz_pX> pmatt;
        multiply(pmatt, kerbas, pmat_r);

        //std::cout << "pmatt" << std::endl << degree_matrix(pmatt) << std::endl;

        return determinant_generic_knowing_degree(det,pmatt,degree);
    }

//	// recursion
//	size_t m1 = m>>1; // m1 ~ m/2
//	size_t m2 = m - m1; // m2 ~ m/2, m1 + m2 = m
//	vector<int> shift( m, 0 );
//	const size_t order = 2*sz-1;
//
//	// left columns of the matrix
//	PolMatDom::PMatrix series( this->field(), m, m1, order );
//	for ( size_t k=0; k<sz; ++k )
//	for ( size_t i=0; i<m; ++i )
//	for ( size_t j=0; j<m1; ++j )
//		series.ref(i,j,k) = pmat.get(i,j,k);
//
//	PolMatDom::PMatrix appbas( this->field(), m, m, 0 );
//	vector<int> rdeg = this->pmbasis( appbas, series, order, shift );
//
//	vector<int>::iterator it = rdeg.begin();
//	advance(it,m1);
//	size_t sz2 = *max_element(it,rdeg.end())+1;
//
//	// product
//	PolMatDom::PMatrix newpmat( this->field(), m2, m2, sz2+sz-1 );
//	// right half columns of pmat_sub
//	PolMatDom::PMatrix pmat_sub( this->field(), m, m2, sz );
//	for ( size_t k=0; k<sz; ++k )
//	for ( size_t i=0; i<m; ++i )
//	for ( size_t j=0; j<m2; ++j )
//		pmat_sub.ref(i,j,k) = pmat.get(i,j+m1,k);
//	// bottom half columns of appbas
//	PolMatDom::PMatrix appbas_sub( this->field(), m2, m, sz2 );
//	for ( size_t k=0; k<sz2; ++k )
//	for ( size_t i=0; i<m2; ++i )
//	for ( size_t j=0; j<m; ++j )
//		appbas_sub.ref(i,j,k) = appbas.get(i+m1,j,k);
//	// compute product
//
//	this->_PMMD.mul( newpmat, appbas_sub, pmat_sub );
//
//	return this->last_diagonal_entry( newpmat );
//
    
}

// TODO first version; improve
// TODO use a more general linsolve instead of specific algo
void determinant_via_linsolve(zz_pX & det, const Mat<zz_pX> & pmat)
{
    const long dim = pmat.NumRows();
    const long degree = dim*deg(pmat); // expected degree of determinant
    Mat<zz_pX> rand_vec;
    random(rand_vec, dim, 1, deg(pmat)+1);

    Mat<zz_pX> sol;
    solve_series_high_precision(sol, pmat, rand_vec, 2*degree+1);

    zz_pX elt;
    reverse(elt, sol[0][0], 2*degree);
    MinPolySeq(det, elt.rep, degree);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
