#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "mat_lzz_pX_extra.h"

NTL_CLIENT

template<typename T>
void printVec(std::vector<T> const &input)
{
	std::cout << '[';
	for (auto elt : input) {
		std::cout << elt << ", ";
	}
	std::cout << ']' << std::endl;
}

/******************************************************/
/* Tests the degree functions for polynomial matrices */
/******************************************************/

int main(){
	zz_p::init(13);
	Mat<zz_pX> pmat;
	pmat.SetDims(2,3);
	pmat[0][0] = random_zz_pX(3);
	pmat[0][1] = random_zz_pX(2);
	pmat[0][2] = random_zz_pX(0);
	pmat[1][0] = random_zz_pX(0);
	pmat[1][1] = random_zz_pX(0);
	pmat[1][2] = random_zz_pX(0);

	cout << pmat << endl;

	cout << "Starting tests:" << endl;

	std::vector<long> degs(pmat.NumRows());
	row_degree(degs,pmat);
	cout << "row degs: ";
	printVec(degs);

	degs.resize(pmat.NumCols());
	col_degree(degs,pmat);
	cout << "col degs: ";
	printVec(degs);

	Mat<long> deg_mat;
	degree_matrix(deg_mat, pmat);
	cout << "deg mat: " << endl << deg_mat << endl;

	Mat<zz_p> lead_mat;
	leading_matrix(lead_mat, pmat);
	cout << "leading mat: " << endl << lead_mat << endl;

	cout << "is reduced: " << boolalpha << is_reduced(pmat) << endl;

	std::vector<long> pivot(pmat.NumRows());
	pivot_index(pivot, pmat, std::vector<long>(), true);
	cout << "row pivot: ";
	printVec(pivot);

	pivot.resize(pmat.NumCols());
	pivot_index(pivot, pmat, std::vector<long>(), false);
	cout << "col pivot: ";
	printVec(pivot);

	cout << endl << "Tests for shifts: " << endl;

	std::vector<long> rs {1,2,3};
	cout << "row shift: ";
	printVec(rs);
	std::vector<long> cs {4,2};
	cout << "col shift: ";
	printVec(cs);

	degree_matrix(deg_mat,pmat,rs,true);
	cout << "degree matrix with row wise shift: " << endl << deg_mat<<endl;
	degree_matrix(deg_mat,pmat,cs,false);
	cout << "degree matrix with col wise shift: " << endl << deg_mat<<endl;

	degs.resize(pmat.NumRows());
	row_degree(degs,pmat,rs);
	cout << "shifted row degree: ";
	printVec(degs);

	degs.resize(pmat.NumCols());
	col_degree(degs,pmat,cs);
	cout << "shifted col degree: ";
	printVec(degs);

	leading_matrix(lead_mat,pmat,rs,true);
	cout << "row shifted leading mat: " << endl << lead_mat << endl;

	leading_matrix(lead_mat,pmat,cs,false);
	cout << "col shifted leading mat: " << endl << lead_mat << endl;

	pivot.resize(pmat.NumRows());
	pivot_index(pivot, pmat, rs, true);
	cout << "row pivot: ";
	printVec(pivot);

	pivot.resize(pmat.NumCols());
	pivot_index(pivot, pmat, cs, false);
	cout << "col pivot: ";
	printVec(pivot);
}
