#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/******************************************************/
/* Tests the degree functions for polynomial matrices */
/******************************************************/

int main(){
	zz_p::init(13);
	Mat<zz_pX> m;
	m.SetDims(2,3);
	m[0][0] = random_zz_pX(3);
	m[0][1] = random_zz_pX(2);
	m[0][2] = random_zz_pX(0);
	m[1][0] = random_zz_pX(0);
	m[1][1] = random_zz_pX(0);
	m[1][2] = random_zz_pX(0);

	cout << m << endl;
	
	cout << "Starting tests:" << endl;
	Vec<long> degs;
	row_degree(degs,m);
	cout << "row degs: " << degs << endl;
	col_degree(degs,m);
	cout << "col degs: " << degs << endl;
	Mat<long> deg_mat;
	degree_matrix(deg_mat, m);
	cout << "deg mat: " << endl << deg_mat << endl;
	Mat<zz_p> lead_mat;
	leading_matrix(lead_mat, m);
	cout << "leading mat: " << endl << lead_mat << endl;
	cout << "is reduced: " << boolalpha << is_reduced(m) << endl;
	Vec<long> pivot;
	pivot_index(pivot, m, Vec<long>(), true);
	cout << "row pivot: " << pivot << endl;
	pivot_index(pivot, m, Vec<long>(), false);
	cout << "col pivot: " << pivot << endl;
	
	cout << endl << "Tests for shifts: " << endl;
	Vec<long> rs;
	rs.append(1); rs.append(2); rs.append(3);
	cout << "row shift: " << rs << endl;
	Vec<long> cs;
	cs.append(4); cs.append(2);
	cout << "col shift: " << cs << endl;
	degree_matrix(deg_mat,m,rs,true);
	cout << "row wise degree after shift: " << endl << deg_mat<<endl;
	degree_matrix(deg_mat,m,cs,false);
	cout << "col wise degree after shift: " << endl << deg_mat<<endl;
	row_degree(degs,m,rs);
	cout << "shifted row degree: " << degs << endl;
	col_degree(degs,m,cs);
	cout << "shifted col degree: " << degs << endl;
	leading_matrix(lead_mat,m,rs,true);
	cout << "row shifted leading mat: " << endl << lead_mat << endl;
	leading_matrix(lead_mat,m,cs,false);
	cout << "col shifted leading mat: " << endl << lead_mat << endl;
	pivot_index(pivot, m, rs, true);
	cout << "row pivot: " << pivot << endl;
	pivot_index(pivot, m, cs, false);
	cout << "col pivot: " << pivot << endl;
}
