#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>

#include "mat_lzz_pX_extra.h"

NTL_CLIENT

template<typename T>
void printVec(std::vector<T> const &input)
{
    std::cout << '[';
    for (auto &elt : input) {
        std::cout << elt << " ";
    }
    std::cout << ']' << std::endl;
}

/******************************************************/
/* Tests the degree functions for polynomial matrices */
/******************************************************/


using namespace std;
int main()
{
    zz_p::init(13);

    Mat<zz_pX> pmat;
    pmat.SetDims(3,4);
    pmat[0][0] = random_zz_pX(4);
    pmat[0][1] = random_zz_pX(1);
    pmat[0][2] = random_zz_pX(2);
    pmat[0][3] = random_zz_pX(0);
    pmat[1][0] = random_zz_pX(3);
    pmat[1][1] = random_zz_pX(1);
    pmat[1][2] = random_zz_pX(1);
    pmat[1][3] = random_zz_pX(0);
    pmat[2][0] = random_zz_pX(4);
    pmat[2][1] = random_zz_pX(2);
    pmat[2][2] = random_zz_pX(3);
    pmat[2][3] = random_zz_pX(3);

    cout << pmat << endl;

    cout << "Starting tests:" << endl;

    VecLong degs(pmat.NumRows());
    row_degree(degs,pmat);
    cout << "row degs: ";
    printVec(degs);

    col_degree(degs,pmat);
    cout << "col degs: ";
    printVec(degs);

    Mat<long> deg_mat;
    degree_matrix(deg_mat, pmat);
    cout << "deg mat: " << endl << deg_mat << endl;

    Mat<zz_p> lead_mat;
    row_leading_matrix(lead_mat, pmat);
    cout << "leading mat: " << endl << lead_mat << endl;

    cout << "is reduced: " << boolalpha << is_row_reduced(pmat) << endl;

    VecLong pivind(pmat.NumRows());
    VecLong pivdeg(pmat.NumRows());
    row_pivots(pivind, pivdeg, pmat);
    cout << "row pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    col_pivots(pivind, pivdeg, pmat);
    cout << "col pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    cout << endl << "Tests for shifts: " << endl;

    VecLong rs {0,2,1,3};
    cout << "row shift: ";
    printVec(rs);
    VecLong cs {4,2,0};
    cout << "col shift: ";
    printVec(cs);

    degree_matrix_rowshifted(deg_mat,pmat,rs);
    cout << "degree matrix with row wise shift: " << endl << deg_mat<<endl;
    degree_matrix_colshifted(deg_mat,pmat,cs);
    cout << "degree matrix with col wise shift: " << endl << deg_mat<<endl;

    row_degree(degs,pmat,rs);
    cout << "shifted row degree: ";
    printVec(degs);

    col_degree(degs,pmat,cs);
    cout << "shifted col degree: ";
    printVec(degs);

    row_leading_matrix(lead_mat,pmat,rs);
    cout << "row shifted leading mat: " << endl << lead_mat << endl;

    col_leading_matrix(lead_mat,pmat,cs);
    cout << "col shifted leading mat: " << endl << lead_mat << endl;

    row_pivots(pivind, pivdeg, pmat, rs);
    cout << "row pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    col_pivots(pivind, pivdeg, pmat, cs);
    cout << "col pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    cout << "weak popov? " << is_row_weak_popov(pmat,rs) << endl;
    cout << "weak ordered popov? " << is_row_ordered_weak_popov(pmat,rs) << endl;

    cout << endl << "testing is_popov" << endl;
    pmat.SetDims(3,3);
    pmat[0][0] = random_zz_pX(5);
    pmat[0][1] = random_zz_pX(2);
    pmat[0][2] = random_zz_pX(1);
    pmat[1][0] = random_zz_pX(3);
    pmat[1][1] = random_zz_pX(3);
    pmat[1][2] = random_zz_pX(1);
    pmat[2][0] = random_zz_pX(2);
    pmat[2][1] = random_zz_pX(2);
    pmat[2][2] = random_zz_pX(0);
    MakeMonic(pmat[0][0]);
    MakeMonic(pmat[1][0]);
    MakeMonic(pmat[1][1]);
    MakeMonic(pmat[2][2]);
    cout << "pmat: " << endl << pmat << endl;
    degree_matrix(deg_mat, pmat);
    cout << "deg mat: " << endl << deg_mat << endl;
    cout << "is_popov: " << boolalpha << is_row_popov(pmat) << endl;

    cout << endl << "Test Left shifts" << endl;
    cout << "left shift operator: " << (pmat << 2) << endl;
    cout << "left shift mutator: " << (pmat <<= 2) << endl;
    cout << "test procedure: " << LeftShift(pmat,2) << endl;
    LeftShift(pmat,pmat,2);
    cout << "test mutator: " << pmat << endl;
    //row
    cout << "row test procedure: " << LeftShiftRow(pmat,0,1) << endl;
    LeftShiftRow(pmat,pmat,0,1);
    cout << "row test mutator: " << pmat << endl;
    //col
    cout << "col test procedure: " << LeftShiftCol(pmat,1,2) << endl;
    LeftShiftCol(pmat,pmat,1,2);
    cout << "col test mutator: " << pmat << endl;

    cout << endl << "Test Right shifts" << endl;
    cout << "left shift operator: " << (pmat >> 2) << endl;
    cout << "left shift mutator: " << (pmat >>= 2) << endl;
    cout << "test procedure: " << RightShift(pmat,2) << endl;
    RightShift(pmat,pmat,2);
    cout << "test mutator: " << pmat << endl;
    //row
    cout << "row test procedure: " << RightShiftRow(pmat,0,1) << endl;
    RightShiftRow(pmat,pmat,0,1);
    cout << "row test mutator: " << pmat << endl;
    //col
    cout << "col test procedure: " << RightShiftCol(pmat,1,2) << endl;
    RightShiftCol(pmat,pmat,1,2);
    cout << "col test mutator: " << pmat << endl;

    cout << endl << "Test trunc" << endl;
    cout << "trunc: " << trunc(pmat,4) << endl;
    trunc(pmat,pmat,4);
    cout << "mutator trunc: " << pmat << endl;

    cout << endl << "Test trunc row" << endl;
    cout << "trunc: " << truncRow(pmat,0,2) << endl;
    truncRow(pmat,pmat,0,2);
    cout << "mutator trunc: " << pmat << endl;

    cout << endl << "Test trunc col" << endl;
    cout << "trunc: " << truncCol(pmat,1,1) << endl;
    truncCol(pmat,pmat,1,1);
    cout << "mutator trunc: " << pmat << endl;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
