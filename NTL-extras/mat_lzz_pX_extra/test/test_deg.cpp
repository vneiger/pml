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
int main(){
    zz_p::init(13);

    Mat<zz_pX> pmat;
    pmat.SetDims(3,4);
    // TODO there is a random function for matrices :-)
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

    std::vector<long> degs(pmat.NumRows());
    row_degree(degs,pmat);
    cout << "row degs: ";
    printVec(degs);

    degs.resize(pmat.NumCols());
    column_degree(degs,pmat);
    cout << "col degs: ";
    printVec(degs);

    Mat<long> deg_mat;
    degree_matrix(deg_mat, pmat);
    cout << "deg mat: " << endl << deg_mat << endl;

    Mat<zz_p> lead_mat;
    leading_matrix(lead_mat, pmat);
    cout << "leading mat: " << endl << lead_mat << endl;

    cout << "is reduced: " << boolalpha << is_reduced(pmat) << endl;

    std::vector<long> pivind(pmat.NumRows());
    std::vector<long> pivdeg(pmat.NumRows());
    pivot_index(pivind, pivdeg, pmat, std::vector<long>(), true);
    cout << "row pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    pivind.resize(pmat.NumCols());
    pivdeg.resize(pmat.NumCols());
    pivot_index(pivind, pivdeg, pmat, std::vector<long>(), false);
    cout << "col pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    cout << endl << "Tests for shifts: " << endl;

    std::vector<long> rs {0,2,1,3};
    cout << "row shift: ";
    printVec(rs);
    std::vector<long> cs {4,2,0};
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
    column_degree(degs,pmat,cs);
    cout << "shifted col degree: ";
    printVec(degs);

    leading_matrix(lead_mat,pmat,rs,true);
    cout << "row shifted leading mat: " << endl << lead_mat << endl;

    leading_matrix(lead_mat,pmat,cs,false);
    cout << "col shifted leading mat: " << endl << lead_mat << endl;

    pivind.resize(pmat.NumRows());
    pivdeg.resize(pmat.NumRows());
    pivot_index(pivind, pivdeg, pmat, rs, true);
    cout << "row pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    pivind.resize(pmat.NumCols());
    pivdeg.resize(pmat.NumCols());
    pivot_index(pivind, pivdeg, pmat, cs, false);
    cout << "col pivot: ";
    printVec(pivind);
    printVec(pivdeg);

    cout << "weak popov? " << is_weak_popov(pmat,rs) << endl;
    cout << "weak ordered popov? " << is_weak_popov(pmat,rs,true,true) << endl;

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
    cout << "is_popov: " << boolalpha << is_popov(pmat) << endl;

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
