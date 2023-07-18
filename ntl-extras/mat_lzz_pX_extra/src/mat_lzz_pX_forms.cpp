#include "lzz_pX_extra.h" // for "is_monic(a)"
#include "mat_lzz_pX_forms.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) ROW/COLUMN DEGREE                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* row degree                                                 */
/*------------------------------------------------------------*/
void row_degree(
                VecLong & rdeg,
                const Mat<zz_pX> & pmat
               )
{ 
    // empty rdeg and reserve space
    rdeg.clear();
    rdeg.reserve(pmat.NumRows());

    // take the max degree in each row of pmat
    for (long i = 0; i < pmat.NumRows(); ++i)
    {
        rdeg.emplace_back(-1);  // rdeg[i] == -1
        for (long j = 0; j < pmat.NumCols(); ++j)
            if (rdeg[i] < deg(pmat[i][j])) 
                rdeg[i] = deg(pmat[i][j]);
    }
}

/*------------------------------------------------------------*/
/* shifted row degree                                         */
/*------------------------------------------------------------*/
void row_degree(
                VecLong & rdeg,
                const Mat<zz_pX> &pmat,
                const VecLong & shift
               )
{
    if ((long)shift.size() != pmat.NumCols())
        throw std::invalid_argument("==row_degree== shift must have length pmat.NumCols()");

    // empty rdeg and reserve space
    rdeg.clear();
    rdeg.reserve(pmat.NumRows());

    // compute minimum shift entry (used for zero entries)
    long min_shift = *std::min_element(shift.begin(),shift.end());

    // take the max shifted degree in each row of pmat
    for (long i = 0; i < pmat.NumRows(); ++i)
    {
        rdeg.emplace_back(min_shift-1);  // rdeg[i] == min(shift)-1
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            long d = deg(pmat[i][j]) + shift[j];
            if (not IsZero(pmat[i][j]) && rdeg[i] < d)
                rdeg[i] = d;
        }
    }
}


/*------------------------------------------------------------*/
/* column degree                                              */
/*------------------------------------------------------------*/
void col_degree(
                VecLong &cdeg,
                const Mat<zz_pX> &pmat
               )
{
    // empty cdeg and reserve space
    cdeg.clear();
    cdeg.reserve(pmat.NumCols());

    // take the max degree in each column of pmat
    for (long j = 0; j < pmat.NumCols(); ++j)
    {
        cdeg.emplace_back(-1);  // cdeg[j] == -1
        for (long i = 0; i < pmat.NumRows(); ++i)
            if (cdeg[j] < deg(pmat[i][j])) 
                cdeg[j] = deg(pmat[i][j]);
    }
} 

/*------------------------------------------------------------*/
/* shifted column degree                                      */
/*------------------------------------------------------------*/
void col_degree(
                VecLong & cdeg,
                const Mat<zz_pX> &pmat,
                const VecLong & shift
               )
{
    if ((long)shift.size() != pmat.NumRows())
        throw std::invalid_argument("==col_degree== shift must have length pmat.NumRows()");

    // empty cdeg and reserve space
    cdeg.clear();
    cdeg.reserve(pmat.NumCols());

    // compute minimum shift entry (used for zero entries)
    long min_shift = *std::min_element(shift.begin(),shift.end());

    // take the max degree in each column of pmat
    for (long j = 0; j < pmat.NumCols(); ++j)
    {
        cdeg.emplace_back(min_shift-1);  // cdeg[j] == min(shift)-1
        for (long i = 0; i < pmat.NumRows(); ++i)
        {
            long d = deg(pmat[i][j]) + shift[i];
            if (not IsZero(pmat[i][j]) && cdeg[j] < d) 
                cdeg[j] = d;
        }
    }
} 




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT INDEX / DEGREE                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* row pivot index/degree                                     */
/*------------------------------------------------------------*/
void row_pivots(
                VecLong & pivind,
                VecLong & pivdeg,
                const Mat<zz_pX> & pmat
               )
{
    // make pivind/pivdeg of the right length, filled with -1
    pivind.clear();
    pivind.resize(pmat.NumRows(), -1);
    pivdeg.clear();
    pivdeg.resize(pmat.NumRows(), -1);

    for (long i = 0; i < pmat.NumRows(); ++i)
    {
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            long d = deg(pmat[i][j]);
            if (d > pivdeg[i])
            {
                pivdeg[i] = d;
                pivind[i] = j;
            }
            else if (d == pivdeg[i] && d>=0)
                pivind[i] = j;
        }
    }
}

/*------------------------------------------------------------*/
/* shifted row pivot index/degree                             */
/*------------------------------------------------------------*/
void row_pivots(
                VecLong & pivind,
                VecLong & pivdeg,
                const Mat<zz_pX> & pmat,
                const VecLong & shift
               )
{
    if ((long)shift.size() != pmat.NumCols())
        throw std::invalid_argument("==row_pivots== shift must have length pmat.NumCols()");

    // compute minimum shift entry (used for zero entries)
    long min_rdeg = *std::min_element(shift.begin(),shift.end()) - 1;

    // empty vectors and reserve space
    pivind.clear();
    pivind.resize(pmat.NumRows(), -1);
    pivdeg.clear();
    pivdeg.resize(pmat.NumRows(), min_rdeg);
    // for the moment, pivdeg will store the shifted row degree

    for (long i = 0; i < pmat.NumRows(); ++i)
    {
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            long d = deg(pmat[i][j]) + shift[j];
            if (d >= shift[j]) // pmat[i][j] nonzero
            {
                if (d > pivdeg[i])
                {
                    pivdeg[i] = d;
                    pivind[i] = j;
                }
                else if (d == pivdeg[i])
                    pivind[i] = j;
            }
        }
        // remove the shift so that pivdeg is not rdeg, but pivdeg
        if (pivind[i]==-1) // zero row
            pivdeg[i] = -1;
        else
            pivdeg[i] = pivdeg[i] - shift[pivind[i]];
    }
}

/*------------------------------------------------------------*/
/* column pivot index/degree                                  */
/*------------------------------------------------------------*/
void col_pivots(
                VecLong & pivind,
                VecLong & pivdeg,
                const Mat<zz_pX> & pmat
               )
{
    // make pivind/pivdeg of the right length, filled with -1
    pivind.clear();
    pivind.resize(pmat.NumCols(), -1);
    pivdeg.clear();
    pivdeg.resize(pmat.NumCols(), -1);

    for (long j = 0; j < pmat.NumCols(); ++j)
    {
        for (long i = 0; i < pmat.NumRows(); ++i)
        {
            long d = deg(pmat[i][j]);
            if (d > pivdeg[j])  // hence pmat[i][j] nonzero
            {
                pivdeg[j] = d;
                pivind[j] = i;
            }
            else if (d == pivdeg[j] && d >= 0)
                pivind[j] = i;
        }
    }
} 

/*------------------------------------------------------------*/
/* shifted column pivot index/degree                          */
/*------------------------------------------------------------*/
void col_pivots(
                VecLong & pivind,
                VecLong & pivdeg,
                const Mat<zz_pX> & pmat,
                const VecLong & shift
               )
{
    if ((long)shift.size() != pmat.NumRows())
        throw std::invalid_argument("==col_pivots== shift must have length pmat.NumRows()");

    // compute minimum shift entry (used for zero entries)
    long min_cdeg = *std::min_element(shift.begin(),shift.end()) - 1;

    // empty vectors and reserve space
    pivind.clear();
    pivind.resize(pmat.NumCols(), -1);
    pivdeg.clear();
    pivdeg.resize(pmat.NumCols(), min_cdeg);
    // for the moment, pivdeg will store the shifted column degree

    for (long j = 0; j < pmat.NumCols(); ++j)
    {
        for (long i = 0; i < pmat.NumRows(); ++i)
        {
            long d = deg(pmat[i][j]) + shift[i];
            if (d >= shift[i]) // pmat[i][j] nonzero
            {
                if (d > pivdeg[j])
                {
                    pivdeg[j] = d;
                    pivind[j] = i;
                }
                else if (d == pivdeg[j])
                    pivind[j] = i;
            }
        }
        // remove the shift so that pivdeg is not cdeg but pivdeg
        if (pivind[j]==-1) // zero column
            pivdeg[j] = -1;
        else
            pivdeg[j] = pivdeg[j] - shift[pivind[j]];
    }
}






/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) DEGREE MATRIX                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* degree matrix                                              */
/*------------------------------------------------------------*/
void degree_matrix(
                   Mat<long> & degmat,
                   const Mat<zz_pX> & pmat
                  )
{
    // set the dimensions of degmat and populate it with degrees
    degmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
            degmat[i][j] = deg(pmat[i][j]);
}

/*------------------------------------------------------------*/
/* row-shifted degree matrix                                  */
/*------------------------------------------------------------*/
void degree_matrix_rowshifted(
                              Mat<long> & degmat,
                              const Mat<zz_pX> & pmat,
                              const VecLong & shift
                             )
{
    if ((long)shift.size() != pmat.NumCols())
        throw std::invalid_argument("==degree_matrix_rowshifted== shift must have length pmat.NumCols()");

    // compute minimum shift entry (used for zero entries)
    long min_shift = *std::min_element(shift.begin(),shift.end());

    // set the dimensions of degmat and populate it with degrees
    degmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            if (IsZero(pmat[i][j]))
                degmat[i][j] = min_shift-1;
            else
                degmat[i][j] = deg(pmat[i][j]) + shift[j];
        }
}

/*------------------------------------------------------------*/
/* column-shifted degree matrix                               */
/*------------------------------------------------------------*/
void degree_matrix_colshifted(
                              Mat<long> & degmat,
                              const Mat<zz_pX> & pmat,
                              const VecLong & shift
                             )
{
    if ((long)shift.size() != pmat.NumRows())
        throw std::invalid_argument("==degree_matrix_colshifted== shift must have length pmat.NumRows()");

    // compute minimum shift entry (used for zero entries)
    long min_shift = *std::min_element(shift.begin(),shift.end());

    // set the dimensions of degmat and populate it with degrees
    degmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            if (IsZero(pmat[i][j]))
                degmat[i][j] = min_shift-1;
            else
                degmat[i][j] = deg(pmat[i][j]) + shift[i];
        }
}





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) LEADING MATRIX                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* row-wise leading matrix of pmat                            */
/*------------------------------------------------------------*/
void row_leading_matrix(
                        Mat<zz_p> & lmat,
                        const Mat<zz_pX> & pmat
                       )
{
    VecLong rdeg;
    row_degree(rdeg,pmat);

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (not IsZero(pmat[r][c]) && deg(pmat[r][c]) == rdeg[r])
                lmat[r][c] = pmat[r][c][deg(pmat[r][c])];
            else
                clear(lmat[r][c]);
}

/*------------------------------------------------------------*/
/* row-wise shifted leading matrix of pmat                    */
/*------------------------------------------------------------*/
void row_leading_matrix(
                        Mat<zz_p> & lmat,
                        const Mat<zz_pX> & pmat,
                        const VecLong & shift
                       )
{
    VecLong rdeg;
    row_degree(rdeg,pmat,shift); // throws if shift doesn't have the right length

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (not IsZero(pmat[r][c]) && deg(pmat[r][c])+shift[c] == rdeg[r])
                lmat[r][c] = pmat[r][c][deg(pmat[r][c])];
            else
                clear(lmat[r][c]);
}

/*------------------------------------------------------------*/
/* column-wise leading matrix of pmat                         */
/*------------------------------------------------------------*/
void col_leading_matrix(
                        Mat<zz_p> & lmat,
                        const Mat<zz_pX> & pmat
                       )
{
    VecLong cdeg;
    col_degree(cdeg,pmat);

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (not IsZero(pmat[r][c]) && deg(pmat[r][c]) == cdeg[c])
                lmat[r][c] = pmat[r][c][deg(pmat[r][c])];
            else
                clear(lmat[r][c]);
}

/*------------------------------------------------------------*/
/* column-wise shifted leading matrix of pmat                 */
/*------------------------------------------------------------*/
void col_leading_matrix(
                        Mat<zz_p> &lmat,
                        const Mat<zz_pX> &pmat,
                        const VecLong & shift
                       )
{
    VecLong cdeg;
    col_degree(cdeg,pmat,shift); // throws if shift doesn't have the right length

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (not IsZero(pmat[r][c]) && deg(pmat[r][c])+shift[r] == cdeg[c])
                lmat[r][c] = pmat[r][c][deg(pmat[r][c])];
            else
                clear(lmat[r][c]);
}




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED REDUCED FORM                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* test row reduced                                           */
/*------------------------------------------------------------*/
bool is_row_reduced(const Mat<zz_pX> & pmat)
{
    Mat<zz_p> lmat;
    row_leading_matrix(lmat,pmat);
    long rank = gauss(lmat);
    return rank == pmat.NumRows();
}

/*------------------------------------------------------------*/
/* test shifted row reduced                                   */
/*------------------------------------------------------------*/
bool is_row_reduced(const Mat<zz_pX> & pmat, const VecLong & shift)
{
    Mat<zz_p> lmat;
    row_leading_matrix(lmat,pmat,shift);
    long rank = gauss(lmat);
    return rank == pmat.NumRows();
}

/*------------------------------------------------------------*/
/* test column reduced                                        */
/*------------------------------------------------------------*/
bool is_col_reduced(const Mat<zz_pX> & pmat)
{
    Mat<zz_p> lmat;
    col_leading_matrix(lmat,pmat);
    long rank = gauss(lmat);
    return rank == pmat.NumCols();
}

/*------------------------------------------------------------*/
/* test shifted column reduced                                */
/*------------------------------------------------------------*/
bool is_col_reduced(const Mat<zz_pX> & pmat, const VecLong & shift)
{
    Mat<zz_p> lmat;
    col_leading_matrix(lmat,pmat,shift);
    long rank = gauss(lmat);
    return rank == pmat.NumCols();
}




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED WEAK POPOV FORM                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* test row-wise weak Popov form                              */
/*------------------------------------------------------------*/
bool is_row_weak_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for pairwise distinct: sort and check no adjacent equal elements
    std::sort(pivind.begin(),pivind.end());
    if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test row-wise shifted weak Popov form                      */
/*------------------------------------------------------------*/
bool is_row_weak_popov(
                       const Mat<zz_pX> &pmat,
                       const VecLong & shift
                      )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for pairwise distinct: sort and check no adjacent equal elements
    std::sort(pivind.begin(),pivind.end());
    if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test column-wise weak Popov form                           */
/*------------------------------------------------------------*/
bool is_col_weak_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for pairwise distinct: sort and check no adjacent equal elements
    std::sort(pivind.begin(),pivind.end());
    if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test column-wise shifted weak Popov form                   */
/*------------------------------------------------------------*/
bool is_col_weak_popov(
                       const Mat<zz_pX> &pmat,
                       const VecLong & shift
                      )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for pairwise distinct: sort and check no adjacent equal elements
    std::sort(pivind.begin(),pivind.end());
    if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test row-wise ordered weak Popov form                      */
/*------------------------------------------------------------*/
bool is_row_ordered_weak_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test row-wise shifted ordered weak Popov form              */
/*------------------------------------------------------------*/
bool is_row_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const VecLong & shift
                              )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test column-wise ordered weak Popov form                   */
/*------------------------------------------------------------*/
bool is_col_ordered_weak_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test column-wise shifted ordered weak Popov form           */
/*------------------------------------------------------------*/
bool is_col_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const VecLong & shift
                              )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    return true;
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED POPOV FORM                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* test row-wise Popov form                                   */
/*------------------------------------------------------------*/
bool is_row_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check pivot index increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    // --> ordered weak Popov form OK
    // now test pivot entries monic
    for (long i = 0; i < pmat.NumRows(); ++i)
        if (not is_monic(pmat[i][pivind[i]]))
            return false;

    // finally test degree of non-pivot entries
    for (long j = 0; j < pmat.NumRows(); ++j)
    {
        for (long i = 0; i < j; ++i)
            if (deg(pmat[i][pivind[j]]) >= pivdeg[j])
                return false;
        for (long i = j+1; i < pmat.NumRows(); ++i)
            if (deg(pmat[i][pivind[j]]) >= pivdeg[j])
                return false;
    }

    return true;
}

/*------------------------------------------------------------*/
/* test row-wise shifted Popov form                           */
/*------------------------------------------------------------*/
bool is_row_popov(const Mat<zz_pX> & pmat, const VecLong & shift)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check pivot index increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    // --> ordered weak Popov form OK
    // now test pivot entries monic
    for (long i = 0; i < pmat.NumRows(); ++i)
        if (not is_monic(pmat[i][pivind[i]]))
            return false;

    // finally test degree of non-pivot entries
    for (long j = 0; j < pmat.NumRows(); ++j)
    {
        for (long i = 0; i < j; ++i)
            if (deg(pmat[i][pivind[j]]) >= pivdeg[j])
                return false;
        for (long i = j+1; i < pmat.NumRows(); ++i)
            if (deg(pmat[i][pivind[j]]) >= pivdeg[j])
                return false;
    }

    return true;
}

/////*------------------------------------------------------------*/
/////* test row-wise Popov form, up to row permutation            */
/////*------------------------------------------------------------*/
////bool is_row_popov_up_to_permutation(const Mat<zz_pX> & pmat);
////
/////*------------------------------------------------------------*/
/////* test row-wise shifted Popov form, up to row permutation    */
/////*------------------------------------------------------------*/
////bool is_row_popov_up_to_permutation(
////                                    const Mat<zz_pX> & pmat,
////                                    const VecLong & shift
////                                    );

/*------------------------------------------------------------*/
/* test column-wise Popov form                                */
/*------------------------------------------------------------*/
bool is_col_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check pivot index increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    // --> ordered weak Popov form OK
    // now test pivot entries monic
    for (long j = 0; j < pmat.NumCols(); ++j)
        if (not is_monic(pmat[pivind[j]][j]))
            return false;

    // finally test degree of non-pivot entries
    for (long i = 0; i < pmat.NumCols(); ++i)
    {
        for (long j = 0; j < i; ++j)
            if (deg(pmat[pivind[i]][j]) >= pivdeg[i])
                return false;
        for (long j = i+1; j < pmat.NumRows(); ++j)
            if (deg(pmat[pivind[i]][j]) >= pivdeg[i])
                return false;
    }

    return true;
}


/*------------------------------------------------------------*/
/* test column-wise shifted Popov form                        */
/*------------------------------------------------------------*/
bool is_col_popov(const Mat<zz_pX> & pmat, const VecLong & shift)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if (std::find(pivind.begin(),pivind.end(),-1) != pivind.end())
        return false;

    // check pivot index increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivind.begin(),pivind.end(),std::greater_equal<long>()) != pivind.end())
        return false;

    // --> ordered weak Popov form OK
    // now test pivot entries monic
    for (long j = 0; j < pmat.NumCols(); ++j)
        if (not is_monic(pmat[pivind[j]][j]))
            return false;

    // finally test degree of non-pivot entries
    for (long i = 0; i < pmat.NumCols(); ++i)
    {
        for (long j = 0; j < i; ++j)
            if (deg(pmat[pivind[i]][j]) >= pivdeg[i])
                return false;
        for (long j = i+1; j < pmat.NumRows(); ++j)
            if (deg(pmat[pivind[i]][j]) >= pivdeg[i])
                return false;
    }

    return true;
}

////*------------------------------------------------------------*/
////* test column-wise Popov form, up to column permutation      */
////*------------------------------------------------------------*/
///bool is_col_popov_up_to_permutation(const Mat<zz_pX> & pmat);
///
////*------------------------------------------------------------*/
////* test column-wise shifted Popov form,                       */
////* up to column permutation                                   */
////*------------------------------------------------------------*/
///bool is_col_popov_up_to_permutation(
///                                    const Mat<zz_pX> & pmat,
///                                    const VecLong & shift
///                                    );





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED FORMS (FORM SELECTOR)                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* Check whether pmat is in the prescribed row-wise form      */
/*------------------------------------------------------------*/

bool is_row_polmatform(const PolMatForm form, const Mat<zz_pX> & pmat)
{
    switch (form)
    {
    case NONE:
        return true;
    case REDUCED:
        return is_row_reduced(pmat);
    case WEAK_POPOV:
        return is_row_weak_popov(pmat);
    case ORD_WEAK_POPOV:
        return is_row_ordered_weak_popov(pmat);
    case POPOV:
        return is_row_popov(pmat);
    default:
        throw std::invalid_argument("==is_row_polmatform== Unknown required polynomial matrix form.");
    }
}

bool is_row_polmatform(
                       const PolMatForm form,
                       const Mat<zz_pX> & pmat,
                       const VecLong &shift
                      )
{
    switch (form)
    {
    case NONE:
        return true;
    case REDUCED:
        return is_row_reduced(pmat,shift);
    case WEAK_POPOV:
        return is_row_weak_popov(pmat,shift);
    case ORD_WEAK_POPOV:
        return is_row_ordered_weak_popov(pmat,shift);
    case POPOV:
        return is_row_popov(pmat,shift);
    default:
        throw std::invalid_argument("==is_row_polmatform== Unknown required polynomial matrix form.");
    }
}

/*------------------------------------------------------------*/
/* Check whether pmat is in the prescribed column-wise form   */
/*------------------------------------------------------------*/

bool is_col_polmatform(const PolMatForm form, const Mat<zz_pX> & pmat)
{
    switch (form)
    {
    case NONE:
        return true;
    case REDUCED:
        return is_col_reduced(pmat);
    case WEAK_POPOV:
        return is_col_weak_popov(pmat);
    case ORD_WEAK_POPOV:
        return is_col_ordered_weak_popov(pmat);
    case POPOV:
        return is_col_popov(pmat);
    default:
        throw std::invalid_argument("==is_col_polmatform== Unknown required polynomial matrix form.");
    }
}

bool is_col_polmatform(
                       const PolMatForm form,
                       const Mat<zz_pX> & pmat,
                       const VecLong &shift
                      )
{
    switch (form)
    {
    case NONE:
        return true;
    case REDUCED:
        return is_col_reduced(pmat,shift);
    case WEAK_POPOV:
        return is_col_weak_popov(pmat,shift);
    case ORD_WEAK_POPOV:
        return is_col_ordered_weak_popov(pmat,shift);
    case POPOV:
        return is_col_popov(pmat,shift);
    default:
        throw std::invalid_argument("==is_col_polmatform== Unknown required polynomial matrix form.");
    }
}






/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* FIND STRONGEST (SHIFTED) FORM                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* Return the strongest row-wise form of pmat among:          */
/* Popov => ordered weak Popov => weak Popov                  */
/*                                       => reduced => none   */
/*------------------------------------------------------------*/
PolMatForm get_row_polmatform(const Mat<zz_pX> & pmat)
{
    if (is_row_popov(pmat))
        return POPOV;
    if (is_row_ordered_weak_popov(pmat))
        return ORD_WEAK_POPOV;
    if (is_row_weak_popov(pmat))
        return WEAK_POPOV;
    if (is_row_reduced(pmat))
        return REDUCED;
    return NONE;
}

PolMatForm get_row_polmatform(
                              const Mat<zz_pX> &pmat,
                              const VecLong &shift
                             )
{
    if (is_row_popov(pmat,shift))
        return POPOV;
    if (is_row_ordered_weak_popov(pmat,shift))
        return ORD_WEAK_POPOV;
    if (is_row_weak_popov(pmat,shift))
        return WEAK_POPOV;
    if (is_row_reduced(pmat,shift))
        return REDUCED;
    return NONE;
}

/*------------------------------------------------------------*/
/* Return the strongest column-wise form of pmat among:       */
/* Popov => ordered weak Popov => weak Popov                  */
/*                                       => reduced => none   */
/*------------------------------------------------------------*/
PolMatForm get_col_polmatform(const Mat<zz_pX> & pmat)
{
    if (is_col_popov(pmat))
        return POPOV;
    if (is_col_ordered_weak_popov(pmat))
        return ORD_WEAK_POPOV;
    if (is_col_weak_popov(pmat))
        return WEAK_POPOV;
    if (is_col_reduced(pmat))
        return REDUCED;
    return NONE;
}

PolMatForm get_col_polmatform(
                              const Mat<zz_pX> &pmat,
                              const VecLong &shift
                             )
{
    if (is_col_popov(pmat,shift))
        return POPOV;
    if (is_col_ordered_weak_popov(pmat,shift))
        return ORD_WEAK_POPOV;
    if (is_col_weak_popov(pmat,shift))
        return WEAK_POPOV;
    if (is_col_reduced(pmat,shift))
        return REDUCED;
    return NONE;
}



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
