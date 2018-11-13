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
                DegVec & rdeg,
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
                        DegVec & rdeg,
                        const Mat<zz_pX> &pmat,
                        const Shift & shift
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
                DegVec &cdeg,
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
                        DegVec & cdeg,
                        const Mat<zz_pX> &pmat,
                        const Shift & shift
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
    // empty vectors and reserve space
    pivind.clear();
    pivind.reserve(pmat.NumRows());
    pivdeg.clear();
    pivdeg.reserve(pmat.NumRows());

    for (long i = 0; i < pmat.NumRows(); ++i)
    {
        pivind.emplace_back(-1);  // pivind[i] == -1
        pivdeg.emplace_back(-1);  // pivdeg[i] == -1
        for (long j = 0; j < pmat.NumCols(); ++j)
            if (deg(pmat[i][j]) > pivdeg[i])  // hence pmat[i][j] nonzero
            {
                pivdeg[i] = deg(pmat[i][j]);
                pivind[i] = j;
            }
            else if (deg(pmat[i][j]) == pivdeg[i] && not IsZero(pmat[i][j]))
                pivind[i] = j;
    }
}

/*------------------------------------------------------------*/
/* shifted row pivot index/degree                             */
/*------------------------------------------------------------*/
void row_pivots(
                        VecLong & pivind,
                        VecLong & pivdeg,
                        const Mat<zz_pX> & pmat,
                        const Shift & shift
                       )
{
    if ((long)shift.size() != pmat.NumCols())
        throw std::invalid_argument("==row_pivots== shift must have length pmat.NumCols()");

    // empty vectors and reserve space
    pivind.clear();
    pivind.reserve(pmat.NumRows());
    pivdeg.clear();
    pivdeg.reserve(pmat.NumRows());

    // compute minimum shift entry (used for zero entries)
    long min_shift = *std::min_element(shift.begin(),shift.end());

    // take the max degree in each row of pmat
    for (long i = 0; i < pmat.NumRows(); ++i)
    {
        pivind.emplace_back(-1);  // pivind[i] == -1
        // for the moment, pivdeg stores the rdeg, pivdeg[i] == min(shift)-1
        pivdeg.emplace_back(min_shift-1);
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            long d = deg(pmat[i][j]) + shift[j];
            if (d > pivdeg[i]) // hence pmat[i][j] nonzero
            {
                pivdeg[i] = d;
                pivind[i] = j;
            }
            if (d == pivdeg[i] && not IsZero(pmat[i][j]))
                pivind[i] = j;
        }
        // remove the shift so that pivdeg is not rdeg but pivdeg
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
    // empty pivind / pivdeg and reserve space
    pivind.clear();
    pivind.reserve(pmat.NumCols());
    pivdeg.clear();
    pivdeg.reserve(pmat.NumCols());

    for (long j = 0; j < pmat.NumCols(); ++j)
    {
        pivind.emplace_back(-1);  // pivind[j] == -1
        pivdeg.emplace_back(-1);  // pivdeg[j] == -1
        for (long i = 0; i < pmat.NumRows(); ++i)
            if (deg(pmat[i][j]) > pivdeg[j])  // hence pmat[i][j] nonzero
            {
                pivdeg[j] = deg(pmat[i][j]);
                pivind[j] = i;
            }
            else if (deg(pmat[i][j]) == pivdeg[j] && not IsZero(pmat[i][j]))
                pivind[j] = i;
    }
} 

/*------------------------------------------------------------*/
/* shifted column pivot index/degree                          */
/*------------------------------------------------------------*/
void col_pivots(
                        VecLong & pivind,
                        VecLong & pivdeg,
                        const Mat<zz_pX> & pmat,
                        const Shift & shift
                       )
{
    if ((long)shift.size() != pmat.NumRows())
        throw std::invalid_argument("==col_pivots== shift must have length pmat.NumRows()");

    // empty vectors and reserve space
    pivind.clear();
    pivind.reserve(pmat.NumCols());
    pivdeg.clear();
    pivdeg.reserve(pmat.NumCols());

    // compute minimum shift entry (used for zero entries)
    long min_shift = *std::min_element(shift.begin(),shift.end());

    // take the max degree in each row of pmat
    for (long j = 0; j < pmat.NumCols(); ++j)
    {
        pivind.emplace_back(-1);  // pivind[j] == -1
        // for the moment, pivdeg stores the cdeg, pivdeg[j] == min(shift)-1
        pivdeg.emplace_back(min_shift-1);
        for (long i = 0; i < pmat.NumRows(); ++i)
        {
            long d = deg(pmat[i][j]) + shift[i];
            if (d > pivdeg[j]) // hence pmat[i][j] nonzero
            {
                pivdeg[j] = d;
                pivind[j] = i;
            }
            if (d == pivdeg[j] && not IsZero(pmat[i][j]))
                pivind[j] = i;
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
                              const Shift & shift
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
                              const Shift & shift
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
    DegVec rdeg;
    row_degree(rdeg,pmat);

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (deg(pmat[r][c]) == rdeg[r])
                lmat[r][c] = pmat[r][c][deg(pmat[r][c])];
            else
                clear(lmat[r][c]);
}

/*------------------------------------------------------------*/
/* row-wise shifted leading matrix of pmat                    */
/*------------------------------------------------------------*/
void row_leading_matrix(
                                Mat<zz_p> &lmat,
                                const Mat<zz_pX> &pmat,
                                const Shift & shift
                               )
{
    DegVec rdeg;
    row_degree(rdeg,pmat,shift);

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (deg(pmat[r][c])+shift[c] == rdeg[r])
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
    DegVec cdeg;
    col_degree(cdeg,pmat);

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (deg(pmat[r][c]) == cdeg[c])
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
                                const Shift & shift
                               )
{
    DegVec cdeg;
    col_degree(cdeg,pmat,shift);

    lmat.SetDims(pmat.NumRows(), pmat.NumCols());
    for (long r = 0; r < lmat.NumRows(); ++r)
        for (long c = 0; c < lmat.NumCols(); ++c)
            if (deg(pmat[r][c])+shift[r] == cdeg[c])
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
bool is_row_reduced(const Mat<zz_pX> & pmat, const Shift & shift)
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
bool is_col_reduced(const Mat<zz_pX> & pmat, const Shift & shift)
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
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
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
                               const Shift & shift
                              )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
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
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
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
                               const Shift & shift
                              )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
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
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivots.begin(),pivots.end(),std::greater_equal<long>()) != pivots.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test row-wise shifted ordered weak Popov form              */
/*------------------------------------------------------------*/
bool is_row_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const Shift & shift
                              )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    row_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivots.begin(),pivots.end(),std::greater_equal<long>()) != pivots.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test column-wise ordered weak Popov form                   */
/*------------------------------------------------------------*/
bool is_row_ordered_weak_popov(const Mat<zz_pX> & pmat)
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat);

    // forbide zero rows
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivots.begin(),pivots.end(),std::greater_equal<long>()) != pivots.end())
        return false;

    return true;
}

/*------------------------------------------------------------*/
/* test column-wise shifted ordered weak Popov form           */
/*------------------------------------------------------------*/
bool is_row_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const Shift & shift
                              )
{
    // retrieve pivot index
    VecLong pivind, pivdeg;
    col_pivots(pivind, pivdeg, pmat, shift);

    // forbide zero rows
    if( std::find(pivind.begin(),pivind.end(),-1) != pivind.end() )
        return false;

    // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
    if (std::adjacent_find(pivots.begin(),pivots.end(),std::greater_equal<long>()) != pivots.end())
        return false;

    return true;
}




bool is_monic(const zz_pX &p){
    return IsOne(LeadCoeff(p));
}

/*------------------------------------------------------------*/
/* TODO comment                                               */
/*------------------------------------------------------------*/
bool is_popov(
              const Mat<zz_pX> &pmat,
              const Shift &shift,
              const bool row_wise,
              const bool up_to_permutation
             )
{
    if (!is_weak_popov(pmat,shift,row_wise,!up_to_permutation))
        return false;


    std::vector<long> pivots;
    DegVec degrees;
    pivots.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
    degrees.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
    pivot_index(pivots,degrees,pmat,shift,row_wise);
    for (unsigned long i = 0; i < pivots.size(); i++){
        long index = pivots[i];
        if (index >= 0){
            if(row_wise){
                if (!is_monic(pmat[i][index]))
                    return false;
                for (long k=0; k < pmat.NumRows(); k++){
                    if (deg(pmat[k][index]) >= degrees[i] && ((unsigned long)k != i))
                        return false;
                }
            }else{ // col-wise
                if (!is_monic(pmat[index][i]))
                    return false;
                for (long k = 0; k < pmat.NumCols(); k++)
                    if (deg(pmat[index][k]) >= degrees[i] && ((unsigned long)k != i))
                        return false;
            }
        }
    }
    return true;
}



PolMatForm get_polmatform(
                          const Mat<zz_pX> &pmat,
                          const Shift &shift,
                          const bool row_wise
                         )
{
    if (is_popov(pmat,shift,row_wise)) {   // TODO waiting for is_popov
        return POPOV;
    }
    if (is_weak_popov(pmat,shift,row_wise,true))
        return ORD_WEAK_POPOV;
    else if (is_weak_popov(pmat,shift,row_wise))
        return WEAK_POPOV;
    else if (is_reduced(pmat,shift,row_wise))
        return REDUCED;
    else
        return NONE;
}

bool is_polmatform(
                   const Mat<zz_pX> &pmat,
                   const PolMatForm form,
                   const Shift &shift,
                   const bool row_wise
                  )
{
    switch (form)
    {
    case NONE: return true;
    case REDUCED: return is_reduced(pmat,shift,row_wise);
    case WEAK_POPOV: return is_weak_popov(pmat,shift,row_wise,false);
    case ORD_WEAK_POPOV: return is_weak_popov(pmat,shift,row_wise,true);
    case POPOV: return is_popov(pmat,shift,row_wise,false);
    default: throw std::invalid_argument("==is_polmatform== Unknown required polynomial matrix form.");
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
