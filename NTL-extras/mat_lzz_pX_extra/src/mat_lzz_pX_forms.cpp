#include "mat_lzz_pX_forms.h"

NTL_CLIENT

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
/* shifted row degree of a matrix                             */
/*------------------------------------------------------------*/
void row_degree(
                DegVec & rdeg,
                const Mat<zz_pX> & pmat,
                const Shift & shift
               )
{ 
    rdeg.resize(pmat.NumRows());

    // retrieve the shifted degree matrix
    Mat<long> degmat;
    degree_matrix(degmat,pmat,shift,true);

    // take the max of each row of degmat
    for (long r = 0; r < pmat.NumRows(); ++r)
    {
        long max_deg = degmat[r][0];
        for (long c = 1; c < pmat.NumCols(); ++c)
            if (max_deg < degmat[r][c]) 
                max_deg = degmat[r][c];
        rdeg[r] = max_deg;
    }
}

/*------------------------------------------------------------*/
/* shifted column degree of a matrix                          */
/*------------------------------------------------------------*/
void column_degree(
                   DegVec &cdeg,
                   const Mat<zz_pX> &pmat,
                   const Shift & shift
                  )
{
    cdeg.resize(pmat.NumCols());

    // retrieve the shifted degree matrix
    Mat<long> degmat;
    degree_matrix(degmat,pmat,shift,false);

    // take the max of each column of degmat
    for (long c = 0; c < pmat.NumCols(); c++)
    {
        long max_deg = degmat[0][c];
        for (long r = 1; r < pmat.NumRows(); r++)
            if (max_deg < degmat[r][c]) 
                max_deg = degmat[r][c];
        cdeg[c] = max_deg;
    }
} 

/*------------------------------------------------------------*/
/* shifted leading matrix                                     */
/*------------------------------------------------------------*/
void leading_matrix(
                    Mat<zz_p> &lmat,
                    const Mat<zz_pX> &pmat,
                    const Shift & shift,
                    const bool row_wise
                   )
{
    // check if shifted + check shift dimension
    bool shifted = check_shift(shift, pmat, row_wise);

    if (row_wise)
    {
        DegVec degree;
        row_degree(degree,pmat,shift);

        lmat.SetDims(pmat.NumRows(), pmat.NumCols());
        for (long r = 0; r < pmat.NumRows(); r++)
        {
            for (long c = 0; c < pmat.NumCols(); c++)
            {
                long d = deg(pmat[r][c]);

                if (shifted)
                    d += shift[c];

                if (d == degree[r])
                    lmat[r][c] = LeadCoeff(pmat[r][c]);
            }
        }
    }
    else
    {
        DegVec degree;
        column_degree(degree,pmat,shift);

        lmat.SetDims(pmat.NumRows(), pmat.NumCols());
        for (long r = 0; r < pmat.NumRows(); r++)
        {
            for (long c = 0; c < pmat.NumCols(); c++)
            {
                long d = deg(pmat[r][c]);

                if (shifted)
                    d += shift[r];

                if (d == degree[c])
                    lmat[r][c] = LeadCoeff(pmat[r][c]);
            }
        }
    }
}

/*------------------------------------------------------------*/
/* test whether pmat is a reduced matrix                      */
/*------------------------------------------------------------*/
bool is_reduced(
                const Mat<zz_pX> & pmat,
                const Shift & shift,
                const bool row_wise
               )
{
    Mat<zz_p> lmat;
    leading_matrix(lmat,pmat,shift,row_wise);
    long rank = gauss(lmat);
    return rank == (row_wise ? pmat.NumRows() : pmat.NumCols());
}

/*------------------------------------------------------------*/
/* TODO comment                                               */
/*------------------------------------------------------------*/
void pivot_index(
                 std::vector<long> & pivind,
                 DegVec & pivdeg,
                 const Mat<zz_pX> & pmat,
                 const Shift & shift,
                 const bool row_wise
                )
{
    // check if shifted + shift dimension
    bool shifted = check_shift(shift, pmat, row_wise);

    // retrieve (shifted) degree matrix
    Mat<long> deg_mat;
    degree_matrix(deg_mat,pmat,shift,row_wise);

    DegVec degree;
    if (row_wise)
    {
        degree.resize(pmat.NumRows());
        row_degree(degree,pmat,shift);
    }
    else
    {
        degree.resize(pmat.NumCols());
        column_degree(degree,pmat,shift);
    }

    long zero_degree = -1;
    if (shifted)
        zero_degree = *std::min_element(shift.begin(),shift.end()) -1;

    if (row_wise)
    {
        if ((long)pivind.size() != pmat.NumRows())
            throw std::invalid_argument("==pivot_index== Provided vector does not have size = NumRows");

        for (long r = 0; r < pmat.NumRows(); ++r)
        {
            if (degree[r] == zero_degree) 
            {
                pivdeg[r] = -1;
                pivind[r] = -1;
            }
            else
            {
                for (long c = 0; c <pmat.NumCols(); ++c)
                {
                    if (deg_mat[r][c] == degree[r]) 
                    {
                        pivdeg[r] = deg(pmat[r][c]);
                        pivind[r] = c;
                    }
                }
            }
        }
    }
    else
    {
        if ((long)pivind.size() != pmat.NumCols())
            throw std::invalid_argument("==pivot_index== Provided vector does not have size = NumCols");

        for(long c = 0; c < pmat.NumCols(); c++)
        {
            if (degree[c] == zero_degree) 
            {
                pivdeg[c] = -1;
                pivind[c] = -1;
            }
            else
            {
                for (long r = 0; r < pmat.NumRows(); r++)
                {
                    if (deg_mat[r][c] == degree[c]) 
                    {
                        pivdeg[c] = deg(pmat[r][c]);
                        pivind[c] = r;
                    }
                }
            }
        }
    }
}

/*------------------------------------------------------------*/
/* TODO comment                                               */
/*------------------------------------------------------------*/
bool is_weak_popov (
                    const Mat<zz_pX> &pmat,
                    const Shift &shift,
                    const bool row_wise,
                    const bool ordered
                   )
{
    //retrieve pivot index
    std::vector<long> pivots;
    DegVec pivdeg;
    pivots.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
    pivdeg.resize(row_wise ? pmat.NumRows() : pmat.NumCols());
    pivot_index(pivots, pivdeg, pmat, shift, row_wise);

    // forbide zero vectors
    if( std::find(pivots.begin(),pivots.end(),-1) != pivots.end() )
        return false;

    if (!ordered)
    { // check for pairwise distinct: sort and check no adjacent equal elements
        std::sort(pivots.begin(),pivots.end());
        if (std::adjacent_find(pivots.begin(),pivots.end()) != pivots.end())
            return false;
    }
    else
    { // check for strictly increasing: no adjacent elements elt1,elt2 with elt1 >= elt2
        if (std::adjacent_find(pivots.begin(),pivots.end(),std::greater_equal<long>()) != pivots.end())
            return false;
    }
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
