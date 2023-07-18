#include "mat_lzz_pX_determinant.h"

#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_linsolve.h"

//#define GENERIC_DET_PROFILE

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

//void determinant_via_diagonal_of_hermite(zz_pX & det, const Mat<zz_pX> & pmat)
//{
//    Vec<zz_pX> diag;
//    diagonal_of_hermite(diag, pmat);
//    set(det); // det = 1
//    for (long i = 0; i < diag.length(); ++i)
//        det *= diag[i];
//}

void determinant_expansion_by_minors(zz_pX & det, const Mat<zz_pX> & pmat)
{
    const long dim = pmat.NumRows();
    if (dim==1)
    {
        det = pmat[0][0];
        return;
    }

    if (dim==2)
    {
        zz_pX buf;
        mul(det, pmat[0][0], pmat[1][1]);
        mul(buf, pmat[0][1], pmat[1][0]);
        sub(det, det, buf);
        return;
    }

    if (dim==3)
    {
        // buffer for multiplications and for 2x2 determinants
        zz_pX buf, det22;

        // initialize det as pmat[2][0] * det(pmat[0,1][1,2])
        mul(det22, pmat[0][1], pmat[1][2]);
        mul(buf, pmat[0][2], pmat[1][1]);
        sub(det22, det22, buf);
        mul(det, pmat[2][0], det22);

        // det += pmat[2][1] * (-det(pmat[0,1][0,2]))
        mul(det22, pmat[0][2], pmat[1][0]);
        mul(buf, pmat[0][0], pmat[1][2]);
        sub(det22, det22, buf);
        mul(buf, pmat[2][1], det22);
        add(det, det, buf);

        // det += pmat[2][2] * det(pmat[0,1][0,1])
        mul(det22, pmat[0][0], pmat[1][1]);
        mul(buf, pmat[0][1], pmat[1][0]);
        sub(det22, det22, buf);
        mul(buf, pmat[2][2], det22);
        add(det, det, buf);

        return;
    }

    if (dim==4)
    {
        // buffer for multiplications and for 3x3 determinants
        zz_pX buf, det33;

        // store the six 2x2 determinants
        Vec<zz_pX> det22(INIT_SIZE, 6);

        // det22[0] = det(pmat[0,1][0,1])
        mul(det22[0], pmat[0][0], pmat[1][1]);
        mul(buf, pmat[0][1], pmat[1][0]);
        sub(det22[0], det22[0], buf);

        // det22[1] = det(pmat[0,1][0,2])
        mul(det22[1], pmat[0][0], pmat[1][2]);
        mul(buf, pmat[0][2], pmat[1][0]);
        sub(det22[1], det22[1], buf);

        // det22[2] = det(pmat[0,1][0,3])
        mul(det22[2], pmat[0][0], pmat[1][3]);
        mul(buf, pmat[0][3], pmat[1][0]);
        sub(det22[2], det22[2], buf);

        // det22[3] = det(pmat[0,1][1,2])
        mul(det22[3], pmat[0][1], pmat[1][2]);
        mul(buf, pmat[0][2], pmat[1][1]);
        sub(det22[3], det22[3], buf);

        // det22[4] = det(pmat[0,1][1,3])
        mul(det22[4], pmat[0][1], pmat[1][3]);
        mul(buf, pmat[0][3], pmat[1][1]);
        sub(det22[4], det22[4], buf);

        // det22[5] = det(pmat[0,1][2,3])
        mul(det22[5], pmat[0][2], pmat[1][3]);
        mul(buf, pmat[0][3], pmat[1][2]);
        sub(det22[5], det22[5], buf);

        // deduce the determinant
        // initialize det as pmat[3][0] * det(pmat[0,1,2][1,2,3])
        mul(det33, pmat[2][1], det22[5]);
        mul(buf, pmat[2][2], det22[4]);
        sub(det33, det33, buf);
        mul(buf, pmat[2][3], det22[3]);
        add(det33, det33, buf);

        mul(det, pmat[3][0], det33);
        NTL::negate(det, det);

        // det +=  pmat[3][1] * det(pmat[0,1,2][0,2,3])
        mul(det33, pmat[2][0], det22[5]);
        mul(buf, pmat[2][2], det22[2]);
        sub(det33, det33, buf);
        mul(buf, pmat[2][3], det22[1]);
        add(det33, det33, buf);

        mul(buf, pmat[3][1], det33);
        add(det, det, buf);

        // det -=  pmat[3][2] * det(pmat[0,1,2][0,1,3])
        mul(det33, pmat[2][0], det22[4]);
        mul(buf, pmat[2][1], det22[2]);
        sub(det33, det33, buf);
        mul(buf, pmat[2][3], det22[0]);
        add(det33, det33, buf);

        mul(buf, pmat[3][2], det33);
        NTL::negate(buf, buf);
        add(det, det, buf);

        // det +=  pmat[3][3] * det(pmat[0,1,2][1,2,3])
        mul(det33, pmat[2][0], det22[3]);
        mul(buf, pmat[2][1], det22[1]);
        sub(det33, det33, buf);
        mul(buf, pmat[2][2], det22[0]);
        add(det33, det33, buf);

        mul(buf, pmat[3][3], det33);
        add(det, det, buf);

        return;
    }

    std::cout << "naive det not implemented for this dimension" << std::endl;
}

void determinant_expansion_by_minors_rec(zz_pX & det, const Mat<zz_pX> & pmat)
{
    const long dim = pmat.NumRows();
    if (dim==1)
    {
        det = pmat[0][0];
        return;
    }

    // dim >= 2
    Mat<zz_pX> buf(INIT_SIZE, dim-1, dim-1);
    zz_pX tmp;
    clear(det);
    for (long k = 0; k < dim; ++k)
    {
        for (long i = 0; i < k; ++i)
            for (long j = 0; j < dim-1; ++j)
                buf[i][j] = pmat[i][j+1];
        for (long i = k+1; i < dim; ++i)
            for (long j = 0; j < dim-1; ++j)
                buf[i-1][j] = pmat[i][j+1];
        determinant_expansion_by_minors_rec(tmp, buf);
        mul(tmp, pmat[k][0], tmp);
        if (k%2 == 0) // k even
            add(det, det, tmp);
        else // k odd
            sub(det, det, tmp);
    }
}

bool determinant_generic_knowing_degree(zz_pX & det, const Mat<zz_pX> & pmat, long degree)
{
#ifdef GENERIC_DET_PROFILE
    double t;
    std::cout << "enter with dim = " << pmat.NumCols() << ", deg = " << deg(pmat) << std::endl;
#endif // GENERIC_DET_PROFILE
    
    const long dim = pmat.NumRows();
    if (dim<=4)
    {
#ifdef GENERIC_DET_PROFILE
            t=GetWallTime();
#endif // GENERIC_DET_PROFILE
        determinant_expansion_by_minors(det, pmat);
#ifdef GENERIC_DET_PROFILE
        std::cout << "\tbase case --> " << GetWallTime()-t << std::endl;
#endif // GENERIC_DET_PROFILE
        return (degree==deg(det));
    }

    const long cdim1 = (dim>>1);  // cdim1 ~ dim/2
    const long cdim2 = dim-cdim1;  // cdim2 ~ dim/2, cdim1+cdim2 = dim

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
    Mat<zz_pX> appbas;
    // degree of kernel basis will be (generically)  D = cdim1 * deg(pmat_l) / (dim - cdim1)
    // --> compute approximants at order deg(pmat_l) + D + 1
    // (cf for example Neiger-Rosenkilde-Solomatov ISSAC 2018, Lemma 4.3)
    long deg_pmat_l = deg(pmat_l);
    long deg_ker = ceil( cdim1 * deg_pmat_l / (double)(dim-cdim1) );
    long order = deg_pmat_l + deg_ker + 1;

    VecLong shift(dim,0);
#ifdef GENERIC_DET_PROFILE
    t = GetWallTime();
#endif // GENERIC_DET_PROFILE
    pmbasis(appbas, pmat_l, order, shift);
#ifdef GENERIC_DET_PROFILE
    t = GetWallTime()-t;
    std::cout << "\tpmbasis order = " << order << " || time " << t << std::endl;
#endif // GENERIC_DET_PROFILE

    // minimal left kernel basis of pmat_r : last rows of app
    Mat<zz_pX> kerbas;
    kerbas.SetDims(cdim2,dim);
    for (long i = 0; i < cdim2; ++i)
        for (long j = 0; j < dim; ++j)
            kerbas[i][j] = appbas[i+cdim1][j];

    // then compute the product
    Mat<zz_pX> pmatt;
#ifdef GENERIC_DET_PROFILE
    t = GetWallTime();
#endif // GENERIC_DET_PROFILE
    multiply(pmatt, kerbas, pmat_r);
#ifdef GENERIC_DET_PROFILE
    t = GetWallTime()-t;
    std::cout << "\tmultiply degrees " << deg(kerbas) << "," << deg(pmat_r) << " || time " << t << std::endl;
#endif // GENERIC_DET_PROFILE

    return determinant_generic_knowing_degree(det,pmatt,degree);
}

// Determinant via linear system solving with random right-hand side
void determinant_via_linsolve(zz_pX & det, const Mat<zz_pX> & pmat)
{
    Vec<zz_pX> u, v;
    random(v, pmat.NumRows(), max(1,deg(pmat)));
    linsolve_via_series(u, det, pmat, v);
}

void determinant_via_evaluation_general(zz_pX & det, const Mat<zz_pX> & pmat)
{
    const long nb_points = pmat.NumRows() * deg(pmat) + 1;
    zz_pX_Multipoint_General ev = get_general_points(nb_points);
    Vec<Mat<zz_p>> evals;
    ev.evaluate_matrix(evals, pmat);
    Vec<zz_p> det_evals(INIT_SIZE, nb_points);
    for (long k = 0; k < nb_points; ++k)
        determinant(det_evals[k], evals[k]);
    ev.interpolate(det, det_evals);
}

void determinant_via_evaluation_geometric(zz_pX & det, const Mat<zz_pX> & pmat)
{
    const long d = deg(pmat);
    const long nb_points = pmat.NumRows() * d + 1;
    zz_pX_Multipoint_Geometric ev = get_geometric_points(nb_points);
    ev.prepare_degree(d);
    Vec<Mat<zz_p>> evals;
    ev.evaluate_matrix(evals, pmat);
    Vec<zz_p> det_evals(INIT_SIZE, nb_points);
    for (long k = 0; k < nb_points; ++k)
        determinant(det_evals[k], evals[k]);
    ev.interpolate(det, det_evals);
}

void determinant_via_evaluation_FFT(zz_pX & det, const Mat<zz_pX> & pmat)
{
    // degree and dimension
    const long d = deg(pmat);
    const long dim = pmat.NumRows();

    // prepare FFT representation
    const long idxk = NextPowerOfTwo(pmat.NumRows() * d + 1);
    const long nb_points = 1 << idxk;
    fftRep R(INIT_SIZE, idxk);

    // matrix of evaluations of pmat: evals[k] contains
    // the evaluation of a at the k-th point
    Vec<Mat<zz_p>> evals(INIT_SIZE, nb_points);
    for (long k = 0; k < nb_points; ++k)
        evals[k].SetDims(dim,dim);

    for (long i = 0; i < dim; ++i)
        for (long j = 0; j < dim; ++j)
        {
            TofftRep(R, pmat[i][j], idxk);
            for (long k = 0; k < nb_points; ++k)
                evals[k][i][j].LoopHole() = R.tbl[0][k];
        }

    zz_p detk;
    for (long k = 0; k < nb_points; ++k)
    {
        determinant(detk, evals[k]);
        R.tbl[0][k] = detk._zz_p__rep;
    }
    FromfftRep(det, R, 0, pmat.NumRows() * d);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
