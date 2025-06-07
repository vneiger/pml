/* 
   Copyright (C) 2024 Gilles Villard

   This file is part of mapml. mapml is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   mapml is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with mapml. If not, see <http://www.gnu.org/licenses/>. */

// *******************************************************


#ifndef MAPML_MATPOLY_EXPORT_C
#define MAPML_MATPOLY_EXPORT_C

#include <string.h>
#include "mapml.h"
#include "conversion.h"

/***********************************************************
 * 
 * 
 *  Reccurence of order m  (m+1 coefficients) 
 * 
 * 
 * 
 * *********************************************************/


ALGEB pm_coeffs(MKernelVector kv, ALGEB *args){

    slong i,k;

    char *str;  // tmp for trace 
    
    // The order of the recurrence
    M_INT m=MapleToInteger32(kv,args[1]);

    // The polynomial coefficients of teh recurrence 
    ALGEB vpol=args[2];

    fmpq_poly_t vp[m+1];  // Gives the polynomial recurrence 
    get_fmpq_poly_array(vp, kv, vpol); 

    // N+1 terms requested 
    M_INT N;
    N=MapleToInteger64(kv,args[3]);

    MapleALGEB_Printf(kv, " \n");
    MapleALGEB_Printf(kv, " ---- %a terms \n",ToMapleInteger(kv,N));

    // The initial conditions 
    //-----------------------

    fmpq_t cinit[m];  
    for (i=0; i<m; i++){

        fmpq_init(cinit[i]);

        fmpq_set_str(cinit[i], MapleToString(kv,args[4+i]), 10);

        //str = fmpq_get_str(NULL, 10, cinit[i]);
        //MapleALGEB_Printf(kv, str);
        //MapleALGEB_Printf(kv, "  \n");
        //MapleALGEB_Printf(kv, "  \n");

    }

    
    // The resulting rationals, terms of the sequence 
    //-----------------------------------------------

    fmpq_t L[N+1];    


    for (i=0; i<m; i++){

        fmpq_init(L[i]);
        fmpq_set(L[i],cinit[i]);

    }


    fmpq_t evalp[m+1];  // the evaluation of the polynomials 

    // Main loop 
    //----------

    fmpq_t res,tt;
    fmpq_init(res);
    fmpq_init(tt);

    fmpz_t e;
    fmpz_init(e);

    for (k=m; k<N+1; k++){

        fmpq_init(L[k]); // Set to 0 

        // Evaluations
        for (i=0; i<m+1; i++){

            fmpz_init_set_si(e,k-m);

            fmpq_poly_evaluate_fmpz(evalp[i], vp[i], e);
            
        }
        
        fmpq_mul(res, L[k-m], evalp[0]);

        for (i=1; i<m; i++){

            fmpq_addmul(res,L[k-m+i],evalp[i]);
            
        }

        fmpq_neg(tt,res);
        fmpq_div(L[k],tt,evalp[m]);

    } // End main loop k 



    ALGEB listcoeff;

    listcoeff= MapleListAlloc(kv, N+1);

    for (k = 1; k <= N+1; k++) {

        str = fmpq_get_str(NULL, 10, L[k-1]);

        MapleListAssign(kv,listcoeff,k, ToMapleString(kv,str));
    }

    return listcoeff;



    // //+++++++++++++++++++++++++++++++++


    // // Dummy print

    // MapleALGEB_Printf(kv, " \n");
    // MapleALGEB_Printf(kv, " ---- \n");


    // for (i=0; i<N+1; i++){
    
    //     str = fmpq_get_str(NULL, 10, L[i]);
    //     MapleALGEB_Printf(kv, str);
    //     MapleALGEB_Printf(kv, "  \n");
    //     MapleALGEB_Printf(kv, " \n");

    // }
    // //MapleALGEB_Printf(kv, " \n");
    // //MapleALGEB_Printf(kv, " ++++  \n");

    // return ToMapleString(kv,str);

}


/**********************************************************
 * 
 * modulo matrix polynomial determinant  
 * 
 *  ALGEB args[1]: matrix polynomial, vector entries 
 *        args[2]: modulus 
 * 
 * Returns the determinant as a list of coefficients 
 * 
 ***********************************************************/


ALGEB pm_determinant(MKernelVector kv, ALGEB *args){

    ALGEB vectmat=args[1];

    ulong modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_mat_t A;

    get_nmod_poly_mat(A, modulus, kv, vectmat);

    nmod_poly_t p;  

    nmod_poly_init(p, modulus); // Inittialization probably not done by flint below ? 

    nmod_poly_mat_det(p, A);
    
    return nmod_poly_to_algeb(kv,p);

}


/*******************************************************************
 * 
 * Approximant row basis for a vector (m x 1 matrix) of derivatives 
 * 
 *  ALGEB args[1]: shift
 *        args[2]: polynomial (vector of coefficients)  
 *        args[3]: order of approximation 
 *        args[4]: modulus 
 *        args[5]: method, "MBasis" or "PMBasis" 
 * 
 *  Returns ALGEB  res:=[M,dct] 
 *    M: a polynomial matrix, list entries, m x m 
 *    dct: list of out defects  
 * 
 *  Eventually, be careful with the sign between either 
 *      the defect (e.g. in gfun) or shifts = -dct in pml
 * 
 *******************************************************************/


ALGEB pm_diff_mbasis(MKernelVector kv, ALGEB *args){


    //double t = 0.0;
    //clock_t tt;

    ulong modulus = MapleToInteger64(kv,args[4]);

    nmod_poly_t p;    // The polynomial that will be differentiated 
    nmod_poly_init(p,modulus);
    get_nmod_poly(p, modulus, kv, args[2]);


    //The dimension is deduced from the size of the shift 
    M_INT m = MapleNumArgs(kv,args[1]);


    // the vector of derivatives 
    nmod_poly_mat_t Vdiff;
    nmod_poly_mat_init(Vdiff, m, 1, modulus);  // No use of vectors for the moment, matrices 

    nmod_poly_set(nmod_poly_mat_entry(Vdiff,0,0),p);

    for (ulong i = 1; i < m; i++)
     nmod_poly_derivative(nmod_poly_mat_entry(Vdiff,i,0),nmod_poly_mat_entry(Vdiff,i-1,0));


    // Resulting basis  
 nmod_poly_mat_t M;
 nmod_poly_mat_init(M, m, m, modulus); 

 ALGEB maple_shift = args[1];
 slong shift[m];
 for (ulong i = 0; i < m; i++) 
    shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

ulong order = MapleToInteger64(kv,args[3]);


    // appromimant computation using PML 

    //tt = clock();


if (strcmp(MapleToString(kv,args[5]),"MBasis") ==0) 
    nmod_poly_mat_mbasis(M, shift, Vdiff, order); 
else 
nmod_poly_mat_pmbasis(M, shift, Vdiff, order); 


   //t = (double)(clock()-tt) / CLOCKS_PER_SEC;


   //MapleALGEB_Printf(kv, MapleToString(kv,args[5]));
   //MapleALGEB_Printf(kv, " ++++  Time diffbasis %f ms\n", ToMapleFloat(kv,t*1000));


   // Construction of the result [M, dct]
ALGEB outdct;

outdct= MapleListAlloc(kv,m);

for (slong i = 1; i <= m; i++) 
    MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));

ALGEB res= MapleListAlloc(kv,2);
MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,M));
MapleListAssign(kv,res,2,outdct);

return res;

}


/**********************************************************
 * 
 * modulo polynomial matrix mbasis  
 * 
 *  ALGEB args[1]: shift
 *        args[2]: matrix polynomial, vector entries 
 *        args[3]: order
 *        args[4]: modulus 
 * 
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

ALGEB pm_matrix_mbasis(MKernelVector kv, ALGEB *args){


   //double t = 0.0;
   //clock_t tt;


   ALGEB vectmat=args[2];

   ulong modulus = MapleToInteger64(kv,args[4]);

   nmod_poly_mat_t A;

   get_nmod_poly_mat(A, modulus, kv, vectmat);

   slong m = A ->r;

   nmod_poly_mat_t M;

   nmod_poly_mat_init(M, m, m, modulus); 

   //slong res_shift[m];

   ALGEB maple_shift = args[1];

   slong shift[m];
   for (ulong i = 0; i < m; i++) 
      shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

  ulong order = MapleToInteger64(kv,args[3]);


  //tt = clock();

  nmod_poly_mat_mbasis(M, shift, A, order);
   //mbasis(M, res_shift, A, order, shift);

  //t = (double)(clock()-tt) / CLOCKS_PER_SEC;
  //MapleALGEB_Printf(kv, " ++++  Time mbasis %f ms\n", ToMapleFloat(kv,t*1000));


  ALGEB outdct;

  outdct= MapleListAlloc(kv,m);

  for (slong i = 1; i <= m; i++) 
    MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));


ALGEB res= MapleListAlloc(kv,2);
MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,M));
MapleListAssign(kv,res,2,outdct);

return res;

}


/**********************************************************
 * 
 * modulo matrix polynomial pmbasis  
 * 
 *  ALGEB args[1]: shift
 *        args[2]: matrix polynomial, vector entries 
 *        args[3]: order
 *        args[4]: modulus 
 * 
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

ALGEB pm_matrix_pmbasis(MKernelVector kv, ALGEB *args){


   //double t = 0.0;
   //clock_t tt;


   ALGEB vectmat=args[2];

   ulong modulus = MapleToInteger64(kv,args[4]);

   nmod_poly_mat_t A;

   get_nmod_poly_mat(A, modulus, kv, vectmat);

   slong m = A ->r;

   nmod_poly_mat_t M;

   nmod_poly_mat_init(M, m, m, modulus); 

   //slong res_shift[m];

   ALGEB maple_shift = args[1];

   slong shift[m];
   for (ulong i = 0; i < m; i++) 
      shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

  ulong order = MapleToInteger64(kv,args[3]);


  //tt = clock();

  nmod_poly_mat_pmbasis(M, shift, A, order);
   //mbasis(M, res_shift, A, order, shift);

  //t = (double)(clock()-tt) / CLOCKS_PER_SEC;
  //MapleALGEB_Printf(kv, " ++++  Time pmbasis %f ms\n", ToMapleFloat(kv,t*1000));


  ALGEB outdct;

  outdct= MapleListAlloc(kv,m);

  for (slong i = 1; i <= m; i++) 
    MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));


ALGEB res= MapleListAlloc(kv,2);
MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,M));
MapleListAssign(kv,res,2,outdct);

return res;

}

//ALGEB coeffs = args[1]; // Work with coeffs in input? 

//     M_INT l = MapleToInteger64(kv,args[2]);

//     MapleALGEB_Printf(kv, " Length %a", ToMapleInteger(kv,l));

//     for (M_INT i=0; i<l; i++) {

//         nmod_poly_set_coeff_ui(p, ), 
//            MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+2)));



/**********************************************************
 * 
 * modulo polynomial matrix weak Popov form   
 * 
 * ++++++++ TODO COMMENT 
 * 
 *  ALGEB args[1]: matrix polynomial, vector entries 
 *        args[2]: shift 
 *        args[3]: modulus 
 * 
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// SEE WHICH ENTRIES

ALGEB pm_weakpopov(MKernelVector kv, ALGEB *args){


   //double t = 0.0;
   //clock_t tt;


    ALGEB vectmat=args[1];

    ulong modulus = MapleToInteger64(kv,args[3]);

    nmod_poly_mat_t A;

    get_nmod_poly_mat(A, modulus, kv, vectmat);


    ALGEB maple_shift = args[2];

    slong m = A ->c;

    slong shift[m];

    for (ulong i = 0; i < m; i++) 
        shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));


    //MapleALGEB_Printf(kv, " ++++ %d\n ", ToMapleInteger(kv,shift[0]));

    slong n = A ->r;
    slong pivind[n];

    nmod_poly_mat_ordered_weak_popov_iter(A, shift, NULL, pivind, NULL, ROW_LOWER);


   // nmod_poly_mat_t M;

   // nmod_poly_mat_init(M, m, m, modulus); 

   // //slong res_shift[m];

   // ALGEB maple_shift = args[1];

   // ulong order = MapleToInteger64(kv,args[3]);

   // tt = clock();

   // nmod_poly_mat_pmbasis(M, shift, A, order);
   // //mbasis(M, res_shift, A, order, shift);

   // t = (double)(clock()-tt) / CLOCKS_PER_SEC;
   // MapleALGEB_Printf(kv, " ++++  Time pmbasis %f ms\n", ToMapleFloat(kv,t*1000));


   // ALGEB outdct;

   // outdct= MapleListAlloc(kv,m);

   // for (slong i = 1; i <= m; i++) 
   //      MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));


    ALGEB res= MapleListAlloc(kv,1);
    MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,A));

    return res;

}

//ALGEB coeffs = args[1]; // Work with coeffs in input? 

//     M_INT l = MapleToInteger64(kv,args[2]);

//     MapleALGEB_Printf(kv, " Length %a", ToMapleInteger(kv,l));

//     for (M_INT i=0; i<l; i++) {

//         nmod_poly_set_coeff_ui(p, ), 
//            MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+2)));



#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
