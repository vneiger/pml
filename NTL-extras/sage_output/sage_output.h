#ifndef SAGE_OUTPUT__H
#define SAGE_OUTPUT__H

#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic functions for output in sage format                  */
/* handles mainly Vec<zz_p>, zz_pX, ZZX                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void sage_init();

/*------------------------------------------------------------*/
/* initializes U<x>=GF(p)[x]                                  */
/*------------------------------------------------------------*/
void sage_init_X();

/*------------------------------------------------------------*/
/* initializes UZ<XX>=Z[XX]                                   */
/*------------------------------------------------------------*/
void sage_init_ZZX();


/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void sage_output(const Vec<ZZ> & v);
void sage_output(const Vec<long> & v);
void sage_output(const Vec<unsigned long> & v);
void sage_output(const Vec<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void sage_assign(const Vec<ZZ> & v, const string & name);
void sage_assign(const Vec<long> & v, const string & name);
void sage_assign(const Vec<unsigned long> & v, const string & name);
void sage_assign(const Vec<zz_p> & v, const string & name);


/*------------------------------------------------------------*/
/* prints a matrix                                            */
/*------------------------------------------------------------*/
void sage_output(const Mat<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a matrix to variable "name"                        */
/*------------------------------------------------------------*/
void sage_assign(const Mat<zz_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a poly (or vector, matrix) with indeterminate "var" */
/*------------------------------------------------------------*/
void sage_output(const ZZX & v, const string & var);
void sage_output(const zz_pX & v, const string & var);
void sage_output(const Vec<zz_pX> & v, const string & var);
void sage_output(const Mat<zz_pX>& a, const string & var);

/*------------------------------------------------------------*/
/* prints a poly (or vector) with indeterminate x             */
/*------------------------------------------------------------*/
void sage_output(const ZZX & v);
void sage_output(const zz_pX & v);
void sage_output(const Vec<zz_pX> & v);
void sage_output(const Mat<zz_pX>& a);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate "var" to variable "name"  */
/*------------------------------------------------------------*/
void sage_assign(const ZZX & v, const string & var, const string & name);
void sage_assign(const zz_pX & v, const string & var, const string & name);
void sage_assign(const Vec<zz_pX> & v, const string & var, const string & name);
void sage_assign(const Mat<zz_pX>& a, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void sage_assign(const ZZX & v, const string & name);
void sage_assign(const zz_pX & v, const string & name);
void sage_assign(const Vec<zz_pX> & v, const string & name);
void sage_assign(const Mat<zz_pX>& a, const string & var);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
