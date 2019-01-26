#ifndef MAGMA_OUTPUT__H
#define MAGMA_OUTPUT__H

#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic functions for output in magma format                 */
/* handles mainly Vec<zz_p>, zz_pX, ZZX                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void magma_init();

/*------------------------------------------------------------*/
/* initializes U<x>=GF(p)[x]                                  */
/*------------------------------------------------------------*/
void magma_init_X();

/*------------------------------------------------------------*/
/* initializes UZ<XX>=Z[XX]                                   */
/*------------------------------------------------------------*/
void magma_init_ZZX();

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ> & v);
void magma_output(const Vec<long> & v);
void magma_output(const Vec<unsigned long> & v);
void magma_output(const Vec<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ> & v, const string & name);
void magma_assign(const Vec<long> & v, const string & name);
void magma_assign(const Vec<unsigned long> & v, const string & name);
void magma_assign(const Vec<zz_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a matrix                                            */
/*------------------------------------------------------------*/
void magma_output(const Mat<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a matrix to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Mat<zz_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a poly (or vector, matrix) with indeterminate "var" */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v, const string & var);
void magma_output(const zz_pX & v, const string & var);
void magma_output(const Vec<zz_pX> & v, const string & var);
void magma_output(const Mat<zz_pX>& a, const string & var);

/*------------------------------------------------------------*/
/* prints a poly (or vector) with indeterminate x             */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v);
void magma_output(const zz_pX & v);
void magma_output(const Vec<zz_pX> & v);
void magma_output(const Mat<zz_pX>& a);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate "var" to variable "name"  */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & var, const string & name);
void magma_assign(const zz_pX & v, const string & var, const string & name);
void magma_assign(const Vec<zz_pX> & v, const string & var, const string & name);
void magma_assign(const Mat<zz_pX>& a, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & name);
void magma_assign(const zz_pX & v, const string & name);
void magma_assign(const Vec<zz_pX> & v, const string & name);
void magma_assign(const Mat<zz_pX>& a, const string & var);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
