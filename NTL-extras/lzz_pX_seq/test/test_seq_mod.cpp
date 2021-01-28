#include "lzz_pX_seq.h"
#include <NTL/lzz_pX.h>
#include <iostream>
#include <sstream>
#include <NTL/vector.h>

NTL_CLIENT
using namespace std;

int main(){
  long p = 9001;
  zz_p::init(p);

  bool verbose = false;
  //bool verbose = true;

  long d;
  if (verbose)
    cout << "Enter d: ";
  cin >> d;

  long ell;
  if (verbose)
    cout << "Enter ell: ";
  cin >> ell;

  long tau;
  if (verbose)
    cout << "Enter tau: ";
  cin >> tau;


  cout << "d: " << d << endl;
  cout << "ell: " << ell << endl;
  cout << "tau: " << tau << endl;

  Vec<Vec<zz_pX>> S;
  for (long i = 0; i < ell; i++){
    if (verbose) cout << "At entry: " << i << endl;
    Vec<zz_pX> cur;
    for (long j = 0; j < tau; j++){
      zz_pX curpol;
      if (verbose)
      	cout << "subentry: " << j << endl;
      cin >> curpol;
      cur.append(curpol);
    }
    S.append(cur);
  }
  Vec<zz_pXY> gens;
  double time = GetWallTime();
  kurakin(d,S,gens);
  time = GetWallTime()-time;
  if(verbose){
    cout << "gens:" << endl;
    for (long i = 0; i < gens.length(); i++)
      cout << gens[i] << endl;
  }
  cout << "kurakin took " << time << endl;
  cout << "checking kurakin" << endl;
  if (check_cancel(S, gens, d))
    cout << "OKAY" << endl;
  else
    cout << "BAD" << endl;
  cout << endl;


  Vec<zz_pXY> gens2;
  time = GetWallTime();
  modified_kurakin(d,S,gens2);
  time = GetWallTime()-time;
  if(verbose){
    cout << "gens2: " << endl;
    for (long i = 0; i < gens2.length(); i++)
      cout << gens2[i] << endl;
  }
  long d_s = 0;
  for (long i = 0; i < gens2.length(); i++)
   if (gens2[i] != zz_pX{0}) d_s++;
  cout << "mo. kurakin took " << time << endl;
  cout << "d* = " << d_s << endl; 
  cout << "checking mo. kurakin" << endl;
  if (check_cancel(S, gens2, d))
    cout << "OKAY" << endl;
  else
    cout << "BAD" << endl;
 cout << endl;

  Vec<zz_pXY> gens3;
  time = GetWallTime();
  berlekamp_massey_pmbasis(d, S, gens3);
  time = GetWallTime() - time;
  if(verbose){
    cout << "gens3: " << endl;
    for (long i = 0; i < gens3.length(); i++)
      cout << gens3[i] << endl;
  }
  cout << "pmbasis took " << time << endl;
  cout << "checking pmbasis" << endl;
  if (check_cancel(S, gens3, d))
    cout << "OKAY" << endl;
  else
    cout << "BAD" << endl;


  Vec<zz_pXY> gens4;
  time = GetWallTime();
  berlekamp_massey_pmbasis_compressed(d, S, gens3);
  time = GetWallTime() - time;
  if(verbose){
    cout << "gens4: " << endl;
    for (long i = 0; i < gens4.length(); i++)
      cout << gens4[i] << endl;
  }
  cout << "pmbasis compressed took " << time << endl;
  cout << "checking pmbasis compressed" << endl;
  if (check_cancel(S, gens4, d))
    cout << "OKAY" << endl;
  else
    cout << "BAD" << endl;
}





