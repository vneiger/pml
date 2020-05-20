#include "lzz_pX_seq.h"
#include <map>
#include <iostream>
#include <cmath>
using namespace std;

const bool verbose = true;

Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pX &a){
  Vec<zz_pX> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = S[i]*a;

  return res;
}

Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pX &a, const long d){
  Vec<zz_pX> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = MulTrunc(S[i],a,d);

  return res;

}

Vec<zz_pX> add(const Vec<zz_pX> &S, const Vec<zz_pX> &T){
  Vec<zz_pX> res;
  long n = min(S.length(), T.length());
  res.SetLength(n);
  for (long i = 0; i < n; i++)
    res[i] = S[i] + T[i];
  return res;
}

long index_non_zero (const Vec<zz_pX> &S){
  if (S.length() == 0) return -1;
  for (long i = 0; i < S.length(); i++)
    if (S[i] != zz_pX(0)) return i;
  return S.length();
}

bool is_zero_seq(const Vec<zz_pX> &S){
  for (long i = 0; i < S.length(); i++)
    if (S[i] != zz_pX(0)) return false;
  return true;
}

long valuation(const zz_pX &a){
  if (deg(a) == -1) return -1;
  for (long i = 0; i <= deg(a); i++){
    if (a[i] != zz_p(0)) return i;
  }
  return 0; // should never happen
}

void shift(Vec<zz_pX> &S, const long s){
  if (s >= S.length()){
    S.SetLength(0);
    return;
  }
  for (long i = 0; i + s < S.length(); i++)
    S[i] = S[i+s];
  S.SetLength(S.length()-s);

}

// solves for a in l = a*r mod x^d
void solve(zz_pX &a, const zz_pX &l, const zz_pX &r, const long d){
  long s = valuation(r);
  auto l_shifted = RightShift(l,s);
  auto r_shifted = RightShift(r,s);

  zz_pX mod;
  SetCoeff(mod,d,1);

  MulMod(a, l_shifted, InvMod(r_shifted, mod), mod);
  cout <<"a at solve: " << a << endl; 
}

void kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens){
  Vec<Vec<zz_pX>> us;
  us.SetLength(d);
  gens.SetLength(d);

  map<long, Vec<zz_pX>> I_seq;
  map<long, zz_pXY> I_fs;

  // phase 1
  zz_pX running;
  SetCoeff(running, 0, 1);
  for (long t = 0; t < d; t++){
    gens[t] = zz_pXY(running);
    us[t] = mulTrunc(S, running, d);
    LeftShift(running, running, 1);

    long k = index_non_zero(us[t]);
    auto it = I_seq.find(k);
    if (it == I_seq.end() && k >= 0){
      I_seq[k] = us[t];
      I_fs[k] = gens[t];
    } 

    if (verbose){
      cout << "at: " << t << endl;
      cout << "f_t: " << gens[t] << endl;
      cout << "u_s[t]" << endl;
      for (long i = 0; i < us[t].length(); i++)
	cout << us[t][i] << endl;
      cout << endl;
    }

    // we don't need to worry about the case where there
    // is already an entry at k, since the previous entry
    // will always have lower valuation
  }

  // phase 2
  for (long s = 1; s < S.length(); s++){
    for (long t = 0; t < d; t++){
      auto &f = gens[t];
      auto &u = us[t];
	
      if (u.length() == 0) continue;
      if (is_zero_seq(u)) continue;

      shift(u,1);
      f = shift_y(f,1);

      while (true){
	if (u.length() == 0) break;
	if (is_zero_seq(u)) break;

	cout << "at (s,t): " << s << ", " << t << endl;
	cout << "f: " << f << endl;
	cout << "u:" << endl;
	for (long i = 0; i < us[t].length(); i++)
	  cout << us[t][i] << endl;

	long k = index_non_zero(u);
	if (k == u.length()) break;

	auto it = I_seq.find(k);
	if (it != I_seq.end() && k >= 0){
	  // means there exists a sequece that
	  // appeared before with the same
	  // number of leading zeros; try
	  // using this to cancel out
	  //the leading non-zero term
	  auto &l = u[k];
	  auto &r = I_seq[k][k];
	  if (valuation(l) < valuation(r))
	    break;
	  else{
	    zz_pX a{1};
	    solve(a,l,r,d);
	    cout << "a: " << a << endl;

	    // update u
	    auto temp = mulTrunc(I_seq[k], -a, d);
	    u = add(u, temp);
	    cout << "u after:" << endl;
	    for (auto & i : u)
	      cout << i << endl;

	    // update f
	    f = trunc_x(f - a * I_fs[k], d);
	    cout << "f after: " << f << endl << endl;
	  }
	}else{
	  // means there exists no prev seq
	  // that has the same number of leading
	  // zeros
	  break;

	}
      }
    }
    // update the ideals
    for (long t = 0; t < d-1; t++){
      auto &u = us[t];
      auto &f = gens[t];
      long k = index_non_zero(u);
      if (k < u.length() && k >= 0){
	auto it = I_seq.find(k);
	if (it != I_seq.end()){
	  if (valuation(u[k]) < valuation(I_seq[k][k])){
	    I_seq[k] = u;
	    I_fs[k] = f;
	  }
	}else{
	  I_seq[k] = u;
	  I_fs[k] = f;
	}
      }
    }

  }
}












