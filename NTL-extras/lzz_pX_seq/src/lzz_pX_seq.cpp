#include "lzz_pX_seq.h"
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <queue>
using namespace std;

const bool verbose = false;

Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pX &a){
  Vec<zz_pX> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = S[i]*a;

  return res;
}

Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pXY &f){
  Vec<zz_pX> res;
  long dy = f.degY();
  res.SetLength(min(0, S.length() - dy));
  for (long i = 0; i < res.length(); i++){
    zz_pX val{0};
    for (long t = 0; t <= dy; t++)
      val += f.rep[t] * S[i+t];
    res[i] = val;
  }

  return res;
}

Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pX &a, const long d){
  Vec<zz_pX> res;
  res.SetLength(S.length());

  for (long i = 0; i < S.length(); i++)
    res[i] = MulTrunc(S[i],a,d);

  return res;

}

Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pXY &f, const long d){
  Vec<zz_pX> res;
  long dy = f.degY();
  res.SetLength(min(0, S.length() - dy));
  for (long i = 0; i < res.length(); i++){
    zz_pX val{0};
    for (long t = 0; t <= dy; t++)
      MulTrunc(val, f.rep[t], S[i+t], d);
    res[i] = val;
  }

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
  if (S.length() == 0) return true;
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

bool check_cancel(const Vec<zz_pX> &S, const Vec<zz_pXY> &gens,
    const long d){
  for (auto &f : gens){
    Vec<zz_pX> res;
    if (d == -1)
      res = mul(S,f);
    else
      res = mulTrunc(S,f,d);
    if (!is_zero_seq(res)) return false;
  }
  return true;
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
    for (long t = 0; t < d; t++){
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

void modified_kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens){
  gens.SetLength(0);
  map<long, Vec<zz_pX>> us; // keeps track of the current sequences
  map<long, zz_pXY> fs; // keeps track of the current poly

  // ideal generated by non-zero elements at index k from previous calcs
  map<long, Vec<zz_pX>> I_seq;
  // polynomial that generated the sequence with non-zero at index k
  map<long, zz_pXY> I_fs;

  vector<bool> bitmap; // bitmap[i] = true if x^t is potentially useful
  bitmap.resize(d,false);

  queue<long> Q;
  Q.push(0);
  bitmap[0] = true;
  while(!Q.empty()){
    long t = Q.front();
    Q.pop();
    zz_pX f_x;
    SetCoeff(f_x,t,1);
    zz_pXY f_xy{f_x};
    auto u = mulTrunc(S, f_x, d);

    if(verbose){
      cout << "t: " << t << endl;
      cout << "u: " << u << endl;
      cout << "f_xy: " << f_xy << endl;
    }

    // add to what we are tracking
    us[t] = u;
    fs[t] = f_xy;

    // update the ideals (checking isn't really needed)
    auto k = index_non_zero(u);
    auto it = I_seq.find(k);
    if (it == I_seq.end() && k >= 0){
      I_seq[k] = u;
      I_fs[k] = f_xy;
    }

    // check if we can find t' st (t'+t) < d and x^t' u[k] = 0
    long tprime = d - valuation(u[k]);
    if (tprime+t < d){ 
      // this means that x^(t'+t) is non-zero and x^(t'+t) *S 
      // will create more zero than prev
      bitmap[tprime+t] = true;
      Q.push(tprime+t);
    }
  }

  for (long s = 1; s < S.length(); s++){
    for (long b = 0; b < d; b++) 
      if (bitmap[b]) Q.push(b);
    while (!Q.empty()){
      auto t = Q.front();
      Q.pop();
      auto &f = fs[t];
      auto &u = us[t];

      auto fc = fs[t];
      auto uc = us[t];

      if(verbose){
	cout << "start at (s,t): (" << s << ", " << t << ")" << endl;
	cout << "u: " << u << endl;
	cout << "f: " << f << endl;
      }

      if (is_zero_seq(u)) continue;
      // shift by y
      f = shift_y(f,1);
      shift(u,1);

      // start the subiterations
      while (true){
	if (u.length() == 0) break;
	long k = index_non_zero(u);
	if (k == u.length()) break;
	long tprime = d - valuation(u[k]);
	if (tprime + t < d && bitmap[tprime+t] == false){
	  bitmap[tprime+t] = true;
	  Q.push(tprime+t);
	  fs[tprime+t] = shift_x(fc, tprime);
	  zz_pX x_tprime;
	  SetCoeff(x_tprime,tprime,1);
	  us[t+tprime] = mulTrunc(uc, x_tprime, d);
	}
	auto it = I_seq.find(k);
	if (it != I_seq.end()){
	  // means there is a prev sequence that starts
	  // with same number of non-zeros
	  auto &l = u[k];
	  auto &r = I_seq[k][k];
	  if (valuation(l) < valuation(r)){
	    // see if there is t' such that x^t' l has
	    // valuation >= r
	    tprime = valuation(r) - valuation(l);
	    if (tprime + t < d && bitmap[tprime+t] == false){
	      bitmap[tprime+t] = true;
	      Q.push(tprime+t);
	      fs[tprime+t] = shift_x(fc, tprime);
	      zz_pX x_tprime;
	      SetCoeff(x_tprime,tprime,1);
	      us[t+tprime] = mulTrunc(uc, x_tprime, d);
	    }
	    break;
	  }else{
	    zz_pX a{1}; 
	    solve(a,l,r,d);
	    if (verbose){
	      cout << "a: " << a << endl;
	    }
	    auto temp = mulTrunc(I_seq[k],-a,d);

	    if (verbose){
	      cout << "temp: " << temp << endl;
	      cout << "u before: " << u << endl;
	    }

	    u = add(u,temp);

	    f = trunc_x(f-a*I_fs[k],d);

	    if (verbose){
	      cout << "u after: " << u << endl;
	      cout << "f after: " << f << endl;
	    }
	  }
	}else{
	  // we already checked if there is t' such that
	  // x^t' u_t[k] = 0 earlier, so we can just break
	  break;
	}
      }
    }
    for (long t = 0; t < d; t++){
      if (bitmap[t] == false) continue;
      auto &u = us[t];
      auto &f = fs[t];
      auto k = index_non_zero(u);
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
  gens.SetLength(d);
  for (long t = 0; t < d; t++){
    if (bitmap[t] == true) gens[t] = fs[t];
  }
}

void fill_in(Vec<zz_pXY> &gens){
  zz_pXY zero{zz_pX{0}};
  for (long t = 1; t < gens.length(); t++){
    if (gens[t] == zero) gens[t] = shift_x(gens[t-1],1);
  }
}








