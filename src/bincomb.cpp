
//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"

using namespace Rcpp;

IntegerMatrix bincomb(const int p){
  int s_2=1,np=1<<p,n_2=np>>1;
  IntegerMatrix x(np,p);
  IntegerVector ones(1,1),zeros(1);
  for(int i=0;i<p;++i,n_2=n_2>>1,s_2=s_2<<1){
    x.column(i)=rep(combine(rep(zeros,n_2),rep(ones,n_2)),s_2);
  }
  return x;
}

RcppExport SEXP Rfast_bincomb(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type x(xSEXP);
    __result = bincomb(x);
    return __result;
END_RCPP
}