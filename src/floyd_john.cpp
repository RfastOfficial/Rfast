//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include "mn.h"

using namespace Rcpp;

NumericVector floyd_john(const int n,NumericVector x){
  NumericVector y=clone(x);
  i4mat_floyd(n,y);
  return y;
}

RcppExport SEXP Rfast_floyd_john(SEXP nSEXP,SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int >::type n(nSEXP);
    traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = floyd_john(n,x);
    return __result;
END_RCPP
}