//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <Rinternals.h>

using namespace Rcpp;

NumericMatrix frame_to_matrix(DataFrame x){
  int i=0,p=x.length(),n=x.nrows();
  NumericMatrix f(n,p);
  NumericVector a;
  DataFrame::iterator xx;
  for(xx=x.begin();xx!=x.end();++xx,++i){
    a=*xx;
    f(_,i)=a;
  }
  return f;
}

RcppExport SEXP Rfast_frame_to_matrix(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< DataFrame >::type x(xSEXP);
    __result = wrap(frame_to_matrix(x));
    return __result;
END_RCPP
}