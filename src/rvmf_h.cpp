
//Author: Manos Papadakis

#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector rvmf_h(double n,double ca,double d1,double x0,double m,double k,double b){
  double ta,u,z,tmp=0;
  NumericVector w(n);
  for(int i=0;i<n;++i) {
    for(ta=-1000.0,u=1.0;ta - ca < log(u);) {
      z = R::rbeta(m,m);
      u = R::runif(0,1);
      tmp = ( 1 - (1 + b) * z ) / ( 1 - (1 - b) * z );
      ta = k * tmp + d1 * log(1 - x0 * tmp);
    }
    w[i]=tmp;
  }
  return w;
}

RcppExport SEXP Rfast_rvmf_h(SEXP nSEXP,SEXP caSEXP,SEXP d1SEXP,SEXP x0SEXP,SEXP mSEXP,SEXP kSEXP,SEXP bSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< double >::type n(nSEXP);
    traits::input_parameter< double >::type ca(caSEXP);
    traits::input_parameter< double >::type d1(d1SEXP);
    traits::input_parameter< double >::type x0(x0SEXP);
    traits::input_parameter< double >::type m(mSEXP);
    traits::input_parameter< double >::type k(kSEXP);
    traits::input_parameter< double >::type b(bSEXP);
    __result = rvmf_h(n,ca,d1,x0,m,k,b);
    return __result;
END_RCPP
}