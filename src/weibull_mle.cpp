//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib.h"
#include "Rfast/templates.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List weibull_mle(NumericVector X, const double tol, const int maxiters){
  int n = X.size(),i=2;
  vec x(X.begin(),n,false);
  vec lx = foreach<std::log,vec>(x);
  vec lx2 = lx%lx;
  vec y = x;
  double mlx = sum(lx)/n, co = sum(y%lx),sy = sum(y), fb = 1+mlx-co/sy,fb2 = -1 -(sum(y%lx2)*sy-co*co)/(sy*sy);
  double b1 = 1,b2 = 1 - fb/fb2;

  while (++i<maxiters && sum(abs(b2 - b1)) > tol) {
    b1 = b2;
    my_pow2(x,y.memptr(),b1,n);
    co = sum(y % lx);
    sy = sum(y);
    fb = 1/b1 + mlx - co/sy;
    fb2 = -1/(b1*b1) - (sum(y % lx2) * sy - co*co)/(sy*sy);
    b2 = b1 - fb/fb2;
  }
  List l;

  l["iters"] = i-1;
  double theta = pow(sy/n,1/b2);
  my_pow2(conv_to<vec>::from(x/theta),y.memptr(),b2,n);
  l["loglik"] = n * log(b2) - n * b2 * log(theta) + (b2 - 1) * n * mlx - sum(y);
  vec param(2);
  param[0] = b2;
  param[1] = theta;

  l["param"] = conv_to<rowvec>::from(param);

  return l;
}

RcppExport SEXP Rfast_weibull_mle(SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type X(XSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = weibull_mle(X,tol,maxiters);
    return __result;
END_RCPP
}
