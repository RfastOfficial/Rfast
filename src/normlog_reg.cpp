//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List normlog_reg(NumericVector Y,NumericMatrix X, const double tol, const int maxiters){
  List l;
  int n = Y.size();

  int pcols = X.ncol();
  mat x(X.begin(),n,pcols,false),xt = x.t();
  colvec y(Y.begin(),n,false);
  colvec b1 = solve(xt * x, xt * log(y + 0.1));
  colvec yhat = (exp(x * b1));

  mat com = x.each_col()%(yhat % yhat);

  mat com2 = x.each_col()%(yhat % y);
  colvec der = conv_to<colvec>::from(sum(com-com2));

  mat der2 = xt*(2 * com - com2);

  colvec b2 = b1 - solve(der2, der);
  int i = 2;

  while(i++<maxiters && sum_with<abs,colvec>(b1-b2) > tol){
    b1 = b2;
    yhat = (exp(x * b1));
    com = x.each_col()%(yhat % yhat);
    com2 = x.each_col()%(yhat % y);
    der = conv_to<colvec>::from(sum(com-com2));
    der2 = xt*(2 * com - com2);
    b2 = b1 - solve(der2, der,solve_opts::fast);
  }

  double deviance = sum_with< square2<double>, vec>(y - yhat);
  double loglik =  - n * 0.5 * (log(deviance/n)+2.83787707);
  l["iters"] = i;
  l["loglik"] = loglik;
  l["deviance"] = deviance;
  l["be"] = b2;
  return l;
}

RcppExport SEXP Rfast_normlog_reg(SEXP YSEXP,SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(normlog_reg(Y,X,tol,maxiters));
  return __result;
  END_RCPP
}
