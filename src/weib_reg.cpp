//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib.h"
#include "Rfast.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List weib_reg(NumericVector Y, NumericMatrix X, const double tol = 1e-07, const int maxiters = 100){
  List l;
  int n = Y.size(), d = X.ncol();
  mat x(X.begin(),n,d,false);
  vec y(Y.begin(),n,false), be(d,fill::zeros),com(n),ini,lam;

  ini = weibull_mle2(y, n, tol, maxiters);

  double sly = sum_with<log,vec>(y), ek = ini[0], lik1 = ini[2], lik2;
  be[0] = log(ini[1]);
  double yhat0 = exp(-be[0]);
  my_pow2(y*yhat0,&com[0],ek,n);
  rowvec sx = sum(x);

  vec logcom = log(com);
  vec comlogcom = com%logcom;

  double derk = n + ek * sly + ek*n*(-be[0]) - sum(comlogcom);
  double derk2 = derk - n - sum(comlogcom%logcom), k;
  k = log(ek) - derk/derk2;


  mat xcom = x.each_col()%com;
  rowvec derb =  sum(xcom)-sx;

  mat derb2 = - ek * cross_x_y<mat,mat,vec>(xcom, x);
  be = be -solve(derb2, derb.t());

  lam = -(x*be);


  vec yhat = exp(lam);
  ek = exp(k);
  my_pow2(y%yhat,&com[0],ek,n);

  lik2 = n*k+(ek-1)*sly+ek*sum(lam)-sum(com);
  int i = 2;

  while (i++<maxiters && lik2-lik1 > tol ) {
    lik1 = lik2;


    logcom = log(com);
    comlogcom = com%logcom;
    derk = n + ek * (sly + sum(lam)) - sum(comlogcom);
    derk2 = derk - n - sum(comlogcom%logcom);
    xcom = x.each_col()%com;
    derb =  sum(xcom)- sx;
    derb2 = -ek * cross_x_y<mat,mat,vec>(xcom, x);
    k = k - derk/derk2;
    be = be - solve(derb2, derb.t());
    lam = -(x*be);
    yhat = exp(lam);
    ek = exp(k);
    my_pow2(y%yhat,&com[0],ek,n);
    lik2 = n*k+(ek-1)*sly+ek*sum(lam)-sum(com);
  }

  l["iters"] = i-1;
  l["loglik"] = n * k + (ek - 1) * sly + ek * sum(lam) - sum(com);
  l["shape"] = ek;
  l["be"] = be;

  return l;
}

RcppExport SEXP Rfast_weib_reg(SEXP YSEXP,SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(weib_reg(Y,X,tol,maxiters));
  return __result;
  END_RCPP
}
